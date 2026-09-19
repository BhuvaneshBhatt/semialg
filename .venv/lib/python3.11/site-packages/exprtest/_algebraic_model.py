"""Exact algebraic recognition and rational-polynomial models."""

from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
from math import prod

import sympy as sp
from sympy.polys.numberfields import to_number_field

from . import _config as cfg
from ._cost import stage_allowed
from ._errors import EXACT_METHOD_ERRORS
from ._memo import minpoly_for
from ._result import Verdict, ZeroClassification


@dataclass(frozen=True)
class AlgebraicExpressionInfo:
    """Structural classification of a closed exact algebraic expression."""

    is_algebraic: bool
    is_simple: bool
    detail: str = ""


@dataclass(frozen=True)
class AlgebraicModel:
    """Polynomial model of an exact algebraic expression.

    ``expression`` is represented by ``numerator / denominator`` after each
    algebraic generator is replaced by the corresponding element of
    ``variables``.  ``minpolys[i]`` is the minimal polynomial over QQ of
    ``generators[i]``, written in ``variables[i]``.
    """

    expression: sp.Expr
    generators: tuple[sp.Expr, ...]
    variables: tuple[sp.Symbol, ...]
    minpolys: tuple[sp.Poly, ...]
    numerator: sp.Expr
    denominator: sp.Expr

    @property
    def degree_product(self) -> int:
        return prod(poly.degree() for poly in self.minpolys) if self.minpolys else 1


def _contains_inexact_number(expr: sp.Expr) -> bool:
    return any(isinstance(atom, (sp.Float,)) for atom in expr.atoms(sp.Number))


def _simple_algebraic_structure(expr: sp.Expr) -> bool:
    """Whether ``expr`` is built from algebraic atoms by cheap operations."""
    if expr.is_Rational:
        return True
    if expr in (sp.I, sp.GoldenRatio):
        return True
    if isinstance(expr, sp.AlgebraicNumber):
        return True
    if getattr(expr, "is_CRootOf", False):
        return True
    if expr.is_Add or expr.is_Mul:
        return all(_simple_algebraic_structure(arg) for arg in expr.args)
    if expr.is_Pow and len(expr.args) == 2:
        base, exponent = expr.args
        if not exponent.is_Rational:
            return False
        if exponent.is_negative and base.is_zero is True:
            return False
        return _simple_algebraic_structure(base)
    return bool(expr.is_number and expr.is_algebraic is True)


def classify_algebraic_expression(expr: sp.Expr) -> AlgebraicExpressionInfo:
    """Classify whether ``expr`` is a closed exact algebraic expression.

    This is intentionally stricter than asking SymPy whether an expression is
    numeric.  Machine/inexact values, infinities, NaN, free symbols, and
    expressions not known to be algebraic are rejected.
    """
    expr = sp.sympify(expr)
    if expr.free_symbols:
        return AlgebraicExpressionInfo(False, False, "expression has free symbols")
    if _contains_inexact_number(expr):
        return AlgebraicExpressionInfo(False, False, "expression contains inexact numbers")
    if expr.has(sp.zoo, sp.oo, -sp.oo, sp.nan):
        return AlgebraicExpressionInfo(False, False, "expression contains a non-finite value")
    try:
        algebraic = expr.is_algebraic is True
    except EXACT_METHOD_ERRORS:
        algebraic = False
    if not algebraic:
        return AlgebraicExpressionInfo(False, False, "expression is not known to be algebraic")
    simple = _simple_algebraic_structure(expr)
    return AlgebraicExpressionInfo(True, simple, "closed exact algebraic expression")


def _collect_generators(expr: sp.Expr, out: list[sp.Expr]) -> None:
    """Collect maximal algebraic atoms needed to rationalize ``expr``."""
    if expr.is_Rational:
        return
    if expr.is_Add or expr.is_Mul:
        for arg in expr.args:
            _collect_generators(arg, out)
        return
    if expr.is_Pow and expr.exp.is_Integer:
        _collect_generators(expr.base, out)
        return

    # Any remaining closed algebraic node can be treated as one generator.
    # Taking maximal nodes is useful for nested radicals: its minpoly already
    # captures the nested algebraic construction, avoiding needless generators.
    if not expr.free_symbols and expr.is_algebraic is True:
        if expr not in out:
            out.append(expr)
        return
    raise ValueError(f"cannot model non-algebraic subexpression: {expr}")


@lru_cache(maxsize=cfg.GENERATOR_CACHE_SIZE)
def extract_algebraic_generators(expr: sp.Expr) -> tuple[sp.Expr, ...]:
    """Return maximal exact algebraic generators occurring in ``expr``."""
    info = classify_algebraic_expression(expr)
    if not info.is_algebraic:
        return ()
    generators: list[sp.Expr] = []
    _collect_generators(sp.sympify(expr), generators)
    return tuple(generators)



def build_algebraic_model(expr: sp.Expr) -> AlgebraicModel | None:
    """Build a rational-polynomial model and generator minimal polynomials.

    Returns ``None`` when exact modeling is unsupported or exceeds the
    configured generator/minimal-polynomial complexity limits.
    """
    expr = sp.sympify(expr)
    if not stage_allowed(expr, "minpoly"):
        return None

    info = classify_algebraic_expression(expr)
    if not info.is_algebraic:
        return None
    if expr.is_Rational:
        return AlgebraicModel(expr, (), (), (), expr, sp.Integer(1))

    try:
        modeled_expr = sp.cancel(expr)
        generators = extract_algebraic_generators(modeled_expr)
        if len(generators) > cfg.ALGEBRAIC_MAX_GENERATORS:
            return None
        variables = tuple(sp.Dummy(f"a{i}") for i in range(len(generators)))
        replacement = dict(zip(generators, variables))
        represented = modeled_expr.xreplace(replacement)

        # After generator replacement the expression is rational in the fresh
        # variables.  cancel() puts it in a numerator/denominator form without
        # making any algebraic assumptions about the generators themselves.
        represented = sp.cancel(represented)
        numerator, denominator = sp.fraction(represented)
        sp.Poly(numerator, *variables, domain=sp.QQ)
        sp.Poly(denominator, *variables, domain=sp.QQ)

        # Exact cancellation can make generators disappear. Drop them before
        # computing minimal polynomials or estimating a primitive-element
        # degree; otherwise dead generators make later stages needlessly
        # expensive.
        used = numerator.free_symbols | denominator.free_symbols
        kept = [(gen, var) for gen, var in zip(generators, variables) if var in used]
        generators = tuple(gen for gen, _ in kept)
        variables = tuple(var for _, var in kept)

        minpolys: list[sp.Poly] = []
        degree_product = 1
        for generator, variable in zip(generators, variables):
            polynomial = minpoly_for(generator, variable)
            degree_product *= polynomial.degree()
            if degree_product > cfg.ALG_MODEL_MAX_DEGREE:
                return None
            minpolys.append(polynomial)

        return AlgebraicModel(
            expression=expr,
            generators=generators,
            variables=variables,
            minpolys=tuple(minpolys),
            numerator=numerator,
            denominator=denominator,
        )
    except EXACT_METHOD_ERRORS:
        return None


def _mod_poly(term: sp.Expr, var: sp.Symbol, divisor: sp.Poly) -> sp.Expr:
    """Reduce an expression modulo ``divisor`` without expanding powers."""
    if var not in term.free_symbols:
        return term
    if term == var:
        return var
    if term.is_Add:
        parts = [_mod_poly(arg, var, divisor) for arg in term.args]
        return sp.Poly(sp.Add(*parts), var, domain="EX").rem(divisor).as_expr()
    if term.is_Mul:
        acc = sp.Integer(1)
        for arg in term.args:
            part = _mod_poly(arg, var, divisor)
            acc = sp.Poly(acc * part, var, domain="EX").rem(divisor).as_expr()
            if acc == 0:
                break
        return acc
    if term.is_Pow and term.exp.is_Integer and term.exp.is_nonnegative:
        power = int(term.exp)
        base = _mod_poly(term.base, var, divisor)
        acc = sp.Integer(1)
        while power:
            if power & 1:
                acc = sp.Poly(acc * base, var, domain="EX").rem(divisor).as_expr()
            power >>= 1
            if power:
                base = sp.Poly(base * base, var, domain="EX").rem(divisor).as_expr()
        return acc
    return sp.Poly(term, var, domain="EX").rem(divisor).as_expr()


def _sequential_minpoly_remainder(polynomial: sp.Expr, model: AlgebraicModel) -> sp.Expr:
    remainder = polynomial
    for variable, defining_poly in zip(model.variables, model.minpolys):
        divisor = sp.Poly(defining_poly.as_expr(), variable, domain="EX")
        remainder = _mod_poly(remainder, variable, divisor)
        if remainder == 0:
            break
    return sp.cancel(remainder)


def _model_denominator_is_proven_nonzero(model: AlgebraicModel) -> bool:
    if model.denominator.is_Rational:
        return model.denominator != 0
    try:
        restored = model.denominator.xreplace(dict(zip(model.variables, model.generators)))
        value = to_number_field(restored)
        if isinstance(value, sp.AlgebraicNumber):
            return not all(c == 0 for c in value.coeffs())
        return value != 0
    except EXACT_METHOD_ERRORS:
        return False


def _generator_degree_product(expr: sp.Expr, limit: int) -> int | None:
    """Estimate the product of generator degrees, stopping after ``limit``."""
    try:
        modeled_expr = sp.cancel(expr)
        generators = extract_algebraic_generators(modeled_expr)
        if len(generators) > cfg.ALGEBRAIC_MAX_GENERATORS:
            return None
        degree_product = 1
        probe = sp.Dummy("a")
        for generator in generators:
            polynomial = minpoly_for(generator, probe)
            degree_product *= int(polynomial.degree())
            if degree_product > limit:
                return degree_product
        return degree_product
    except EXACT_METHOD_ERRORS:
        return None


def algebraic_relation_reduction_test(expr: sp.Expr) -> ZeroClassification:
    """Try to prove algebraic zero by minimal-polynomial remainder reduction.

    A zero remainder is a proof.  A nonzero remainder is deliberately returned
    as UNKNOWN because relations between distinct algebraic generators may not
    be represented by their separate univariate minimal polynomials.
    """
    if not stage_allowed(sp.sympify(expr), "minpoly"):
        return ZeroClassification(Verdict.UNKNOWN, "algebraic-remainder", detail="exact-method budget exceeded")

    model = build_algebraic_model(expr)
    if model is None:
        return ZeroClassification(
            Verdict.UNKNOWN,
            "algebraic-remainder",
            detail="could not build an algebraic generator model",
        )
    if not model.generators:
        verdict = Verdict.ZERO_PROVEN if model.numerator == 0 else Verdict.NONZERO_PROVEN
        return ZeroClassification(
            verdict,
            "algebraic-remainder",
            detail="exact rational expression",
            evidence="exact-algebraic-reduction",
        )
    try:
        remainder = _sequential_minpoly_remainder(model.numerator, model)
    except EXACT_METHOD_ERRORS as exc:
        return ZeroClassification(Verdict.UNKNOWN, "algebraic-remainder", detail=str(exc))
    if remainder == 0:
        if not _model_denominator_is_proven_nonzero(model):
            return ZeroClassification(
                Verdict.UNKNOWN,
                "algebraic-remainder",
                detail="numerator reduced to zero but denominator nonvanishing was not proved",
            )
        return ZeroClassification(
            Verdict.ZERO_PROVEN,
            "algebraic-remainder",
            detail=(f"numerator reduced to zero modulo {len(model.minpolys)} "
                    "generator minimal polynomial(s), with denominator proven nonzero"),
            evidence="minimal-polynomial-remainder",
        )
    return ZeroClassification(
        Verdict.UNKNOWN,
        "algebraic-remainder",
        detail="nonzero remainder is inconclusive because generator dependencies may remain",
    )


def common_number_field_test(expr: sp.Expr, model: AlgebraicModel | None = None) -> ZeroClassification:
    """Decide a closed algebraic expression in one exact number field.

    The conversion is gated by a degree-product budget because primitive
    element construction can grow rapidly.  Success yields a proof of either
    zero or nonzero.
    """
    if not stage_allowed(sp.sympify(expr), "minpoly"):
        return ZeroClassification(Verdict.UNKNOWN, "algebraic-common-field", detail="exact-method budget exceeded")

    expr = sp.sympify(expr)
    info = classify_algebraic_expression(expr)
    if not info.is_algebraic:
        return ZeroClassification(Verdict.UNKNOWN, "algebraic-common-field", detail=info.detail)
    model = model or build_algebraic_model(expr)
    degree_product = (model.degree_product if model is not None else
                      _generator_degree_product(expr, cfg.ALG_COMMON_MAX_DEGREE))
    if degree_product is None:
        return ZeroClassification(
            Verdict.UNKNOWN,
            "algebraic-common-field",
            detail="could not establish a safe primitive-element complexity bound",
        )
    if degree_product > cfg.ALG_COMMON_MAX_DEGREE:
        return ZeroClassification(
            Verdict.UNKNOWN,
            "algebraic-common-field",
            detail=(f"degree product {degree_product} exceeds common-field budget "
                    f"{cfg.ALG_COMMON_MAX_DEGREE}"),
        )
    try:
        value = to_number_field(expr)
        # AlgebraicNumber(0).is_zero is not consistently True across SymPy
        # versions, so inspect its exact coefficient representation as well.
        coeffs = tuple(value.coeffs()) if isinstance(value, sp.AlgebraicNumber) else ()
        zero = value.as_expr() == 0 if isinstance(value, sp.AlgebraicNumber) else value == 0
        if coeffs and all(c == 0 for c in coeffs):
            zero = True
        return ZeroClassification(
            Verdict.ZERO_PROVEN if zero else Verdict.NONZERO_PROVEN,
            "algebraic-common-field",
            detail=("common number-field representation is exactly zero" if zero
                    else "common number-field representation is exactly nonzero"),
            evidence="common-number-field",
        )
    except EXACT_METHOD_ERRORS as exc:
        return ZeroClassification(Verdict.UNKNOWN, "algebraic-common-field", detail=str(exc))





def clear_model_cache() -> None:
    """Clear cached algebraic generator metadata."""
    extract_algebraic_generators.cache_clear()
