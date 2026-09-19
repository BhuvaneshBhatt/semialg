"""Certified localization and contraction primitives for GTZ decomposition.

This module implements the first two structural stages of the
Gianni--Trager--Zacharias primary-decomposition strategy over ``QQ``:

* choose and certify a maximal independent set ``U`` for ``QQ[x]/I``;
* extend ``I`` to the zero-dimensional ideal ``I QQ(U)[X]``;
* contract localized ideals by denominator clearing plus exact saturation;
* compute a finite saturation exponent by the ascending colon chain
  ``I:h^m`` and exact Groebner-basis stabilization.

All search data is independently replayable.  No finite iteration bound is
used for saturation stabilization; termination follows from Noetherianity.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass
from itertools import combinations

import sympy as sp

from .equality_ideal import EqualityIdealContext
from .groebner_utils import compute_groebner_basis
from .ideal_ops import canonical_qq_basis_uncached, qq_ideal_equal_uncached


def _qq_generators(
    generators: Sequence[sp.Expr], variables: Sequence[sp.Symbol]
) -> tuple[sp.Expr, ...]:
    vars_ = tuple(variables)
    try:
        polys = tuple(sp.Poly(sp.expand(f), *vars_, domain=sp.QQ) for f in generators)
    except (sp.PolynomialError, sp.polys.polyerrors.CoercionFailed) as exc:
        raise ValueError("GTZ localization currently requires QQ polynomial generators") from exc
    return tuple(p.as_expr() for p in polys if not p.is_zero)


def _canonical_basis(
    generators: Sequence[sp.Expr], variables: Sequence[sp.Symbol]
) -> tuple[sp.Expr, ...]:
    return canonical_qq_basis_uncached(generators, variables)


def _ideal_equal(
    left: Sequence[sp.Expr], right: Sequence[sp.Expr], variables: Sequence[sp.Symbol]
) -> bool:
    return qq_ideal_equal_uncached(left, right, variables)


def _leading_supports(context: EqualityIdealContext) -> tuple[frozenset[int], ...]:
    if context.groebner_basis is None:
        return tuple()
    order = context.groebner_basis.order
    out = []
    for poly in context.groebner_basis.polys:
        monom = poly.LM(order=order)
        support = frozenset(i for i, e in enumerate(monom.exponents) if int(e) > 0)
        if support:
            out.append(support)
    return tuple(out)


def _standard_monomial_count_from_basis(
    basis: sp.polys.polytools.GroebnerBasis,
    variables: Sequence[sp.Symbol],
) -> int | None:
    """Return finite quotient dimension from pure-power bounds when zero-dimensional."""
    vars_ = tuple(variables)
    if not vars_:
        return 1
    leading = [tuple(int(e) for e in p.LM(order=basis.order).exponents) for p in basis.polys]
    pure_bounds: list[int] = []
    for i in range(len(vars_)):
        powers = [exp[i] for exp in leading if exp[i] > 0 and sum(exp) == exp[i]]
        if not powers:
            return None
        pure_bounds.append(min(powers))
    count = 0
    for exponent in __import__("itertools").product(*(range(b) for b in pure_bounds)):
        if not any(all(e >= g for e, g in zip(exponent, lm, strict=True)) for lm in leading):
            count += 1
    return count


@dataclass(frozen=True)
class IndependentLocalizationCertificate:
    variables: tuple[sp.Symbol, ...]
    source_generators: tuple[sp.Expr, ...]
    dimension: int
    independent_variables: tuple[sp.Symbol, ...]
    dependent_variables: tuple[sp.Symbol, ...]
    leading_supports: tuple[frozenset[int], ...]
    localized_basis: tuple[sp.Expr, ...]
    quotient_dimension: int


@dataclass(frozen=True)
class IndependentLocalizationResult:
    independent_variables: tuple[sp.Symbol, ...]
    dependent_variables: tuple[sp.Symbol, ...]
    localized_basis: tuple[sp.Expr, ...]
    coefficient_field: object
    quotient_dimension: int
    certificate: IndependentLocalizationCertificate
    complete: bool = True
    method: str = "gtz_independent_localization"


def certified_independent_localization(
    generators: Sequence[sp.Expr], variables: Sequence[sp.Symbol]
) -> IndependentLocalizationResult:
    """Choose a maximal independent set and localize to dimension zero.

    For ``d = dim(QQ[x]/I)``, a size-``d`` subset ``U`` is certified independent
    when no leading monomial of an exact Groebner basis lies entirely in
    ``QQ[U]``.  Maximality then follows from ``|U| = d``.  The complement ``X``
    is localized over ``QQ(U)`` and its exact Groebner basis must be
    zero-dimensional.
    """
    vars_ = tuple(variables)
    gens = _qq_generators(tuple(generators), vars_)
    context = EqualityIdealContext.build(gens, vars_)
    if context.inconsistent:
        raise ValueError("the unit ideal has no GTZ independent localization")
    dimension = int(context.dimension)
    supports = _leading_supports(context)

    chosen: tuple[sp.Symbol, ...] | None = None
    index = {v: i for i, v in enumerate(vars_)}
    valid_subsets = []
    for subset in combinations(vars_, dimension):
        positions = frozenset(index[v] for v in subset)
        if not all(not support.issubset(positions) for support in supports):
            continue
        # Prefer parameter sets that minimize coefficient-field growth after
        # localization.  This is only a search heuristic; independence and
        # zero-dimensionality remain certified exactly below.
        parameter_degree = 0
        parameter_terms = 0
        for generator in gens:
            poly = sp.Poly(generator, *vars_, domain=sp.QQ)
            for monom, _coeff in poly.terms():
                degree = sum(monom[index[v]] for v in subset)
                parameter_degree = max(parameter_degree, int(degree))
                if degree:
                    parameter_terms += 1
        score = (parameter_degree, parameter_terms, tuple(index[v] for v in subset))
        valid_subsets.append((score, tuple(subset)))
    if valid_subsets:
        chosen = min(valid_subsets, key=lambda item: item[0])[1]
    if chosen is None:
        raise ArithmeticError("could not certify a maximal independent variable set")
    dependent = tuple(v for v in vars_ if v not in chosen)

    field = sp.QQ.frac_field(*chosen) if chosen else sp.QQ
    if dependent:
        localized = compute_groebner_basis(
            gens, dependent, order="grevlex", domain=field, modular=None
        )
        unit = len(localized.polys) == 1 and localized.polys[0].is_one
        if unit:
            quotient_dimension = 0
        else:
            quotient_dimension = _standard_monomial_count_from_basis(localized, dependent)
            if quotient_dimension is None:
                raise ArithmeticError("localized ideal was not certified zero-dimensional")
        localized_basis = tuple(sp.cancel(p.as_expr()) for p in localized.polys)
    else:
        localized_basis = tuple(sp.Integer(0) for _ in ())
        quotient_dimension = 1

    certificate = IndependentLocalizationCertificate(
        variables=vars_,
        source_generators=gens,
        dimension=dimension,
        independent_variables=chosen,
        dependent_variables=dependent,
        leading_supports=supports,
        localized_basis=localized_basis,
        quotient_dimension=quotient_dimension,
    )
    result = IndependentLocalizationResult(
        chosen,
        dependent,
        localized_basis,
        field,
        quotient_dimension,
        certificate,
    )
    if not verify_independent_localization_certificate(certificate):
        raise ArithmeticError("internal GTZ localization certificate verification failed")
    return result


def verify_independent_localization_certificate(
    certificate: IndependentLocalizationCertificate,
) -> bool:
    """Replay a GTZ independent-localization certificate from exact source data."""
    try:
        vars_ = tuple(certificate.variables)
        gens = _qq_generators(certificate.source_generators, vars_)
        context = EqualityIdealContext.build(gens, vars_)
        if context.inconsistent or int(context.dimension) != int(certificate.dimension):
            return False
        supports = _leading_supports(context)
        if supports != tuple(certificate.leading_supports):
            return False
        independent = tuple(certificate.independent_variables)
        dependent = tuple(certificate.dependent_variables)
        if set(independent).intersection(dependent):
            return False
        if tuple(v for v in vars_ if v in independent) != independent:
            return False
        if tuple(v for v in vars_ if v not in independent) != dependent:
            return False
        if len(independent) != context.dimension:
            return False
        positions = frozenset(vars_.index(v) for v in independent)
        if any(support.issubset(positions) for support in supports):
            return False
        field = sp.QQ.frac_field(*independent) if independent else sp.QQ
        if dependent:
            localized = compute_groebner_basis(
                gens, dependent, order="grevlex", domain=field, modular=None
            )
            basis = tuple(sp.cancel(p.as_expr()) for p in localized.polys)
            if basis != tuple(certificate.localized_basis):
                return False
            unit = len(localized.polys) == 1 and localized.polys[0].is_one
            count = 0 if unit else _standard_monomial_count_from_basis(localized, dependent)
            return count is not None and count == certificate.quotient_dimension
        return certificate.localized_basis == tuple() and certificate.quotient_dimension == 1
    except (ArithmeticError, ValueError, TypeError, sp.PolynomialError):
        return False


def _clear_fraction_field_denominators(
    generators: Sequence[sp.Expr],
    independent: Sequence[sp.Symbol],
    dependent: Sequence[sp.Symbol],
) -> tuple[tuple[sp.Expr, ...], sp.Expr]:
    """Clear coefficient denominators from ``QQ(U)[X]`` generators exactly."""
    all_vars = (*independent, *dependent)
    cleared: list[sp.Expr] = []
    denominator_product = sp.Integer(1)
    for generator in generators:
        expr = sp.cancel(generator)
        num, den = sp.fraction(expr)
        den = (
            sp.Poly(den, *independent, domain=sp.QQ).as_expr()
            if independent
            else sp.sympify(sp.QQ.convert(den))
        )
        num = sp.expand(num)
        if den == 0:
            raise ZeroDivisionError("zero coefficient denominator in localized ideal")
        denominator_product = sp.expand(denominator_product * den)
        cleared.append(sp.Poly(num, *all_vars, domain=sp.QQ).as_expr())
    if independent:
        _, primitive = sp.Poly(denominator_product, *independent, domain=sp.QQ).primitive()
        denominator_product = sp.expand(primitive.as_expr())
    else:
        denominator_product = sp.Integer(1)
    return tuple(cleared), denominator_product


def _saturate_exact(
    generators: Sequence[sp.Expr], splitter: sp.Expr, variables: Sequence[sp.Symbol]
) -> tuple[sp.Expr, ...]:
    """Compute ``I:splitter**infinity`` by Rabinowitsch elimination."""
    vars_ = tuple(variables)
    if sp.expand(splitter) == 0:
        raise ValueError("cannot saturate by zero")
    u = sp.Dummy("gtz_sat")
    gb = sp.groebner(
        (*generators, 1 - u * sp.expand(splitter)),
        u,
        *vars_,
        order="lex",
        domain=sp.QQ,
    )
    eliminated = tuple(
        sp.expand(p.as_expr()) for p in gb.polys if u not in p.as_expr().free_symbols
    )
    return _canonical_basis(eliminated, vars_)


@dataclass(frozen=True)
class LocalizationContractionCertificate:
    variables: tuple[sp.Symbol, ...]
    independent_variables: tuple[sp.Symbol, ...]
    dependent_variables: tuple[sp.Symbol, ...]
    localized_generators: tuple[sp.Expr, ...]
    canonical_localized_basis: tuple[sp.Expr, ...]
    cleared_generators: tuple[sp.Expr, ...]
    denominator_product: sp.Expr
    contracted_generators: tuple[sp.Expr, ...]


@dataclass(frozen=True)
class LocalizationContractionResult:
    generators: tuple[sp.Expr, ...]
    denominator_product: sp.Expr
    certificate: LocalizationContractionCertificate
    complete: bool = True
    method: str = "gtz_localization_contraction"


def contract_localized_ideal(
    localized_generators: Sequence[sp.Expr],
    independent_variables: Sequence[sp.Symbol],
    dependent_variables: Sequence[sp.Symbol],
) -> LocalizationContractionResult:
    """Contract an ideal from ``QQ(U)[X]`` to ``QQ[U,X]`` exactly.

    Coefficient denominators are cleared, and contraction is computed as the
    saturation of the cleared ideal by their product.
    """
    independent = tuple(independent_variables)
    dependent = tuple(dependent_variables)
    variables = (*independent, *dependent)
    localized = tuple(sp.cancel(g) for g in localized_generators)
    field = sp.QQ.frac_field(*independent) if independent else sp.QQ
    if dependent:
        local_gb = sp.groebner(localized, *dependent, order="grevlex", domain=field)
        canonical_localized = tuple(sp.cancel(p.as_expr()) for p in local_gb.polys)
    else:
        canonical_localized = localized
    cleared, denominator = _clear_fraction_field_denominators(
        canonical_localized, independent, dependent
    )
    if denominator == 1:
        contracted = _canonical_basis(cleared, variables)
    else:
        contracted = _saturate_exact(cleared, denominator, variables)
    certificate = LocalizationContractionCertificate(
        variables,
        independent,
        dependent,
        localized,
        canonical_localized,
        cleared,
        denominator,
        contracted,
    )
    result = LocalizationContractionResult(contracted, denominator, certificate)
    if not verify_localization_contraction_certificate(certificate):
        raise ArithmeticError("internal GTZ contraction certificate verification failed")
    return result


def verify_localization_contraction_certificate(
    certificate: LocalizationContractionCertificate,
) -> bool:
    """Replay denominator clearing, saturation, and localized re-extension exactly."""
    try:
        independent = tuple(certificate.independent_variables)
        dependent = tuple(certificate.dependent_variables)
        variables = tuple(certificate.variables)
        if tuple(certificate.variables) != (*independent, *dependent):
            return False
        if len(set(variables)) != len(variables) or not all(
            isinstance(v, sp.Symbol) for v in variables
        ):
            return False
        if set(independent).intersection(dependent):
            return False
        localized = tuple(sp.cancel(g) for g in certificate.localized_generators)
        field = sp.QQ.frac_field(*independent) if independent else sp.QQ
        if dependent:
            local_gb = sp.groebner(localized, *dependent, order="grevlex", domain=field)
            canonical_localized = tuple(sp.cancel(p.as_expr()) for p in local_gb.polys)
        else:
            canonical_localized = localized
        if canonical_localized != tuple(certificate.canonical_localized_basis):
            return False
        cleared, denominator = _clear_fraction_field_denominators(
            canonical_localized, independent, dependent
        )
        if cleared != tuple(certificate.cleared_generators):
            return False
        if sp.expand(denominator - certificate.denominator_product) != 0:
            return False
        if denominator == 0 or denominator.free_symbols - set(independent):
            return False
        for generator in cleared:
            try:
                sp.Poly(generator, *variables, domain=sp.QQ)
            except (sp.PolynomialError, sp.polys.polyerrors.CoercionFailed):
                return False
        expected = (
            _canonical_basis(cleared, certificate.variables)
            if denominator == 1
            else _saturate_exact(cleared, denominator, certificate.variables)
        )
        if expected != tuple(certificate.contracted_generators):
            return False

        # Re-extension must recover the supplied localized ideal exactly.
        if dependent:
            lhs = sp.groebner(localized, *dependent, order="grevlex", domain=field)
            rhs = sp.groebner(expected, *dependent, order="grevlex", domain=field)
            return tuple(sp.cancel(p.as_expr()) for p in lhs.polys) == tuple(
                sp.cancel(p.as_expr()) for p in rhs.polys
            )
        return True
    except (
        ArithmeticError,
        ValueError,
        TypeError,
        ZeroDivisionError,
        sp.PolynomialError,
        sp.polys.polyerrors.CoercionFailed,
    ):
        return False


def _ideal_intersection(
    left: Sequence[sp.Expr], right: Sequence[sp.Expr], variables: Sequence[sp.Symbol]
) -> tuple[sp.Expr, ...]:
    """Return the exact intersection of two ideals by elimination."""
    vars_ = tuple(variables)
    t = sp.Dummy("gtz_ideal_intersection")
    generators = tuple(t * sp.expand(g) for g in left) + tuple(
        (1 - t) * sp.expand(g) for g in right
    )
    gb = sp.groebner(generators, t, *vars_, order="lex", domain=sp.QQ)
    eliminated = tuple(
        sp.expand(p.as_expr()) for p in gb.polys if t not in p.as_expr().free_symbols
    )
    return _canonical_basis(eliminated, vars_)


def _ideal_intersection_with_principal(
    generators: Sequence[sp.Expr], h: sp.Expr, variables: Sequence[sp.Symbol]
) -> tuple[sp.Expr, ...]:
    """Return ``I intersect <h>`` by exact elimination."""
    vars_ = tuple(variables)
    t = sp.Dummy("gtz_intersection")
    ideal_generators = tuple(t * sp.expand(g) for g in generators)
    gb = sp.groebner(
        (*ideal_generators, (1 - t) * sp.expand(h)),
        t,
        *vars_,
        order="lex",
        domain=sp.QQ,
    )
    eliminated = tuple(
        sp.expand(p.as_expr()) for p in gb.polys if t not in p.as_expr().free_symbols
    )
    return _canonical_basis(eliminated, vars_)


def _colon_once(
    generators: Sequence[sp.Expr], h: sp.Expr, variables: Sequence[sp.Symbol]
) -> tuple[sp.Expr, ...]:
    """Return ``I:<h>`` from ``(I intersect <h>)/h`` exactly."""
    vars_ = tuple(variables)
    intersection = _ideal_intersection_with_principal(generators, h, vars_)
    quotients = []
    divisor = sp.Poly(sp.expand(h), *vars_, domain=sp.QQ)
    for generator in intersection:
        quotient, remainder = sp.div(
            sp.Poly(generator, *vars_, domain=sp.QQ), divisor, domain=sp.QQ
        )
        if not remainder.is_zero:
            raise ArithmeticError("principal intersection generator was not divisible by h")
        quotients.append(quotient.as_expr())
    return _canonical_basis(quotients, vars_)


@dataclass(frozen=True)
class SaturationStabilizationCertificate:
    variables: tuple[sp.Symbol, ...]
    source_generators: tuple[sp.Expr, ...]
    splitter: sp.Expr
    exponent: int
    colon_chain: tuple[tuple[sp.Expr, ...], ...]
    saturated_generators: tuple[sp.Expr, ...]
    companion_generators: tuple[sp.Expr, ...]


@dataclass(frozen=True)
class SaturationStabilizationResult:
    generators: tuple[sp.Expr, ...]
    exponent: int
    companion_generators: tuple[sp.Expr, ...]
    certificate: SaturationStabilizationCertificate
    complete: bool = True
    method: str = "gtz_colon_stabilization"


def saturation_stabilization(
    generators: Sequence[sp.Expr], splitter: sp.Expr, variables: Sequence[sp.Symbol]
) -> SaturationStabilizationResult:
    """Compute the least detected ``m`` with ``I:h^m = I:h^(m+1)`` exactly.

    The ascending chain is followed without a fixed attempt limit.  Once two
    consecutive canonical Groebner bases agree, the stable ideal equals
    ``I:h**infinity``.  Noetherianity guarantees termination.
    """
    vars_ = tuple(variables)
    gens = _qq_generators(tuple(generators), vars_)
    h = sp.Poly(sp.expand(splitter), *vars_, domain=sp.QQ).as_expr()
    if h == 0:
        raise ValueError("saturation splitter must be nonzero")
    current = _canonical_basis(gens, vars_)
    chain = [current]
    exponent = 0
    while True:
        nxt = _colon_once(current, h, vars_)
        exponent += 1
        chain.append(nxt)
        if nxt == current:
            # I:h^(m-1) = I:h^m, so stabilization exponent is m-1.
            stable_exponent = exponent - 1
            stable = current
            break
        current = nxt
    rabinowitsch = _saturate_exact(gens, h, vars_)
    if stable != rabinowitsch:
        raise ArithmeticError("colon-chain stabilization disagrees with exact saturation")
    power = sp.expand(h**stable_exponent)
    companion = _canonical_basis((*gens, power), vars_)
    reconstructed = _ideal_intersection(stable, companion, vars_)
    source_basis = _canonical_basis(gens, vars_)
    if reconstructed != source_basis:
        raise ArithmeticError("GTZ saturation split failed exact source reconstruction")
    certificate = SaturationStabilizationCertificate(
        vars_, gens, h, stable_exponent, tuple(chain), stable, companion
    )
    result = SaturationStabilizationResult(stable, stable_exponent, companion, certificate)
    if not verify_saturation_stabilization_certificate(certificate):
        raise ArithmeticError("internal saturation-stabilization certificate failed")
    return result


def verify_saturation_stabilization_certificate(
    certificate: SaturationStabilizationCertificate,
) -> bool:
    """Replay the colon chain, saturation, and GTZ split identity exactly."""
    try:
        vars_ = tuple(certificate.variables)
        gens = _qq_generators(certificate.source_generators, vars_)
        h = sp.Poly(certificate.splitter, *vars_, domain=sp.QQ).as_expr()
        if h == 0 or certificate.exponent < 0:
            return False
        current = _canonical_basis(gens, vars_)
        replay = [current]
        for _ in range(certificate.exponent + 1):
            nxt = _colon_once(current, h, vars_)
            replay.append(nxt)
            current = nxt
        if tuple(replay) != tuple(certificate.colon_chain):
            return False
        # The last equality certifies stabilization at exponent m.
        if replay[-1] != replay[-2]:
            return False
        stable = replay[-1]
        if stable != tuple(certificate.saturated_generators):
            return False
        if stable != _saturate_exact(gens, h, vars_):
            return False
        companion = _canonical_basis((*gens, sp.expand(h**certificate.exponent)), vars_)
        if companion != tuple(certificate.companion_generators):
            return False
        reconstructed = _ideal_intersection(stable, companion, vars_)
        return reconstructed == _canonical_basis(gens, vars_)
    except (
        ArithmeticError,
        ValueError,
        TypeError,
        sp.PolynomialError,
        sp.polys.polyerrors.CoercionFailed,
    ):
        return False


__all__ = [
    "IndependentLocalizationCertificate",
    "IndependentLocalizationResult",
    "LocalizationContractionCertificate",
    "LocalizationContractionResult",
    "SaturationStabilizationCertificate",
    "SaturationStabilizationResult",
    "certified_independent_localization",
    "contract_localized_ideal",
    "saturation_stabilization",
    "verify_independent_localization_certificate",
    "verify_localization_contraction_certificate",
    "verify_saturation_stabilization_certificate",
]
