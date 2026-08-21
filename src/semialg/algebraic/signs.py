from __future__ import annotations

from collections.abc import Sequence

import sympy as sp

from ..exact_arithmetic import exact_sign
from .cache import CACHE, expr_key, sample_expr_key
from .roots import refine_isol_intv
from .samples import AlgebraicRoot, RationalSample, Sample, sample_to_expr


def _check_samples(samples: Sequence[Sample]) -> None:
    for sample in samples:
        if not isinstance(sample, (RationalSample, AlgebraicRoot)):
            raise TypeError("sign_at_sample requires explicit semialg Sample objects")


def _sign_at_opaque_last_root(
    expr: sp.Expr,
    gens: Sequence[sp.Symbol],
    samples: Sequence[Sample],
) -> int | None:
    """Evaluate a polynomial sign at an opaque final algebraic sample.

    When an algebraic-coefficient fiber is represented by its defining
    polynomial plus an isolating interval, converting the sample to a SymPy
    algebraic expression may be impossible.  The sign of another polynomial is
    nevertheless determined exactly: detect a common root with a gcd, otherwise
    refine the sample interval until the target polynomial has no root there and
    evaluate its constant sign at a rational midpoint.
    """

    if not samples or not isinstance(samples[-1], AlgebraicRoot):
        return None
    root = samples[-1]
    if root.root_expr is not None:
        return None
    var = gens[len(samples) - 1]
    substitutions = {
        gens[index]: sample_to_expr(sample) for index, sample in enumerate(samples[:-1])
    }
    target_expr = sp.expand(expr.subs(substitutions))
    try:
        target = sp.Poly(target_expr, var, extension=True)
        defining_expr = root.polynomial.as_expr().subs(root.variable, var)
        defining = sp.Poly(defining_expr, var, extension=True)
    except (sp.PolynomialError, ValueError, TypeError):
        return None
    if target.degree() <= 0:
        try:
            return exact_sign(target.as_expr())
        except ValueError:
            return None
    try:
        common = sp.gcd(target, defining)
        if (
            common.degree() > 0
            and int(common.count_roots(root.interval.left, root.interval.right)) > 0
        ):
            return 0
    except (
        NotImplementedError,
        sp.PolynomialError,
        sp.polys.polyerrors.DomainError,
        ValueError,
        TypeError,
    ):
        return None

    refined = root
    for _ in range(48):
        left, right = refined.interval.left, refined.interval.right
        try:
            if int(target.count_roots(left, right)) == 0:
                midpoint = sp.Rational(left + right, 2)
                return exact_sign(target.eval(midpoint))
        except (
            NotImplementedError,
            sp.PolynomialError,
            sp.polys.polyerrors.DomainError,
            ValueError,
            TypeError,
        ):
            return None
        refined = refine_isol_intv(refined, steps=2)
    return None


def sign_at_sample(poly: sp.Poly | sp.Expr, samples: Sequence[Sample]) -> int:
    """Return the exact sign of ``poly`` at an explicit algebraic sample tuple.

    Callers must pass ``RationalSample`` or ``AlgebraicRoot`` objects.
    The implementation first uses exact substitution through SymPy algebraic
    expressions and then falls back to high-precision evaluation only after
    symbolic sign checks are inconclusive.
    """

    _check_samples(samples)
    if isinstance(poly, sp.Poly):
        expr = poly.as_expr()
        gens = poly.gens
    else:
        expr = sp.expand(poly)
        gens = tuple(sorted(expr.free_symbols, key=lambda sym: sym.name))
    if len(samples) > len(gens):
        raise ValueError("more sample coordinates were supplied than polynomial generators")
    cache_key = (
        expr_key(sp.sympify(expr)),
        tuple(sp.srepr(g) for g in gens[: len(samples)]),
        tuple(sample_expr_key(s) for s in samples),
    )
    cached = CACHE.signs.get(cache_key)
    if cached is not None:
        CACHE.stats.sign_hits += 1
        return cached
    CACHE.stats.sign_misses += 1
    opaque_result = _sign_at_opaque_last_root(expr, gens, samples)
    if opaque_result is not None:
        CACHE.signs.put(cache_key, opaque_result)
        return opaque_result
    substitutions = {gens[index]: sample_to_expr(sample) for index, sample in enumerate(samples)}
    value = sp.simplify(sp.expand(expr.subs(substitutions)))
    result: int | None = None
    if value == 0:
        result = 0
    else:
        sign = sp.sign(value)
        if sign in (-1, 0, 1):
            result = int(sign)
        elif value.is_positive:
            result = 1
        elif value.is_negative:
            result = -1
    if result is None:
        # Keep the exact API exact.  Algebraic-number conversion often exposes
        # a RootOf/explicit algebraic representation whose sign SymPy can prove.
        # If it cannot, decline rather than guessing from fixed precision.
        try:
            algebraic = sp.polys.numberfields.to_number_field(value).to_root()
            if algebraic == 0 or algebraic.is_zero is True:
                result = 0
            elif algebraic.is_positive is True:
                result = 1
            elif algebraic.is_negative is True:
                result = -1
            else:
                exact_sign = sp.sign(algebraic)
                if exact_sign in (-1, 0, 1):
                    result = int(exact_sign)
        except (TypeError, ValueError, NotImplementedError, sp.PolynomialError):
            pass
    if result is not None:
        CACHE.signs.put(cache_key, result)
        return result
    raise ValueError(f"could not determine sign of {sp.sstr(expr)} at sample {samples!r}")
