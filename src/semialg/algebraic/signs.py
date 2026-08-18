from __future__ import annotations

from collections.abc import Sequence

import sympy as sp

from .samples import AlgebraicRoot, RationalSample, Sample, sample_to_expr


def _check_samples(samples: Sequence[Sample]) -> None:
    for sample in samples:
        if not isinstance(sample, (RationalSample, AlgebraicRoot)):
            raise TypeError("sign_at_sample requires explicit semialg Sample objects")


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
    substitutions = {gens[index]: sample_to_expr(sample) for index, sample in enumerate(samples)}
    value = sp.simplify(sp.expand(expr.subs(substitutions)))
    if value == 0:
        return 0
    sign = sp.sign(value)
    if sign in (-1, 0, 1):
        return int(sign)
    if value.is_positive:
        return 1
    if value.is_negative:
        return -1
    numeric = sp.N(value, 160)
    if numeric == 0:
        return 0
    if numeric > 0:
        return 1
    if numeric < 0:
        return -1
    raise ValueError(f"could not determine sign of {sp.sstr(expr)} at sample {samples!r}")
