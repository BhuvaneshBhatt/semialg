"""Generic algebraic fiber degree of polynomial and rational maps."""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass

import sympy as sp

from .algebraic.rational_univariate.quotient import (
    _leading_exponent_grevlex,
    _standard_exponent_count,
)
from .normalization import normalize_variables


@dataclass(frozen=True)
class ParametricMapDegree:
    """Generic complex fiber degree of a rational map in characteristic zero."""

    degree: int | None
    generically_finite: bool
    source_dimension: int
    target_dimension: int
    equations: tuple[sp.Expr, ...]


def parametric_map_degree(
    mapping: Sequence[object],
    parameters: Sequence[sp.Symbol | str],
) -> ParametricMapDegree:
    """Return the generic algebraic fiber degree of a rational map.

    The computation forms ``F(u) = F(t)`` with ``t`` treated as generic
    coefficients.  For a zero-dimensional generic fiber, the quotient-algebra
    dimension is the degree.  This is an algebraic generic degree; inequalities
    restricting a real parameter domain are intentionally not folded into it.
    """

    params = tuple(normalize_variables(parameters, append_context_symbols=False))
    maps = tuple(map(sp.sympify, mapping))
    if not params:
        return ParametricMapDegree(1, True, 0, len(maps), ())
    copies = tuple(sp.Dummy(f"u{i + 1}", real=True) for i in range(len(params)))
    repl = dict(zip(params, copies, strict=True))
    equations: list[sp.Expr] = []
    for expr in maps:
        difference = sp.together(expr.xreplace(repl) - expr)
        numerator, denominator = sp.fraction(difference)
        if denominator == 0:
            raise ValueError("invalid rational map")
        try:
            sp.Poly(numerator, *copies)
        except (sp.PolynomialError, TypeError, ValueError) as exc:
            raise NotImplementedError(
                "map degree requires rational functions of the parameters"
            ) from exc
        equations.append(sp.expand(numerator))
    try:
        coefficient_domain = sp.QQ.frac_field(*params)
        basis = sp.groebner(equations, *copies, order="grevlex", domain=coefficient_domain)
    except (sp.PolynomialError, TypeError, ValueError, NotImplementedError):
        return ParametricMapDegree(None, False, len(params), len(maps), tuple(equations))
    if not basis.is_zero_dimensional:
        return ParametricMapDegree(None, False, len(params), len(maps), tuple(equations))
    leading = tuple(_leading_exponent_grevlex(poly) for poly in basis.polys)
    degree = _standard_exponent_count(leading, len(copies))
    return ParametricMapDegree(
        int(degree),
        True,
        len(params),
        len(maps),
        tuple(equations),
    )


__all__ = ["ParametricMapDegree", "parametric_map_degree"]
