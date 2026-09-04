"""Exact specialized decision procedures for polynomial negativity/nonnegativity.

The backend implements the critical-value core of the Zeng decision strategy:
cheap exact witnesses are tried first, then coercive polynomials are reduced to
an exact zero-dimensional gradient system.  A result is marked complete only
when these hypotheses prove that every global minimum is represented by the
critical system. Unsupported cases remain explicitly incomplete.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass

import sympy as sp

from ._linear_relations import certified_sign
from .algebraic.rational_univariate import RationalUnivariateError
from .solve.zero_dimensional import is_zero_dimensional, solve_zero_dimensional_system
from .witness_heuristics import _certified_negative_point, find_negative_witness_fast


@dataclass(frozen=True)
class PolynomialNegativityResult:
    variables: tuple[sp.Symbol, ...]
    has_negative_point: bool | None
    point: tuple[tuple[sp.Symbol, sp.Expr], ...] | None
    complete: bool
    method: str
    notes: tuple[str, ...] = ()

    @property
    def assignment(self) -> dict[sp.Symbol, sp.Expr] | None:
        return dict(self.point) if self.point is not None else None

    @property
    def nonnegative(self) -> bool | None:
        if not self.complete or self.has_negative_point is None:
            return None
        return not self.has_negative_point

    def require_decision(self) -> bool:
        if not self.complete or self.has_negative_point is None:
            raise NotImplementedError("polynomial negativity decision is incomplete")
        return self.has_negative_point


def _coercive_even_leading_form(polynomial: sp.Expr, variables: tuple[sp.Symbol, ...]) -> bool:
    """Certify coercivity from a positive even-monomial leading form.

    This intentionally recognizes a strict subset of coercive polynomials: all
    top-degree monomials have nonnegative rational coefficients and even
    exponents, and every variable has a positive pure top-degree monomial.
    """
    poly = sp.Poly(sp.expand(polynomial), *variables, domain=sp.QQ)
    degree = int(poly.total_degree())
    if degree <= 0 or degree % 2:
        return False
    top = [(exp, sp.Rational(coeff)) for exp, coeff in poly.terms() if sum(exp) == degree]
    if not top or any(coeff < 0 for _, coeff in top):
        return False
    if any(any(int(power) % 2 for power in exponent) for exponent, _ in top):
        return False
    for index in range(len(variables)):
        pure = tuple(degree if j == index else 0 for j in range(len(variables)))
        if not any(exp == pure and coeff > 0 for exp, coeff in top):
            return False
    return True


def zeng_negative_point(
    polynomial: sp.Expr,
    variables: Sequence[sp.Symbol],
    *,
    random_lines: int = 8,
    seed: int = 1234,
) -> PolynomialNegativityResult:
    """Decide whether a rational polynomial is negative somewhere when certified.

    The implementation uses the exact critical-value reduction underlying the
    Zeng semidefinite-polynomial decision method.  On the currently certified
    coercive fragment, global minima are attained and the gradient ideal is
    solved exactly by the package's RUR backend.  Fast odd-degree/ray and
    random-line searches may prove the positive (negative-point) side earlier;
    they never prove nonnegativity.
    """
    vars_ = tuple(variables)
    if random_lines < 0:
        raise ValueError("random_lines must be nonnegative")
    if not vars_:
        sign = certified_sign(sp.sympify(polynomial))
        return PolynomialNegativityResult(
            vars_, sign == -1, tuple() if sign == -1 else None, sign is not None, "zeng_constant"
        )
    expr = sp.expand(polynomial)
    try:
        poly = sp.Poly(expr, *vars_, domain=sp.QQ)
    except (sp.PolynomialError, ValueError, TypeError) as exc:
        raise RationalUnivariateError("Zeng decision requires a rational polynomial") from exc

    if poly.is_zero:
        return PolynomialNegativityResult(vars_, False, None, True, "zeng", ("zero_polynomial",))
    constant_sign = certified_sign(expr) if not expr.free_symbols.intersection(vars_) else None
    if constant_sign is not None:
        point = tuple((v, sp.Integer(0)) for v in vars_) if constant_sign == -1 else None
        return PolynomialNegativityResult(vars_, constant_sign == -1, point, True, "zeng_constant")

    fast = find_negative_witness_fast(expr, vars_, random_lines=random_lines, seed=seed)
    if fast.found:
        return PolynomialNegativityResult(vars_, True, fast.point, True, f"zeng+{fast.method}")

    if not _coercive_even_leading_form(expr, vars_):
        return PolynomialNegativityResult(
            vars_, None, None, False, "zeng", ("coercivity_not_certified",)
        )

    gradient = tuple(sp.expand(sp.diff(expr, variable)) for variable in vars_)
    gradient = tuple(g for g in gradient if g != 0)
    if not gradient:
        sign = certified_sign(expr.subs({v: 0 for v in vars_}))
        return PolynomialNegativityResult(vars_, sign == -1, None, sign is not None, "zeng")
    if not is_zero_dimensional(gradient, vars_):
        return PolynomialNegativityResult(
            vars_, None, None, False, "zeng", ("positive_dimensional_gradient_locus",)
        )
    try:
        points = solve_zero_dimensional_system(gradient, vars=vars_, real=True)
    except RationalUnivariateError:
        return PolynomialNegativityResult(
            vars_, None, None, False, "zeng", ("critical_rur_failure",)
        )

    undecidable = False
    for point in points:
        verified = _certified_negative_point(expr, vars_, point)
        if verified is not None:
            return PolynomialNegativityResult(
                vars_,
                True,
                verified,
                True,
                "zeng_critical_values",
                (f"critical_points={len(points)}",),
            )
        assignment = dict(zip(vars_, point, strict=True))
        if certified_sign(expr.subs(assignment)) is None:
            undecidable = True
    if undecidable:
        return PolynomialNegativityResult(
            vars_, None, None, False, "zeng_critical_values", ("critical_value_sign_undecidable",)
        )
    return PolynomialNegativityResult(
        vars_, False, None, True, "zeng_critical_values", (f"critical_points={len(points)}",)
    )


def polynomial_nonnegative(
    polynomial: sp.Expr,
    variables: Sequence[sp.Symbol],
    **kwargs,
) -> bool:
    """Return a certified global nonnegativity decision or raise if unsupported."""
    result = zeng_negative_point(polynomial, variables, **kwargs)
    return not result.require_decision()


def find_negative_point(
    polynomial: sp.Expr,
    variables: Sequence[sp.Symbol],
    **kwargs,
) -> dict[sp.Symbol, sp.Expr] | None:
    """Return a certified negative point or ``None`` for certified nonnegativity."""
    result = zeng_negative_point(polynomial, variables, **kwargs)
    negative = result.require_decision()
    return result.assignment if negative else None


__all__ = [
    "PolynomialNegativityResult",
    "find_negative_point",
    "polynomial_nonnegative",
    "zeng_negative_point",
]
