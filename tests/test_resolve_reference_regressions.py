"""Exact-real Resolve regressions adapted from the supplied historical corpus.

Only exact real semialgebraic cases are retained.  Numeric-precision,
complex-default, integer/prime-domain, and system-option cases are outside this
package's Resolve contract and are intentionally excluded.
"""

import sympy as sp

from semialg import resolve_formula
from semialg.formula import ParsedPrenexFormula, parse_formula


def _exists(expr: sp.Expr, *variables: sp.Symbol):
    return ParsedPrenexFormula(
        vars=tuple(variables),
        quantifiers=tuple(("exists", var) for var in variables),
        matrix=parse_formula(expr),
        matrix_expr=expr,
    )


def test_resolve_exact_circle_has_point_below_diagonal() -> None:
    x, y = sp.symbols("x y", real=True)
    expr = sp.And(sp.Eq(x**2 + y**2, 1), x + y < 1)
    assert resolve_formula(_exists(expr, x, y)) is True


def test_resolve_sum_of_squares_cannot_be_negative() -> None:
    x, y = sp.symbols("x y", real=True)
    assert resolve_formula(_exists(x**2 + y**2 < 0, x, y)) is False


def test_resolve_exact_univariate_open_interval() -> None:
    x = sp.Symbol("x", real=True)
    assert resolve_formula(_exists(sp.And(x**2 < 2, x > 1), x)) is True


def test_resolve_parabola_meets_unit_disk() -> None:
    x, y = sp.symbols("x y", real=True)
    expr = sp.And(x**2 + y**2 <= 1, sp.Eq(y, x**2))
    assert resolve_formula(_exists(expr, x, y)) is True


def test_resolve_quartic_constraints_have_real_witness() -> None:
    x, y = sp.symbols("x y", real=True)
    expr = sp.And(sp.Eq(x**4 + y**4, 1), y**2 < x**2 + 2, x > 0, y > sp.Rational(1, 2))
    assert resolve_formula(_exists(expr, x, y)) is True


def test_resolve_even_power_zero_boundary_distinguishes_strictness() -> None:
    x, y = sp.symbols("x y", real=True)
    nonnegative = (x - y + 1) ** 2 + (2 * x - 3 * y - 3) ** 4
    assert resolve_formula(_exists(nonnegative < 0, x, y)) is False
    assert resolve_formula(_exists(nonnegative <= 0, x, y)) is True


def test_resolve_exact_annulus_boundary_is_nonempty() -> None:
    x, y = sp.symbols("x y", real=True)
    radius2 = x**2 + y**2
    assert resolve_formula(_exists(sp.And(radius2 >= 2, radius2**2 <= 4), x, y)) is True
