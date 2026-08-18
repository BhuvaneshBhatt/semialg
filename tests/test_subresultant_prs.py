from __future__ import annotations

import sympy as sp

from semialg import principal_subresultant_coefficients, subresultant_prs
from semialg.algebraic import SubresultantPRSResult


def test_subresultant_prs_uses_exact_univariate_polys():
    x = sp.symbols("x")
    result = subresultant_prs(x**3 - 2 * x + 1, x**2 - 1, x, domain=sp.QQ)

    assert isinstance(result, SubresultantPRSResult)
    assert result.variable == x
    assert result.polynomials[0] == sp.Poly(x**3 - 2 * x + 1, x, domain=sp.QQ)
    assert result.polynomials[1] == sp.Poly(x**2 - 1, x, domain=sp.QQ)
    assert sp.factor(result.resultant) == sp.resultant(x**3 - 2 * x + 1, x**2 - 1, x)
    assert result.principal_coefficients == tuple(
        sp.factor(poly.LC()) for poly in result.polynomials
    )


def test_principal_subresultant_coefficients_matches_result_object():
    x = sp.symbols("x")
    f = x**4 - 1
    g = x**3 - x
    result = subresultant_prs(f, g, x)
    assert principal_subresultant_coefficients(f, g, x) == result.principal_coefficients


def test_subresultant_prs_detects_common_factor_by_zero_resultant():
    x = sp.symbols("x")
    result = subresultant_prs((x - 1) * (x + 2), (x - 1) * (x + 3), x)
    assert result.resultant == 0
    assert any(
        poly.degree() == 1 and sp.rem(poly.as_expr(), x - 1, x) == 0 for poly in result.polynomials
    )
