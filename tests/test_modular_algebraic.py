import sympy as sp

from semialg.algebraic.modular import factor_univariate_qq_modular


def test_modular_factorization_exact_reconstruction_with_rational_unit():
    # (3/5) * (z-1)^2 * (z^2+z+1)
    z = sp.Symbol("z")
    expr = sp.expand(sp.Rational(3, 5) * (z - 1) ** 2 * (z**2 + z + 1))
    poly = sp.Poly(expr, z, domain=sp.QQ)
    coeffs = tuple(reversed(poly.all_coeffs()))
    result = factor_univariate_qq_modular(coeffs)
    assert result.certified
    assert sorted(m for _f, m in result.factors) == [1, 2]
