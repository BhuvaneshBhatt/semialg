import sympy as sp

from semialg import semialgebraic_minimize
from semialg.interval_decomposition import finite_real_roots

x = sp.Symbol("x", real=True)


def test_algebraic_coefficients():
    polynomial = (x - sp.sqrt(2)) ** 2 * (x + sp.sqrt(3))
    assert finite_real_roots(polynomial, x) == (-sp.sqrt(3), sp.sqrt(2))


def test_real_branches():
    assert finite_real_roots((x - sp.sqrt(2)) * (x * x + 1), x) == (sp.sqrt(2),)
    assert finite_real_roots(sp.S.Zero, x) == ()
    assert finite_real_roots(sp.S(3), x) == ()


def test_algebraic_constraint_endpoint():
    result = semialgebraic_minimize(x, x * x >= 2, (x,), return_result=True)
    assert result.value == -sp.oo and result.attained is False
