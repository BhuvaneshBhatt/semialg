import sympy as sp

from semialg import equivalent, is_satisfiable
from semialg.algebraic.cache import clear_algebraic_caches, expr_key
from semialg.cad_algorithms.performance_cache import clear_cad_caches


def _pollute_caches():
    u, v = sp.symbols("u v", real=True)
    assert is_satisfiable(sp.And(u**2 + v**2 < 3, u > 1), [u, v])
    assert equivalent(sp.Eq((u - 1) * (u + 1), 0), sp.Or(sp.Eq(u, 1), sp.Eq(u, -1)), [u])


def test_exact_decision_is_invariant_under_cache_pollution():
    x = sp.Symbol("x", real=True)
    formula = sp.And(x**2 - 2 < 0, x > 0)
    clear_algebraic_caches()
    clear_cad_caches()
    cold = is_satisfiable(formula, [x])
    warm = is_satisfiable(formula, [x])
    _pollute_caches()
    polluted = is_satisfiable(formula, [x])
    clear_algebraic_caches()
    clear_cad_caches()
    recomputed = is_satisfiable(formula, [x])
    assert cold == warm == polluted == recomputed is True


def test_cache_keys_preserve_assumption_distinct_symbol_identity():
    xr = sp.Symbol("x", real=True)
    xi = sp.Symbol("x", integer=True)
    assert xr != xi
    assert expr_key(xr + 1) != expr_key(xi + 1)


def test_cache_key_construction_is_structural_not_factorizing():
    x, y = sp.symbols("x y")
    expr = (x + y + 1) ** 12 - (x + y - 1) ** 12
    assert expr_key(expr) is expr
