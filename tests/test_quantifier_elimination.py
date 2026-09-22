import sympy as sp

from semialg import Exists, ForAll, project_region, quantifier_eliminate


def test_affine_existential_presolve():
    x, y = sp.symbols("x y", real=True)
    result = quantifier_eliminate(Exists(x, sp.And(sp.Eq(x, y + 1), x > 0)), return_result=True)
    assert result.method == "affine-presolve"
    assert sp.simplify(sp.Equivalent(result.formula, y > -1)) is sp.true


def test_linear_projection_uses_fourier_motzkin():
    x, y = sp.symbols("x y", real=True)
    result = project_region(sp.And(x >= y, x <= 1), (x,), return_result=True)
    assert result.method == "fourier-motzkin"
    assert sp.simplify(sp.Equivalent(result.formula, y <= 1)) is sp.true


def test_quadratic_existential_projection():
    x, a = sp.symbols("x a", real=True)
    result = quantifier_eliminate(Exists(x, sp.Eq(x**2, a)), return_result=True)
    # Specialist portfolio may use VS or CAD depending on planner details.
    assert x not in result.formula.free_symbols
    assert sp.simplify(sp.Equivalent(result.formula, a >= 0)) is sp.true
    assert result.certified


def test_universal_quantifier():
    x, a = sp.symbols("x a", real=True)
    out = quantifier_eliminate(ForAll(x, x**2 + a >= 0))
    assert x not in out.free_symbols
    for value, expected in [(-1, False), (0, True), (2, True)]:
        assert bool(out.subs(a, value)) is expected


def test_explicit_prefix_and_nested_prefix_are_equivalent():
    x, y, a = sp.symbols("x y a", real=True)
    matrix = sp.Implies(x >= a, y**2 >= 0)
    explicit = quantifier_eliminate(matrix, (("exists", x), ("forall", y)))
    nested = quantifier_eliminate(Exists(x, ForAll(y, matrix)))
    assert sp.simplify(sp.Equivalent(explicit, nested)) is sp.true


def test_projection_of_hyperbola_keeps_real_nonzero_condition():
    x, y = sp.symbols("x y", real=True)
    out = project_region(sp.Eq(x * y, 1), (y,))
    assert y not in out.free_symbols
    for value, expected in [(0, False), (1, True), (-2, True)]:
        assert bool(out.subs(x, value)) is expected


def test_independent_incidence_blocks_are_eliminated_separately():
    x, y, a, b = sp.symbols("x y a b", real=True)
    result = quantifier_eliminate(
        sp.And(x**2 >= a, y**2 >= b),
        (("exists", x), ("exists", y)),
        return_result=True,
    )
    assert result.method == "independent-blocks"
    assert len(result.variable_blocks) >= 2
    assert not ({x, y} & result.formula.free_symbols)
    assert result.formula is sp.true or sp.simplify(result.formula) is sp.true


def test_groebner_consequences_enable_safe_affine_presolve():
    x, y, a = sp.symbols("x y a", real=True)
    # Subtracting the two equations yields x-y=0; neither input equality is
    # itself globally affine-solvable for x because of the quadratic term.
    matrix = sp.And(sp.Eq(x**2 + x + a, 0), sp.Eq(x**2 + y + a, 0))
    result = quantifier_eliminate(matrix, (("exists", x),), variables=(a, y, x), return_result=True)
    assert x not in result.formula.free_symbols
    # The real projection is y^2 + y + a == 0, not merely a Zariski closure
    # guessed from an elimination ideal.
    expected = sp.Eq(y**2 + y + a, 0)
    for av, yv in [(0, 0), (0, -1), (1, 0), (-2, 1)]:
        assert bool(result.formula.subs({a: av, y: yv})) == bool(expected.subs({a: av, y: yv}))
    assert result.metadata.get("groebner_equalities_added") or "groebner" in result.method
