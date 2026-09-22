import sympy as sp

from semialg import implicitize_polynomial_map, singular_locus, zariski_closure


def test_implicitize_parabola_polynomial_map():
    t, x, y = sp.symbols("t x y", real=True)
    equations = implicitize_polynomial_map((t, t**2), (t,), image_variables=(x, y))
    assert len(equations) == 1
    assert sp.Poly(equations[0], x, y).monic() == sp.Poly(x**2 - y, x, y).monic()


def test_implicitize_twisted_cubic():
    t = sp.symbols("t", real=True)
    x, y, z = sp.symbols("x y z", real=True)
    result = implicitize_polynomial_map(
        (t, t**2, t**3), (t,), image_variables=(x, y, z), return_result=True
    )
    basis = sp.groebner(result.equations, x, y, z, order="lex")
    for invariant in (x**2 - y, x * y - z, y**2 - x * z):
        assert sp.expand(basis.reduce(invariant)[1]) == 0
    assert not any(t in equation.free_symbols for equation in result.equations)


def test_implicitization_respects_algebraic_parameter_domain():
    s, t = sp.symbols("s t", real=True)
    x, y = sp.symbols("x y", real=True)
    equations = implicitize_polynomial_map(
        (s, t), (s, t), image_variables=(x, y), domain_equations=(s * t - 1,)
    )
    basis = sp.groebner(equations, x, y, order="lex")
    assert sp.expand(basis.reduce(x * y - 1)[1]) == 0


def test_zariski_closure_returns_exact_formula_and_metadata():
    t, x, y = sp.symbols("t x y", real=True)
    result = zariski_closure((t, t**2), (t,), image_variables=(x, y), return_result=True)
    assert result.variables == (x, y)
    assert sp.simplify(result.formula.subs({x: 2, y: 4})) is sp.true
    assert sp.simplify(result.formula.subs({x: 2, y: 5})) is sp.false


def test_dominant_polynomial_map_has_full_ambient_zariski_closure():
    s, t = sp.symbols("s t", real=True)
    x, y = sp.symbols("x y", real=True)
    assert zariski_closure((s, t), (s, t), image_variables=(x, y)) is sp.true


def test_singular_locus_of_implicitized_cusp():
    t, x, y = sp.symbols("t x y", real=True)
    equations = implicitize_polynomial_map((t**2, t**3), (t,), image_variables=(x, y))
    locus = singular_locus(equations, (x, y))
    assert sp.simplify(locus.subs({x: 0, y: 0})) is sp.true
    assert sp.simplify(locus.subs({x: 1, y: 1})) is sp.false


def test_nonpolynomial_mapping_is_rejected():
    t = sp.symbols("t", real=True)
    try:
        implicitize_polynomial_map((sp.sin(t),), (t,))
    except ValueError as exc:
        assert "polynomial" in str(exc)
    else:
        raise AssertionError("expected non-polynomial map rejection")
