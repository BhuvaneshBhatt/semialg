import sympy as sp

from semialg import singular_locus, tangent_cone, tangent_space


def test_singular_locus_of_cusp_contains_only_origin_over_reals():
    x, y = sp.symbols("x y", real=True)
    locus = singular_locus(y**2 - x**3, [x, y])
    assert sp.simplify(locus.subs({x: 0, y: 0})) is sp.true
    assert sp.simplify(locus.subs({x: 1, y: 1})) is sp.false


def test_tangent_space_at_cusp_origin_is_full_zariski_tangent_plane():
    x, y = sp.symbols("x y", real=True)
    result = tangent_space(y**2 - x**3, {x: 0, y: 0}, [x, y])
    assert result.dimension == 2
    assert result.jacobian == sp.zeros(1, 2)


def test_tangent_space_at_regular_parabola_point():
    x, y = sp.symbols("x y", real=True)
    result = tangent_space(y - x**2, {x: 0, y: 0}, [x, y])
    assert result.dimension == 1
    assert result.jacobian == sp.Matrix([[0, 1]])


def test_tangent_cone_of_cusp_is_double_x_axis():
    x, y = sp.symbols("x y", real=True)
    result = tangent_cone(y**2 - x**3, {x: 0, y: 0}, [x, y])
    dx, dy = result.direction_variables
    assert result.certified is True
    assert sp.expand(result.initial_forms[0] - dy**2) == 0
    assert dx in result.direction_variables


def test_tangent_cone_is_exact_for_multigenerator_cancellation():
    x, y = sp.symbols("x y", real=True)
    result = tangent_cone(
        (x**2 + y**3, x**2 - y**3),
        {x: 0, y: 0},
        [x, y],
    )
    dx, dy = result.direction_variables
    basis = sp.groebner(result.ideal_generators, dx, dy, order="lex")
    assert basis.reduce(dx**2)[1] == 0
    assert basis.reduce(dy**3)[1] == 0
    assert result.certified is True
    assert result.method == "saturated_m_adic_deformation"


def test_tangent_cone_depends_on_ideal_not_supplied_generators():
    x, y = sp.symbols("x y", real=True)
    first = tangent_cone(
        (x**2 + y**3, x**2 - y**3),
        {x: 0, y: 0},
        [x, y],
    )
    second = tangent_cone(
        (x**2 + y**3, 2 * x**2, 2 * y**3),
        {x: 0, y: 0},
        [x, y],
    )
    assert first.ideal_generators == second.ideal_generators
