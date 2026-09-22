import sympy as sp

from semialg import Box, Polytope
from semialg.topology.semialgebraic import (
    component_decomposition,
    dimension_strata,
    hardt_trivialization,
    triangulate_region,
)


def test_square_triangulation():
    square = Box(((0, 1), (0, 1)))
    triangulation = triangulate_region(square)
    assert triangulation.certified
    assert triangulation.dimension == 2
    assert triangulation.complex.simplex_count == 2
    assert sum(piece.measure() for piece in triangulation.complex.simplices()) == 1


def test_cube_pulling_triangulation():
    cube = Box(((0, 1), (0, 1), (0, 1)))
    triangulation = triangulate_region(cube)
    assert triangulation.dimension == 3
    assert triangulation.complex.simplex_count == 6
    assert sum(piece.measure() for piece in triangulation.complex.simplices()) == 1


def test_polytope_triangulation_uses_exact_vertices():
    polytope = Polytope(((0, 0), (2, 0), (2, 1), (0, 1)))
    triangulation = triangulate_region(polytope)
    assert set(triangulation.complex.vertices) == set(polytope.vertices)
    assert triangulation.complex.simplex_count == 2


def test_real_line_compact_triangulation():
    x = sp.symbols("x", real=True)
    region = sp.Or(sp.And(0 <= x, x <= 1), sp.And(2 <= x, x <= 3))
    triangulation = triangulate_region(region, (x,))
    assert triangulation.dimension == 1
    assert triangulation.complex.simplex_count == 2
    assert triangulation.complex.vertices == (
        (sp.Integer(0),),
        (sp.Integer(1),),
        (sp.Integer(2),),
        (sp.Integer(3),),
    )


def test_dimension_strata_mixed_dimensions():
    x, y = sp.symbols("x y", real=True)
    region = sp.Or(sp.And(0 < x, x < 1, 0 < y, y < 1), sp.Eq(x, 2) & sp.Eq(y, 0))
    decomposition = dimension_strata(region, (x, y))
    assert decomposition.dimension == 2
    assert {stratum.dimension for stratum in decomposition.strata} == {0, 2}


def test_component_decomposition_returns_samples():
    x = sp.symbols("x", real=True)
    region = sp.Or(sp.And(-2 <= x, x <= -1), sp.And(1 <= x, x <= 3))
    decomposition = component_decomposition(region, (x,))
    assert decomposition.count == 2
    samples = sorted(point[x] for point in decomposition.sample_points)
    assert samples[0] < 0 < samples[1]
    assert all(component.dimension == 1 for component in decomposition.components)


def test_hardt_quadratic_family():
    a, y = sp.symbols("a y", real=True)
    decomposition = hardt_trivialization(sp.Eq(y**2, a), (a,), y)
    assert decomposition.certified
    assert len(decomposition.strata) == 2
    by_components = {stratum.fiber_component_count: stratum for stratum in decomposition.strata}
    assert set(by_components) == {1, 2}
    assert sp.simplify(by_components[1].condition.subs(a, 0)) is sp.true
    assert sp.simplify(by_components[2].condition.subs(a, 1)) is sp.true
    assert sp.simplify(decomposition.image_condition.subs(a, -1)) is sp.false
    assert sp.simplify(decomposition.empty_fiber_condition.subs(a, -1)) is sp.true


def test_hardt_band_family_has_constant_fiber_dimension():
    a, y = sp.symbols("a y", real=True)
    region = sp.And(a > 0, -a < y, y < a)
    decomposition = hardt_trivialization(region, (a,), y)
    assert len(decomposition.strata) == 1
    stratum = decomposition.strata[0]
    assert stratum.fiber_dimension == 1
    assert stratum.fiber_component_count == 1
    piece = stratum.fiber_pieces[0]
    assert sp.simplify(piece.model_coordinate(0).subs(a, 2) - sp.Rational(1, 2)) == 0


def test_hardt_closed_band_cell_pieces_form_one_fiber_component():
    a, y = sp.symbols("a y", real=True)
    region = sp.And(a > 0, -a <= y, y <= a)
    decomposition = hardt_trivialization(region, (a,), y)
    assert len(decomposition.strata) == 1
    stratum = decomposition.strata[0]
    assert stratum.fiber_piece_count == 3
    assert stratum.fiber_component_count == 1
