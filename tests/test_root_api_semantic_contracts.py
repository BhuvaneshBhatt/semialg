"""Independent semantic hardening for the final root-API adequacy gaps."""

from __future__ import annotations

import sympy as sp

import semialg as sa


def _triangle_area(a, b, c):
    return abs(sp.Matrix([[b[0] - a[0], b[1] - a[1]], [c[0] - a[0], c[1] - a[1]]]).det()) / 2


def _tetra_volume(points):
    a, b, c, d = map(sp.Matrix, points)
    return abs(sp.Matrix.hstack(b - a, c - a, d - a).det()) / 6


def test_prove_positive_positive_scaling_metamorphic_invariant():
    x = sp.symbols("x", real=True)
    base = x**2 + 1
    assert sa.prove_positive(base, (x,)) is True
    assert sa.prove_positive(7 * base, (x,)) is True
    # A negative scale reverses the relation and must not inherit positivity.
    assert sa.prove_positive(-7 * base, (x,)) is False


def test_zariski_closure_parabola_matches_elimination_oracle():
    t, x, y = sp.symbols("t x y", real=True)
    closure = sa.zariski_closure((t, t**2), (t,), image_variables=(x, y))
    # Independent elimination by substitution gives y - x**2 = 0.
    assert sa.equivalent(closure, sp.Eq(y - x**2, 0), (x, y)) is True


def test_function_domain_equivalent_rational_presentations_agree():
    x = sp.symbols("x", real=True)
    expanded = 1 / (x**2 - 1)
    factored = 1 / ((x - 1) * (x + 1))
    left = sa.function_domain(expanded, (x,))
    right = sa.function_domain(factored, (x,))
    assert sa.equivalent(left, right, (x,)) is True
    assert sa.is_satisfiable(sp.And(left, sp.Eq(x, 1)), (x,)) is False
    assert sa.is_satisfiable(sp.And(left, sp.Eq(x, -1)), (x,)) is False


def test_region_moment_translation_covariance_closed_form_oracle():
    x = sp.symbols("x", real=True)
    # Integral_[a,a+2] x dx = 2a + 2, here a=3.
    shifted = sa.region_moment(sp.And(x >= 3, x <= 5), (x,), powers=(1,))
    assert sp.simplify(shifted - 8) == 0
    base_mass = sa.region_moment(sp.And(x >= 0, x <= 2), (x,))
    base_first = sa.region_moment(sp.And(x >= 0, x <= 2), (x,), powers=(1,))
    assert sp.simplify(shifted - (base_first + 3 * base_mass)) == 0


def test_distance_to_region_translation_metamorphic_invariant():
    x = sp.symbols("x", real=True)
    d1 = sa.distance_to_region((5,), sp.And(x >= 0, x <= 2), (x,))
    d2 = sa.distance_to_region((12,), sp.And(x >= 7, x <= 9), (x,))
    assert d1 == d2 == 3
    assert sa.distance_to_region((1,), sp.And(x >= 0, x <= 2), (x,)) == 0


def test_polygonal_region_fill_orientation_representation_contract():
    square = [(0, 0), (2, 0), (2, 2), (0, 2)]
    positive = sa.polygonal_region_from_paths([square], fill_rule="positive")
    reversed_positive = sa.polygonal_region_from_paths(
        [list(reversed(square))], fill_rule="positive"
    )
    reversed_negative = sa.polygonal_region_from_paths(
        [list(reversed(square))], fill_rule="negative"
    )
    assert len(positive.components) == 1
    assert len(reversed_positive.components) == 0
    assert len(reversed_negative.components) == 1


def test_tetrahedralize_cell_preserves_cube_volume_oracle():
    cube = [(i, j, k) for k in (0, 1) for j in (0, 1) for i in (0, 1)]
    result = sa.tetrahedralize_cell(cube)
    total = sum(_tetra_volume([cube[i] for i in tet]) for tet in result.tetrahedra)
    assert sp.simplify(total - 1) == 0
    assert all(_tetra_volume([cube[i] for i in tet]) > 0 for tet in result.tetrahedra)


def test_topology_summary_circle_known_invariants_oracle():
    x, y = sp.symbols("x y", real=True)
    summary = sa.topology_summary(sp.Eq(x**2 + y**2, 1), (x, y))
    assert summary.dimension == 1
    assert summary.connected_components == 1
    assert summary.euler_characteristic == 0
    assert summary.betti_numbers[:2] == (1, 1)


def test_triangulate_polytope_preserves_square_area_oracle():
    vertices = [(0, 0), (2, 0), (2, 1), (0, 1)]
    decomposition = sa.triangulate_polytope(sa.Polytope(vertices), strategy="pulling")
    area = sum(
        _triangle_area(*(decomposition.vertices[i] for i in simplex))
        for simplex in decomposition.simplices
    )
    assert sp.simplify(area - 2) == 0
    # Every simplex is nondegenerate and uses valid decomposition vertices.
    assert all(
        _triangle_area(*(decomposition.vertices[i] for i in s)) > 0 for s in decomposition.simplices
    )


def test_discretize_solution_unbounded_interval_boundary_clipping():
    x = sp.symbols("x", real=True)
    solution = sa.solve_semialgebraic(x > 1, (x,), count=0)
    data = sa.discretize_solution(solution, bounds=((-4, 6),))
    assert data.dimension == 1
    assert data.source == "interval-components"
    assert data.segments == (((sp.Integer(1),), (sp.Integer(6),)),)


def test_sign_vector_positive_scaling_and_order_metamorphic_contract():
    x = sp.symbols("x", real=True)
    polys = (x - 1, x, x + 1)
    point = {x: sp.Rational(1, 2)}
    signs = sa.sign_vector(polys, point)
    scaled = sa.sign_vector(tuple(5 * p for p in polys), point)
    reversed_signs = sa.sign_vector(tuple(reversed(polys)), point)
    assert signs == scaled == (-1, 1, 1)
    assert reversed_signs == tuple(reversed(signs))
