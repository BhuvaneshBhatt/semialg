"""Independent behavioral smoke coverage for APIs added during unified milestones."""

import sympy as sp

import semialg


def test_algebraic_statistics_secondary_calls():
    x, y, t = sp.symbols("x y t", real=True)
    assert semialg.certified_radicalization((x**2,), (x,))
    assert semialg.irreducible_components((x * y,), (x, y))
    assert semialg.reduced_component_singular_loci((x * y,), (x, y)) is not None
    assert semialg.minimal_prime_intersections((x * y,), (x, y)) is not None
    assert semialg.local_dimension_strata((x * y,), (x, y)) is not None
    assert semialg.stratified_singular_geometry((x * y,), (x, y)) is not None
    assert semialg.local_branch_geometry((x * y,), {x: 0, y: 0}, (x, y)) is not None
    assert semialg.implicitize_polynomial_map((t, t**2), (t,), image_variables=(x, y))
    assert semialg.zariski_closure((t, t**2), (t,), image_variables=(x, y)) is not None


def test_semialgebraic_constraint_secondary_calls():
    x, y = sp.symbols("x y", real=True)
    region = sp.And(x >= 0, y >= 0, x + y <= 1)
    assert semialg.polynomial_constraints(region, (x, y)).clauses
    assert semialg.active_constraints(region, (0, 0), (x, y)).constraints
    assert semialg.relative_interior(region, (x, y)) is not None
    assert semialg.relative_boundary(region, (x, y)) is not None
    assert semialg.semialgebraic_tangent_cone(x >= 0, (0,), (x,)) is not None
    cert = semialg.nonnegative_combination_certificate(2 * x, (x,), (x,), allow_sos=False)
    assert cert is not None and semialg.verify_nonnegative_combination_certificate(cert, (x,))
    assert semialg.implied_polynomial_inequality(region, x >= 0, (x, y)).certified
    assert semialg.redundant_polynomial_inequalities(sp.And(x >= 0, 2 * x >= 0), (x,))
    assert semialg.component_constraint_descriptions(sp.And(sp.Eq(x * y, 0), x >= 0), (x, y))


def test_parameterization_and_qe_secondary_calls():
    t, x, y = sp.symbols("t x y", real=True)
    assert semialg.parameterization_geometry((t**2, t**3), sp.true, (t,)) is not None
    assert semialg.parameterization_critical_locus((t**2, t**3), sp.true, (t,)) is not None
    assert semialg.parameterization_critical_values((t**2, t**3), sp.true, (t,)) is not None
    assert semialg.quantifier_eliminate(semialg.Exists(x, sp.Eq(x, y))) is sp.true
    assert semialg.project_region(sp.And(x >= y, x <= 1), (x,)) is not None


def test_restored_geometry_secondary_calls():
    component = semialg.PolygonalComponent(
        ((0, 0), (0, 4), (4, 4), (4, 0)),
        (((1, 1), (2, 1), (2, 2), (1, 2)),),
    )
    region = semialg.PolygonalSet((component,))
    assert semialg.polygon_vertices(region)
    assert semialg.outer_polygons(region)
    assert semialg.inner_polygons(region)
    dedup = semialg.deduplicate_indexed_vertices(((0, 0), (1, 0), (0, 0)), ((0, 1), (2, 1)))
    assert len(dedup.vertices) == 2
    cube = semialg.Hexahedron(
        ((0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0), (0, 0, 1), (1, 0, 1), (1, 1, 1), (0, 1, 1))
    )
    boundary = semialg.convert_region(cube, "boundary")
    assert semialg.polyhedron_vertices(boundary)
    assert semialg.polyhedron_face_indices(boundary)
    assert semialg.outer_polyhedra(boundary)
    assert semialg.inner_polyhedra(boundary) == ()
    assert semialg.canonicalize_region([(0, 0), (1, 0), (0, 1)]) is not None
    assert semialg.polygonal_region_from_paths([((0, 0), (1, 0), (1, 1), (0, 1))]) is not None
