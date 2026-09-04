import sympy as sp

from semialg import (
    SemialgebraicRegion,
    as_cad_region,
    is_regular_closed_region,
    is_regular_open_region,
)
from semialg.cad_algorithms.structured_cells import extract_structured_cad_cells
from semialg.cad_region import (
    cad_combine,
    cad_extend,
    cad_signature,
    evaluate_delineable_curve,
    evaluate_delineable_surfaces,
    triangulate_cad_cell,
)


def test_cad_signature_is_stable_and_informative():
    x = sp.symbols("x", real=True)
    region = as_cad_region(SemialgebraicRegion(sp.And(x >= 0, x <= 1), (x,)))
    sig = cad_signature(region)
    assert sig.ambient_dimension == 1
    assert sig.region_dimension == 1
    assert sig.selected_cell_count == 3
    assert sig.total_leaf_count >= sig.selected_cell_count
    assert sig == region.signature()


def test_cad_boolean_combine_reuses_shared_decomposition():
    x = sp.symbols("x", real=True)
    region = as_cad_region(SemialgebraicRegion(sp.And(x >= 0, x <= 1), (x,)))
    combined = cad_combine(region, region, op="intersection")
    assert combined.result.cad is region.result.cad
    assert combined.result.diagnostics["cad_reused"] is True
    assert combined.result.cell_set.indices == region.result.cell_set.indices


def test_cad_extend_reuses_only_existing_sign_invariant_boundaries():
    x = sp.symbols("x", real=True)
    region = as_cad_region(SemialgebraicRegion(sp.And(x >= 0, x <= 1), (x,)))
    reused = cad_extend(region, x >= 0)
    refined = cad_extend(region, x < sp.Rational(1, 2))
    assert reused.result.cad is region.result.cad
    assert reused.result.diagnostics["cad_reused"] is True
    assert refined.result.diagnostics["cad_reused"] is False


def test_regular_open_and_closed_region_operations():
    x = sp.symbols("x", real=True)
    closed = SemialgebraicRegion(sp.And(x >= 0, x <= 1), (x,))
    opened = SemialgebraicRegion(sp.And(x > 0, x < 1), (x,))
    point = SemialgebraicRegion(sp.Eq(x, 0), (x,))
    assert closed.is_regular_closed()
    assert opened.is_regular_open()
    assert is_regular_closed_region(closed)
    assert is_regular_open_region(opened)
    assert not point.is_regular_closed()
    assert closed.interior_closure().equals_region(closed)
    assert opened.closure_interior().equals_region(opened)


def test_triangulate_bounded_structured_cad_cell():
    x, y = sp.symbols("x y", real=True)
    expr = sp.And(x > 0, x < 1, y > 0, y < x + 1)
    decomp = extract_structured_cad_cells(expr, (x, y))
    cell = decomp.full_dimensional_cells[0]
    mesh = triangulate_cad_cell(cell, slices=4)
    assert len(mesh.coordinates) == 10
    assert len(mesh.simplices) == 8
    assert all(len(simplex) == 3 for simplex in mesh.simplices)


def test_evaluate_delineable_curve_root_indices():
    x, y, z = sp.symbols("x y z", real=True)
    points = evaluate_delineable_curve(
        ((z**2 - y, 1),),
        (y**2 - x, 2),
        (x, y, z),
        (1, 4),
    )
    assert points[0] == (1.0, 1.0, -1.0)
    assert points[1][0:2] == (4.0, 2.0)
    assert abs(points[1][2] + 2**0.5) < 1e-12


def test_evaluate_delineable_surfaces_and_merge():
    x, y, z = sp.symbols("x y z", real=True)
    coords, inds = evaluate_delineable_surfaces(
        ((z**2 - x - y, 2), (z - sp.sqrt(x + y), 1)),
        (x, y, z),
        ((1, 0), (4, 0)),
        merge_tolerance=1e-10,
    )
    assert len(coords) == 2
    assert inds[0] == inds[1]
    assert coords == ((1.0, 0.0, 1.0), (4.0, 0.0, 2.0))


def test_structured_cell_extracts_delineable_boundary_pairs():
    x, y = sp.symbols("x y", real=True)
    expr = sp.And(x > 0, x < 1, y**2 < x)
    decomp = extract_structured_cad_cells(expr, (x, y))
    cell = decomp.full_dimensional_cells[0]
    pairs = cell.boundary_pairs(variable=y)
    assert len(pairs) == 2
    assert {index for _, index in pairs} == {1, 2}
    assert all(
        sp.expand(poly - (x - y**2)) == 0 or sp.expand(poly - (y**2 - x)) == 0 for poly, _ in pairs
    )
    assert all(desc.fiber_variable == y for desc in cell.boundary_descriptors if desc.level == 2)


def test_higher_dimensional_cad_cell_tetrahedral_meshing():
    x, y, z = sp.symbols("x y z", real=True)
    expr = sp.And(x > 0, x < 1, y > 0, y < 1, z > 0, z < 1)
    decomp = extract_structured_cad_cells(expr, (x, y, z))
    cell = decomp.full_dimensional_cells[0]
    mesh = triangulate_cad_cell(cell, slices=1, inset=0)
    assert mesh.dimension == 3
    assert len(mesh.coordinates) == 8
    assert len(mesh.simplices) == 6
    assert all(len(simplex) == 4 for simplex in mesh.simplices)


def test_adjacency_aware_meshing_shares_common_boundary_vertices():
    from semialg.cad_region import triangulate_cad_cells

    x, y = sp.symbols("x y", real=True)
    expr = sp.And(x > 0, x < 2, sp.Ne(x, 1), y > 0, y < 1)
    decomp = extract_structured_cad_cells(expr, (x, y))
    cells = decomp.full_dimensional_cells
    assert len(cells) == 2
    local = [triangulate_cad_cell(cell, slices=2, inset=0) for cell in cells]
    mesh = triangulate_cad_cells(cells, slices=2, tolerance=1e-12)
    assert len(mesh.coordinates) < sum(len(item.coordinates) for item in local)
    shared = [owners for owners in mesh.vertex_cells if len(owners) > 1]
    assert len(shared) >= 2
    assert len(mesh.simplex_cells) == len(mesh.simplices)
