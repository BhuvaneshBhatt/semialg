import sympy as sp

from semialg import SemialgebraicRegion, as_cad_region, simplify_region
from semialg.cad_algorithms.cells import extract_structured_cad_cells
from semialg.cad_region import triangulate_cad_cells


def test_exact_projection_image_preimage_methods():
    x, y = sp.symbols("x y", real=True)
    disk = SemialgebraicRegion(x**2 + y**2 <= 1, (x, y))
    proj = disk.project((y,))
    assert proj.variables == (x,)
    assert proj.equals_region(SemialgebraicRegion(sp.And(x >= -1, x <= 1), (x,)))

    interval = SemialgebraicRegion(sp.And(x >= 0, x <= 1), (x,))
    image = interval.image(x**2, variables=(y,))
    assert image.equals_region(SemialgebraicRegion(sp.And(y >= 0, y <= 1), (y,)))

    target = SemialgebraicRegion(sp.And(y >= -1, y <= 1), (y,))
    pre = target.preimage(x, (x,))
    assert pre.equals_region(SemialgebraicRegion(sp.And(x >= -1, x <= 1), (x,)))


def test_reusable_point_location_and_sign_vector():
    x = sp.symbols("x", real=True)
    region = SemialgebraicRegion(sp.And(x >= 0, x <= 2), (x,)).as_cad_region()
    located = region.locate_point((sp.Rational(1, 2),))
    assert located.selected is True
    assert located.dimension == 1
    assert located.cell_index
    assert located.sign_vector
    assert region.sign_vector((sp.Rational(1, 2),)) == located.sign_vector


def test_cell_complex_and_euler_characteristic():
    x = sp.symbols("x", real=True)
    region = SemialgebraicRegion(sp.And(x >= 0, x <= 1), (x,))
    complex_ = region.cell_complex()
    assert complex_.f_vector == (2, 1)
    assert complex_.euler_characteristic() == 1
    assert region.euler_characteristic() == 1
    sector = next(i for i, c in enumerate(complex_.cells) if c.dimension == 1)
    assert len(complex_.boundary_of(sector)) == 2
    assert complex_.incidence_matrix(1).shape == (2, 1)


def test_singular_regular_locus_and_local_dimension():
    x, y = sp.symbols("x y", real=True)
    cusp_region = SemialgebraicRegion(y**2 <= x**3, (x, y))
    singular = cusp_region.singular_locus()
    assert singular.contains((0, 0))
    assert not singular.contains((1, 0))
    assert cusp_region.regular_locus().contains((1, 0))

    segment = SemialgebraicRegion(sp.And(x >= 0, x <= 1), (x,))
    assert segment.local_dimension((0,)) == 1
    assert segment.local_dimension((sp.Rational(1, 2),)) == 1
    assert segment.local_dimension((2,)) == -1


def test_conforming_adjacency_mesh_reports_shared_faces():
    x, y = sp.symbols("x y", real=True)
    formula = sp.And(
        x >= 0,
        x <= 2,
        y >= 0,
        y <= 1,
        sp.Or(x <= 1, x >= 1, evaluate=False),
        evaluate=False,
    )
    region = SemialgebraicRegion(formula, (x, y)).as_cad_region()
    decomp = extract_structured_cad_cells(region.result)
    cells = tuple(c for c in decomp.full_dimensional_cells if c.bounded)
    assert len(cells) >= 2
    mesh = triangulate_cad_cells(cells, slices=2, require_conforming=True)
    assert mesh.conforming is True
    assert mesh.shared_vertex_count > 0


def test_region_introspection_and_exact_simplification():
    x, a = sp.symbols("x a", real=True)
    region = SemialgebraicRegion(sp.And(x >= 0, x <= a), (x,))
    assert region.region_variables("coordinates") == (x,)
    assert region.parameters == (a,)
    assert set(region.region_variables("all")) == {x, a}

    nested = SemialgebraicRegion(sp.Or(sp.And(x >= 0, x <= 1), sp.And(x >= 0, x <= 2)), (x,))
    simplified = simplify_region(nested, exact=True)
    assert simplified.equals_region(SemialgebraicRegion(sp.And(x >= 0, x <= 2), (x,)))


def test_exact_region_measure_and_integration_hooks():
    x, y = sp.symbols("x y", real=True)
    box = SemialgebraicRegion(sp.And(x >= 0, x <= 1, y >= 0, y <= 2), (x, y))
    assert sp.simplify(box.measure() - 2) == 0
    assert sp.simplify(box.integrate(x + y) - 3) == 0

    cad_region = as_cad_region(box)
    pieces = cad_region.integrate(x + y, evaluate=False)
    assert pieces
    assert all(piece.certified_bounds for piece in pieces)
