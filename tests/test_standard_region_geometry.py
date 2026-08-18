import sympy as sp

from semialg import BallRegion, BoxRegion, PolygonRegion, discretize_region_geometry


def test_discretize_explicit_standard_regions():
    x, y = sp.symbols("x y", real=True)
    data = discretize_region_geometry(BoxRegion([(0, 1), (0, 2)]), variables=[x, y])
    assert data.dimension == 2
    assert len(data.polygons) == 1

    poly = discretize_region_geometry(PolygonRegion([(0, 0), (1, 0), (0, 1)]), variables=[x, y])
    assert poly.source == "standard_region:PolygonRegion"
    assert len(poly.polygons[0]) == 3

    disk = discretize_region_geometry(BallRegion((0, 0), 1), variables=[x, y], samples_per_curve=8)
    assert len(disk.polygons) == 1
    assert len(disk.polygons[0]) == 8
