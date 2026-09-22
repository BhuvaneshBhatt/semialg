import semialg
from semialg import Box, FinitePointSet, Interval, Parallelogram, TetrahedralComplex


def test_canonical_shape_names_are_runtime_class_names():
    for cls in (FinitePointSet, Interval, Box, Parallelogram, TetrahedralComplex):
        assert cls.__name__ == cls.__qualname__
        assert not cls.__name__.endswith("Region")


def test_removed_shape_region_names_are_not_public():
    removed = (
        "Point" + "Region",
        "Interval" + "Region",
        "Box" + "Region",
        "Simplex" + "Region",
        "Polygon" + "Region",
        "Tetrahedron" + "Region",
        "Polyhedron" + "Region",
        "Parallelogram" + "Region",
        "Parallelepiped" + "Region",
        "Prism" + "Region",
        "Pyramid" + "Region",
        "Ball" + "Region",
        "Sphere" + "Region",
        "SphericalShell" + "Region",
        "Cylinder" + "Region",
        "Cone" + "Region",
        "Stadium" + "Region",
        "Capsule" + "Region",
        "Conic" + "Region",
    )
    assert all(not hasattr(semialg, name) for name in removed)
