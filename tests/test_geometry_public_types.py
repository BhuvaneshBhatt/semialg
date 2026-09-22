"""Construction contracts for canonical geometry public types."""

import semialg


def test_canonical_geometry_types_construct():
    assert semialg.Box(((0, 1), (0, 1))).dimension() == 2
    assert semialg.FinitePointSet([(0, 0), (1, 1)]).dimension() == 0
    assert semialg.Interval(0, 1).dimension() == 1
    assert semialg.Parallelogram((0, 0), ((1, 0), (0, 1))).dimension() == 2
    assert (
        semialg.TetrahedralComplex([[(0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1)]]).dimension() == 3
    )
    assert semialg.Zonotope((0, 0), ((1, 0), (0, 1))).dimension() == 2
