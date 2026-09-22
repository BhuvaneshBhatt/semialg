"""Executable companion to the semialg documentation example."""

import sympy as sp

from semialg import (
    centroid,
    covariance_matrix,
    inertia_tensor,
    region_measure,
    semialgebraic_measure,
)

x, y = sp.symbols("x y", real=True)
disk = x**2 + y**2 <= 1

area = semialgebraic_measure(disk, [x, y])
center = centroid(disk, [x, y])
cov = covariance_matrix(disk, [x, y])
inertia = inertia_tensor(disk, [x, y])

print(area, center)
print(cov)
print(inertia)

assert area == sp.pi
assert center == {x: 0, y: 0}
assert cov == sp.eye(2) / 4
assert inertia == sp.eye(2) * sp.pi / 4

segment = sp.Eq(y, 0) & (x >= 0) & (x <= 1)
assert region_measure(segment, [x, y]) == 0
assert region_measure(segment, [x, y], measure_dimension="intrinsic") == 1
