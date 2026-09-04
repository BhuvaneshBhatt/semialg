# Exact measure, centroid, covariance, and inertia

Compute several geometric statistics of the unit disk from exact region integrals.

## Problem and interpretation

Symmetry forces the centroid to the origin and the off-diagonal moments to
vanish. The uniform probability covariance is

\[
\frac14 I,
\]

while the unit-density planar inertia tensor about the origin is

\[
\frac{\pi}{4}I.
\]

The example demonstrates how measure and moment APIs share the same exact
integration machinery but return mathematically different normalizations.

## Executable example

```python
import sympy as sp

from semialg import (
    centroid,
    covariance_matrix,
    inertia_tensor,
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
```

The matching executable file is `examples/gallery/06_measure_centroid_and_moments.py`; the documentation
test suite runs every gallery script.
