# Connected components and Euler characteristic

Use CAD semantics to recover components and basic exact topology.

## Problem and interpretation

The two-ray set has two connected components. The disk has Euler characteristic
\(1\), whereas the closed annulus has Euler characteristic \(0\).

These operations use exact CAD-cell topology rather than a sampled graph.
They are therefore insensitive to numerical meshing resolution, though the
cost can be much higher than numerical topology for large problems.

## Executable example

```python
import sympy as sp

from semialg import connected_components, euler_characteristic, is_interior_disjoint

x, y = sp.symbols("x y", real=True)

two_rays = (x <= -1) | (x >= 1)
components = connected_components(two_rays, [x])
print(components)
assert len(components) == 2
assert is_interior_disjoint((x >= 0) & (x <= 1), (x >= 1) & (x <= 2), [x])

disk = x**2 + y**2 <= 1
annulus = (x**2 + y**2 >= 1) & (x**2 + y**2 <= 4)

assert euler_characteristic(disk, [x, y]) == 1
assert euler_characteristic(annulus, [x, y]) == 0
```

The matching executable file is `examples/gallery/05_topology_and_components.py`; the documentation
test suite runs every gallery script.
