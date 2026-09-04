# Minkowski sums and support functions

Combine exact set operations with convex-geometric queries.

## Problem and interpretation

The Minkowski sum of \([0,1]\) and \([0,2]\) is \([0,3]\). For the unit disk,
the support function in direction \((3,4)\) is the Euclidean norm of the
direction, namely \(5\).

The two operations illustrate different implementation routes: Minkowski sum
is fundamentally an image/projection problem, while a support function is a
global optimization problem.

## Executable example

```python
import sympy as sp

from semialg import minkowski_sum, support_function

x, y = sp.symbols("x y", real=True)

left = (x >= 0) & (x <= 1)
right = (x >= 0) & (x <= 2)
summed = minkowski_sum(left, right, [x])

print(summed)
assert summed == (x >= 0) & (x <= 3)

disk = x**2 + y**2 <= 1
h = support_function(disk, [3, 4], [x, y])
print(h)
assert h == 5
```

The matching executable file is `examples/gallery/08_convex_geometry_operations.py`; the documentation
test suite runs every gallery script.
