# Distances and exact bounding boxes

Derive metric and coordinate bounds from semialgebraic sets.

## Problem and interpretation

Distance queries reduce to exact optimization of squared Euclidean distance.
Bounding boxes are coordinate-wise range problems. Keeping these operations as
first-class APIs avoids requiring users to formulate the corresponding
optimization problems by hand.

## Executable example

```python
import sympy as sp

from semialg import bounding_box, distance_between_regions

x, y = sp.symbols("x y", real=True)

left = (x >= -1) & (x <= 0)
right = (x >= 2) & (x <= 3)

distance = distance_between_regions(left, right, [x])
assert distance == 2

rectangle = (x**2 <= 4) & (y**2 <= 9)
box = bounding_box(rectangle, [x, y])

print(distance)
print(box)
assert box == {x: (-2, 2), y: (-3, 3)}
```

The matching executable file is `examples/gallery/09_distance_and_bounding_geometry.py`; the documentation
test suite runs every gallery script.
