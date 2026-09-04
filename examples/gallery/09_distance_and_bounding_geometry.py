"""Executable companion to the semialg documentation example."""

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
