"""Executable companion to the semialg documentation example."""

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
