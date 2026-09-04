"""Executable companion to the semialg documentation example."""

import sympy as sp

from semialg import argmin_set, is_equal

x, y = sp.symbols("x y", real=True)
box = (x >= -1) & (x <= 1) & (y >= 0) & (y <= 1)

minimizers = argmin_set(x**2, box, [x, y])
expected = sp.Eq(x, 0) & (y >= 0) & (y <= 1)

print(minimizers)
assert is_equal(minimizers, expected, [x, y])
