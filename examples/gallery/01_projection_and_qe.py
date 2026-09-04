"""Executable companion to the semialg documentation example."""

import sympy as sp

from semialg import is_equal, semialgebraic_projection

x, y = sp.symbols("x y", real=True)
region = (y >= x**2) & (y <= 1)

projection = semialgebraic_projection(
    region,
    eliminate=[y],
    variables=[x, y],
)

print(projection)
assert is_equal(projection, (x >= -1) & (x <= 1), [x])
