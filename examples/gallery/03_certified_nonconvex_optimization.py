"""Executable companion to the semialg documentation example."""

import sympy as sp

from semialg import semialgebraic_minimize

x, y = sp.symbols("x y", real=True)
disk = x**2 + y**2 <= 1
objective = (x - 2) ** 2 + y**2

result = semialgebraic_minimize(objective, disk, [x, y], return_result=True)

print(result)
assert result.value == 1
assert result.points == ({x: 1, y: 0},)
assert result.attained is True
assert result.certified is True
