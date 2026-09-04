"""Executable companion to the semialg documentation example."""

import sympy as sp

from semialg import function_range

x, y = sp.symbols("x y", real=True)
disk = x**2 + y**2 <= 1

value_set = function_range(x + y, disk, [x, y])

print(value_set)
t = next(iter(value_set.free_symbols - {x, y}))
assert sp.simplify(value_set.subs(t, sp.sqrt(2))) is sp.true
assert sp.simplify(value_set.subs(t, 2)) is sp.false
