"""Closure membership distinguishes exclusion from inability to approach."""

import sympy as sp

from semialg import point_in_closure

x, y = sp.symbols("x y", real=True)
punctured_disk = sp.And(x**2 + y**2 < 1, sp.Ne(x**2 + y**2, 0))
assert point_in_closure(punctured_disk, (0, 0), (x, y)) is True
assert point_in_closure(punctured_disk, (2, 0), (x, y)) is False
