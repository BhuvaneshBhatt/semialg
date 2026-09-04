"""Executable companion to the semialg documentation example."""

import sympy as sp

from semialg import connected_components, euler_characteristic

x, y = sp.symbols("x y", real=True)

two_rays = (x <= -1) | (x >= 1)
components = connected_components(two_rays, [x])
print(components)
assert len(components) == 2

disk = x**2 + y**2 <= 1
annulus = (x**2 + y**2 >= 1) & (x**2 + y**2 <= 4)

assert euler_characteristic(disk, [x, y]) == 1
assert euler_characteristic(annulus, [x, y]) == 0
