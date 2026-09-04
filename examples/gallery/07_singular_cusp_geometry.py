"""Executable companion to the semialg documentation example."""

import sympy as sp

from semialg import singular_locus, tangent_cone, tangent_space

x, y = sp.symbols("x y", real=True)
f = y**2 - x**3

singular = singular_locus([f], [x, y])
space = tangent_space([f], {x: 0, y: 0}, [x, y])
cone = tangent_cone([f], {x: 0, y: 0}, [x, y])

print(singular)
print(space)
print(cone)

assert space.dimension == 2
assert cone.certified is True
assert cone.ideal_generators == (cone.direction_variables[1] ** 2,)
