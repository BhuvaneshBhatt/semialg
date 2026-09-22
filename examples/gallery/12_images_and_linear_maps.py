"""Executable companion to the semialg documentation example."""

import sympy as sp

from semialg import contains_point, region_image

x, y, u, v = sp.symbols("x y u v", real=True)
disk = x**2 + y**2 <= 1

image = region_image(
    disk,
    [x + y, x - y],
    variables=[x, y],
    image_variables=[u, v],
)

print(image)
assert contains_point(image, {u: 0, v: 0}, [u, v])
assert contains_point(image, {u: 1, v: 1}, [u, v])
assert not contains_point(image, {u: 2, v: 0}, [u, v])
