# Polynomial images and transformed regions

Compute the exact image of a disk under an invertible linear map.

## Problem and interpretation

The linear map

\[
(u,v)=(x+y,x-y)
\]

scales Euclidean norm by \(\sqrt2\), so the unit disk maps to

\[
u^2+v^2\le2.
\]

The raw reconstructed image formula may contain exact algebraic root-function
boundaries rather than printing in this simplest form. The executable checks
representative interior, boundary, and exterior points. Mathematically the image
is the disk `u**2 + v**2 <= 2`; in general, exact set equality should be checked
semantically rather than by comparing printed formulas.

## Executable example

```python
import sympy as sp

from semialg import contains_point, semialgebraic_image

x, y, u, v = sp.symbols("x y u v", real=True)
disk = x**2 + y**2 <= 1

image = semialgebraic_image(
    [x + y, x - y],
    disk,
    variables=[x, y],
    image_variables=[u, v],
)

print(image)
assert contains_point(image, {u: 0, v: 0}, [u, v])
assert contains_point(image, {u: 1, v: 1}, [u, v])
assert not contains_point(image, {u: 2, v: 0}, [u, v])
```

The matching executable file is `examples/gallery/12_images_and_linear_maps.py`; the documentation
test suite runs every gallery script.
