# Projection as quantifier elimination

Eliminate a coordinate from a parabolic strip and recover the exact projected interval.

## Problem and interpretation

The region is

\[
S=\{(x,y):x^2\le y\le1\}.
\]

Projecting onto the \(x\)-axis asks exactly the quantified question

\[
\exists y\;(x^2\le y\le1).
\]

Such a \(y\) exists iff \(x^2\le1\), hence the projected set is \([-1,1]\).

This is a good example of why a dedicated projection API is preferable to manually assembling a quantified formula: the mathematical intent is explicit, while the backend remains free to use presolve and complete QE/CAD as needed.

## Executable example

```python
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
```

The matching executable file is `examples/gallery/01_projection_and_qe.py`; the documentation
test suite runs every gallery script.
