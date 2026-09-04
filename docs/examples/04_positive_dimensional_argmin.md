# An optimizer locus with positive dimension

Return the entire minimizer set rather than a single witness.

## Problem and interpretation

Minimizing \(x^2\) over the rectangle does not select one point. Every point on

\[
\{(0,y):0\le y\le1\}
\]

is a global minimizer.

This is exactly the situation in which `argmin_set` is preferable to an
optimization result containing one or several representative optimizer points.
The returned semialgebraic formula describes the complete optimizer locus.

## Executable example

```python
import sympy as sp

from semialg import argmin_set, is_equal

x, y = sp.symbols("x y", real=True)
box = (x >= -1) & (x <= 1) & (y >= 0) & (y <= 1)

minimizers = argmin_set(x**2, box, [x, y])
expected = sp.Eq(x, 0) & (y >= 0) & (y <= 1)

print(minimizers)
assert is_equal(minimizers, expected, [x, y])
```

The matching executable file is `examples/gallery/04_positive_dimensional_argmin.py`; the documentation
test suite runs every gallery script.
