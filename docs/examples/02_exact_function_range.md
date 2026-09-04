# Exact function range on a disk

Compute every attainable value of a linear polynomial on the unit disk.

## Problem and interpretation

For \(f(x,y)=x+y\), Cauchy–Schwarz gives the familiar bound
\(|f|\le\sqrt2\), and `semialg` reconstructs that range exactly as a
semialgebraic condition on its generated value symbol.

The important distinction is that `function_range` computes the **whole image**
of the region under \(f\), not merely its two extreme values. For a disconnected
domain or nonlinear map, the range need not be a single interval.

## Executable example

```python
import sympy as sp

from semialg import function_range

x, y = sp.symbols("x y", real=True)
disk = x**2 + y**2 <= 1

value_set = function_range(x + y, disk, [x, y])

print(value_set)
t = next(iter(value_set.free_symbols - {x, y}))
assert sp.simplify(value_set.subs(t, sp.sqrt(2))) is sp.true
assert sp.simplify(value_set.subs(t, 2)) is sp.false
```

The matching executable file is `examples/gallery/02_exact_function_range.py`; the documentation
test suite runs every gallery script.
