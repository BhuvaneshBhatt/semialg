# Closure and admissible approaches

A deleted point can still be an admissible limit point. This example uses `point_in_closure` on a punctured disk and distinguishes the excluded-but-approachable origin from a genuinely separated point.

```python
import sympy as sp

from semialg import point_in_closure

x, y = sp.symbols("x y", real=True)
punctured_disk = sp.And(x**2 + y**2 < 1, sp.Ne(x**2 + y**2, 0))
assert point_in_closure(punctured_disk, (0, 0), (x, y)) is True
assert point_in_closure(punctured_disk, (2, 0), (x, y)) is False
```
