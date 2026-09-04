# Tutorial: exact optimization over a region

A region formula can be supplied directly to the optimization API:

```python
import sympy as sp
from semialg import SemialgebraicRegion, semialgebraic_minimize

x, y = sp.symbols("x y", real=True)
region = SemialgebraicRegion(
    sp.And(x >= 0, y >= 0, x + y <= 3),
    (x, y),
)

result = semialgebraic_minimize(x + 2*y, region.formula, region.variables, return_result=True)
```

For bounded closed affine polytopes with an affine objective, semialg uses the exact vertex specialization before the general KKT/CAD pipeline. More complicated problems retain exact comparison and can use CAD for global certification.

When a public result carries replay data, use:

```python
from semialg import replay_certificate

replay = replay_certificate(result)
```

See the certificate matrix for the distinction between `verified=True`, `False`, and `None` (insufficient retained replay payload).
