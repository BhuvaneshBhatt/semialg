# Tutorial: coordinates versus parameters

A symbolic region distinguishes ambient coordinates from parameters:

```python
import sympy as sp
from semialg import SemialgebraicRegion

x, r = sp.symbols("x r", real=True)
interval = SemialgebraicRegion(sp.And(x >= 0, x <= r), (x,))

assert interval.variables == (x,)
assert interval.parameters == (r,)
assert interval.region_variables("all") == (x, r)
```

This distinction matters when performing parameter stratification, projection, or image computations: `r` is not accidentally treated as a coordinate to eliminate.

For a particular parameter value, specialize the formula first and then construct the corresponding region/CAD. For a global parameter analysis, use semialg's parameter-condition and stratification APIs so degree drops and exceptional algebraic boundaries remain explicit.
