# Root and parameter classification

Root classification is a natural use of CAD/QE because root-existence and root-counting questions can be expressed as real first-order formulas.

## Public API

- `classify_real_roots`
- `solvability_conditions`
- `root_count_conditions`

## Examples

```python
import sympy as sp
from semialg import solvability_conditions, root_count_conditions

x, a, b = sp.symbols("x a b", real=True)

solvability_conditions(sp.Eq(x**2 + a*x + b, 0), [x], [a, b])
# a**2 - 4*b >= 0

root_count_conditions(x**2 + a*x + b, x, [a, b])
# 2 roots if a**2 - 4*b > 0
# 1 root if a**2 - 4*b == 0
# 0 roots if a**2 - 4*b < 0
```

## Current scope

The current implementation has strong support for common univariate parameterized families, especially linear and quadratic cases. Higher-degree parameter families may use CAD-backed classification or conservative fallback behavior depending on the case.
