# Symbolic optimization

`semialg` includes initial exact global optimization routines over supported semialgebraic domains.

## Public API

- `semialgebraic_minimize`
- `semialgebraic_maximize`
- `function_range`

## Examples

```python
import sympy as sp
from semialg import semialgebraic_minimize, semialgebraic_maximize

x, y = sp.symbols("x y", real=True)

semialgebraic_minimize(x**2, [x >= 2], [x])
# value 4 at x = 2

semialgebraic_minimize(x, [x > 0], [x])
# infimum 0, not attained

semialgebraic_maximize(x*y, [x >= 0, y >= 0, x + y <= 1], [x, y])
# value 1/4 at (1/2, 1/2)
```

## Relationship to quantifier elimination

Global optimization can be expressed by projecting objective-level sets. For a polynomial objective `f` and constraints `C`, the set of achievable values below `t` is described by:

```text
exists x. C(x) and f(x) <= t
```

The current optimizer includes exact univariate support and selected low-dimensional polynomial cases. A more complete optimizer would use full objective-value projection, KKT/critical-point methods, and possibly SOS/SDP certificate backends.
