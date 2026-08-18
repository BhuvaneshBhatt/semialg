# Solving and sampling semialgebraic systems

`semialg` provides wrappers for deciding feasibility and producing representative sample points.

## `solve_semialgebraic`

```python
import sympy as sp
from semialg import solve_semialgebraic

x, y = sp.symbols("x y", real=True)

sol = solve_semialgebraic(
    [x**2 + y**2 <= 1, x > 0, y > 0],
    [x, y],
)

sol.satisfiable
# True

sol.sample
# one satisfying sample point, when available
```

A `SemialgebraicSolution` contains:

- `formula`
- `variables`
- `satisfiable`
- `sample`
- `samples`
- `method`
- `diagnostics`

## Sampling and sign evaluation

```python
from semialg import sample_point, sample_points, sign_at, sign_vector

pt = sample_point(sp.And(x > 0, x < 1), [x])
sign_at(x - sp.Rational(1, 2), pt)

sign_vector([x, x - 1], {x: sp.Rational(1, 2)})
# (+, -) in the package's sign representation
```

## Current scope

The solver emphasizes exact feasibility and representative samples. It does not yet guarantee a human-minimal symbolic description of every solution set.
