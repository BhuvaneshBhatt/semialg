# Decision procedures and inequality proving

The decision layer answers exact real-domain questions about formulas.

## Satisfiability and equivalence

```python
import sympy as sp
from semialg import is_satisfiable, is_tautology, implies, equivalent

x, y = sp.symbols("x y", real=True)

is_satisfiable(sp.And(x**2 + y**2 <= 1, x > 0, y > 0), [x, y])
# True

is_tautology(sp.Or(x < 0, x >= 0), [x])
# True

implies(x > 1, x**2 > 1, [x])
# True

equivalent(x**2 <= 1, sp.And(x >= -1, x <= 1), [x])
# True
```

## Inequality proving

The inequality provers reduce sign claims to satisfiability checks:

- `prove_nonnegative(f)` checks unsatisfiability of `f < 0`.
- `prove_positive(f)` checks unsatisfiability of `f <= 0`.
- `prove_nonpositive(f)` checks unsatisfiability of `f > 0`.
- `prove_negative(f)` checks unsatisfiability of `f >= 0`.

```python
from semialg import prove_positive, prove_nonnegative

prove_nonnegative((x - 1)**2, [x])
# True

prove_positive(x**2 + 1, [x])
# True

prove_nonnegative(x*y, [x, y], assumptions=sp.And(x >= 0, y >= 0))
# True
```

## Notes

These functions work over the real domain and are intended for polynomial and semialgebraic formulas. They are often used internally by simplification, region predicates, and optimization routines.
