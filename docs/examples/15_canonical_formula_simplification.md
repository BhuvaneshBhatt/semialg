# Canonical and stable formula simplification

Normalize equivalent polynomial inequalities and zero sets to deterministic forms and verify simplification is idempotent.

## Executable example

```python
import sympy as sp

from semialg import simplify_boole
from semialg.simplify import simplify_semialgebraic_formula

x, y = sp.symbols("x y", real=True)

first = simplify_boole((2 * x - 2 * y > 0) & (x > -1), [x, y], semantic=False)
second = simplify_boole((y - x < 0) & (x > -1), [x, y], semantic=False)
zero_set = simplify_semialgebraic_formula(
    sp.Eq(6 * (x - 1) ** 4 * (x + 2) ** 2, 0), implication_minimize=False
)

print(first)
print(second)
print(zero_set)

assert first == second == ((x > -1) & (x - y > 0))
assert zero_set == sp.Eq(x**2 + x - 2, 0)
assert simplify_boole(first, [x, y], semantic=False) == first
```

The matching executable file is `examples/gallery/15_canonical_formula_simplification.py`; the documentation test suite runs every gallery script.
