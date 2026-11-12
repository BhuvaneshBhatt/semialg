# Getting started

Install semialg with algroots 0.4.1 or newer. SymPy 1.14 or newer is required by the exact algebra backend. Start with the executable README example, then follow [worked applications](worked_applications.md).

<!-- examples: executable -->
```python
import sympy as sp
from semialg.qe import quantifier_eliminate
from semialg.decision import equivalent

x, a = sp.symbols("x a", real=True)
condition = quantifier_eliminate(sp.Eq(x**2, a), [("exists", x)], variables=(a, x))
assert equivalent(condition, a >= 0, (a,))
```

The result is an exact parameter condition, it is neither a numerical sample nor a statement only about generic coefficients. Continue with [API selection](guides/choosing_an_api.md), [result meanings](concepts/exactness_and_certification.md) and [limits](limitations.md).
