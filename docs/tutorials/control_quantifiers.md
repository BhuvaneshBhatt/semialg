# Adaptive versus fixed control

Use this workflow to decide whether a response may depend on a disturbance. A prefix is written from the outermost quantifier to the innermost. The actuator can copy each disturbance in [-1,1], but a single fixed setting cannot copy every disturbance. Quantifiers range over the reals; explicit guards express the operating domain.

<!-- semialg-exec -->
```python
"""Adaptive and fixed control choices have different quantified specifications."""

import sympy as sp

from semialg import quantifier_eliminate

u, d = sp.symbols("u d", real=True)
disturbance = sp.And(d >= -1, d <= 1)
response = sp.And(u >= -1, u <= 1, sp.Eq(u, d))
specification = sp.Or(sp.Not(disturbance), response)
adaptive = quantifier_eliminate(specification, (("forall", d), ("exists", u)))
fixed = quantifier_eliminate(specification, (("exists", u), ("forall", d)))
assert adaptive is sp.true
assert fixed is sp.false
```

Run the [standalone example](https://github.com/BhuvaneshBhatt/semialg/blob/main/examples/control_quantifiers.py) against an installed wheel. See the [replay matrix](../quality/certificate_matrix.md) and [failure contracts](../guides/errors_and_failure_modes.md) for the scope of each guarantee.
