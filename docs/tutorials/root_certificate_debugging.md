# Debug a root certificate

Use this workflow when a root count seems inconsistent with an isolating interval. Counts are of distinct roots, so the repeated root at 1 counts once. Open endpoints remove both endpoint roots. Replay checks the retained arithmetic; expected-input arguments also prevent a valid proof for another problem from being reused.

<!-- semialg-exec -->
```python
"""Root counts, endpoint conventions and expected-statement binding."""

import pickle
from dataclasses import replace

import sympy as sp

from semialg.algebraic import certify_polynomial_root_interval, verify_root_interval_certificate
from semialg.cache_control import clear_caches

x = sp.Symbol("x")
polynomial = x * (x - 1) ** 2
closed = certify_polynomial_root_interval(polynomial, 0, 1, var=x)
opened = certify_polynomial_root_interval(
    polynomial, 0, 1, var=x, include_left=False, include_right=False
)
assert closed.root_count == 2
assert opened.root_count == 0
clear_caches(collect=False)
restored = pickle.loads(pickle.dumps(closed))
assert verify_root_interval_certificate(restored, polynomial=polynomial, var=x, left=0, right=1)
assert not verify_root_interval_certificate(replace(restored, root_count=1))
assert not verify_root_interval_certificate(restored, polynomial=x * (x - 2), var=x)
```

Run the [standalone example](https://github.com/BhuvaneshBhatt/semialg/blob/main/examples/root_certificate_debugging.py) against an installed wheel. See the [replay matrix](../quality/certificate_matrix.md) and [failure contracts](../guides/errors_and_failure_modes.md) for the scope of each guarantee.
