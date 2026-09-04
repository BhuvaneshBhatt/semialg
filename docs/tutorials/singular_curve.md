# Tutorial: singular curve, local dimension, and CAD structure

Use the cusp

```text
y^2 = x^3.
```

```python
import sympy as sp
from semialg import SemialgebraicRegion

x, y = sp.symbols("x y", real=True)
cusp = SemialgebraicRegion(sp.Eq(y**2, x**3), (x, y))
```

The origin is algebraically singular:

```python
singular = cusp.singular_locus()
assert singular.contains((0, 0))
assert not singular.contains((1, 1))
```

The set nevertheless has local dimension one at the cusp:

```python
assert cusp.local_dimension((0, 0)) == 1
assert cusp.local_dimension((1, 1)) == 1
assert cusp.local_dimension((-1, 0)) == -1
```

`local_dimension()` is computed from the exact CAD complex: it takes the maximum dimension of selected cells whose closures contain the query point.

This is intentionally distinct from corner/manifold classification. A rectangular corner can have local dimension two while not being an algebraic hypersurface singularity under `singular_locus()`'s current definition.
