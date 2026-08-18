# Region operations and predicates

`semialg` represents semialgebraic regions using SymPy Boolean formulas over real variables.

## Boolean region operations

```python
from semialg import region_union, region_intersection, region_difference, region_complement
```

These functions construct symbolic Boolean formulas for union, intersection, difference, and complement.

## Closure, interior, boundary, dimension, components

```python
import sympy as sp
from semialg import region_closure, region_interior, region_boundary, region_dimension

x, y = sp.symbols("x y", real=True)

region_closure(sp.And(x > 0, x < 1), [x])
# (x >= 0) & (x <= 1)

region_boundary(sp.And(x > 0, x < 1), [x])
# Eq(x, 0) | Eq(x, 1)

region_boundary(x**2 + y**2 < 1, [x, y])
# Eq(x**2 + y**2 - 1, 0) in supported cases

region_dimension(x**2 + y**2 < 1, [x, y])
# 2
```

## Region predicates

```python
from semialg import region_subset, region_equal, region_disjoint, region_bounded, region_closed, region_compact

region_subset(x > 1, x > 0, [x])
# True

region_equal(x**2 <= 1, sp.And(x >= -1, x <= 1), [x])
# True

region_disjoint(x < 0, x > 0, [x])
# True

region_bounded(sp.And(x >= 0, x <= 1), [x])
# True
```

## Current scope

Some region operations are exact for common 1D and 2D semialgebraic forms and use CAD/QE-backed predicates where available. Full topological normalization of arbitrary semialgebraic sets remains future work.
