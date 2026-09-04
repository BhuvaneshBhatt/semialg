# Symbolic region conditions and active-boundary strata

Region predicates are structural by default: they do not unexpectedly launch CAD/QE. Exact elimination is opt-in, and parameter reality can be made explicit when a symbolic workflow needs those assumptions in the returned formula.

```python
import sympy as sp

from semialg import ParametricRegion, region_active_boundary_strata
from semialg.symbolic_regions import (
    RegionElement,
    SemialgebraicRegion,
    as_semialgebraic_region,
    region_element_conditions,
)

x, y, t = sp.symbols("x y t", real=True)
a = sp.Symbol("a")

mapped = as_semialgebraic_region(ParametricRegion((t,), ((t, 0, 1),), (t, t**2)), (x, y))
membership = RegionElement((x, y), mapped)
symbolic = region_element_conditions(membership)
eliminated = region_element_conditions(membership, eliminate=True)
assert symbolic == membership
assert eliminated != membership

parameter_region = SemialgebraicRegion(sp.And(x >= 0, x <= a), (x,))
parameter_atom = RegionElement((x,), parameter_region)
explicit_domains = region_element_conditions(parameter_atom, real_parameters=True)
assert explicit_domains.has(sp.Contains(a, sp.S.Reals, evaluate=False))

square = sp.And(x >= 0, x <= 1, y >= 0, y <= 1)
strata = region_active_boundary_strata(square, (x, y))
assert sum(stratum.active_count == 1 for stratum in strata) == 4
assert sum(stratum.active_count == 2 for stratum in strata) == 4

print("symbolic membership:", symbolic)
print("with explicit real parameters:", explicit_domains)
print("eliminated membership:", eliminated)
print("active boundary strata:", len(strata))
```

The mapped-region example deliberately contains an existential parameter after coercion, so the default `region_element_conditions(...)` leaves the symbolic membership atom intact. `eliminate=True` requests exact QE/CAD instead. The separate algebraic-parameter example shows `real_parameters=True`, which adds an explicit `Contains(a, Reals)` guard.

`region_active_boundary_strata` classifies boundary pieces by the recognized inequality residuals that are active. For the square this gives four one-active edge strata and four two-active corner strata. This is a defining-constraint stratification, not a claim of Whitney regularity; redundant inequalities may refine the returned strata.

The matching executable file is `examples/gallery/17_symbolic_region_conditions.py`.
