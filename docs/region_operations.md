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
from semialg.reasoning import (
    region_bounded, region_closed, region_compact,
    region_disjoint, region_equal, region_subset,
)

region_subset(x > 1, x > 0, [x])
# True

region_equal(x**2 <= 1, sp.And(x >= -1, x <= 1), [x])
# True

region_disjoint(x < 0, x > 0, [x])
# True

region_bounded(sp.And(x >= 0, x <= 1), [x])
# True
```

## CAD-semantic topology

`region_closure`, `region_interior`, and `region_boundary` use an adapted complete CAD for composite Boolean semialgebraic formulas. Closure is computed from CAD-cell incidence, interior as the complement of the closure of the complement, and boundary as the intersection of the two closures. This correctly removes internal decomposition seams; for example, the shared point in `[0, 1] ∪ [1, 2]` is an interior point rather than a boundary point.

Atomic polynomial relations use the same CAD semantics as composite Boolean formulas. A compact atom-wise formula is retained only when its truth set is verified against the resulting CAD cells. This matters for cases such as `x**2*y > 0`, whose closure is not obtained by merely replacing `>` with `>=`. The atom-wise transformation remains available only with `strategy="syntactic"`.

### Exact dimension

`region_dimension` is an exact operation: it returns the maximum Euclidean dimension of a selected cell in a complete adapted CAD. The former equation-count/interior heuristic is no longer used by the public API.

## Symbolic membership and relation conditions

For symbolic region predicates, condition conversion is deliberately non-CAD by default.

```python
from semialg.symbolic_regions import region_element_conditions, region_relation_conditions

region_element_conditions(expr)                       # structural lowering only
region_element_conditions(expr, eliminate=True)       # request exact QE/CAD
region_element_conditions(expr, real_parameters=True) # add Contains(p, Reals) guards

region_relation_conditions(expr)
region_relation_conditions(expr, eliminate=True)
```

`real_parameters=True` makes real-domain assumptions for algebraic-level region parameters explicit in the returned Boolean formula.

## Singular, nonsmooth, and active-boundary structure

`region_singular_locus` uses reduced-real component-relative Jacobian regularity. A smooth component contributes no intrinsic singularities, while intersections of distinct certified components contribute union-induced singularities even when the branches are individually smooth or tangent. Inequality boundaries are analyzed only on exact boundary CAD cells where the corresponding residual is actually active, and relative to each certified equality component rather than a global maximum dimension. The final locus is restricted to the actual topological boundary, so redundant or inactive inequalities cannot manufacture singularities. If a required equidimensional decomposition cannot be certified complete, the formula-only API raises rather than falling back to a global-rank approximation; `region_singular_locus_result` exposes the incomplete status and any certified lower-bound singular subset explicitly. Internally, exact incidence strata distinguish pairwise from higher-order component intersections and exclude components not in the recorded incidence set. `region_nonsmooth_locus` additionally detects transverse intersections on realized multi-active inequality strata, such as polygon corners and polyhedral ridges. `region_active_boundary_strata` reuses the rich boundary CAD metadata and returns pairwise-disjoint exact strata classified by the realized active residual set.

```python
from semialg import region_active_boundary_strata, region_nonsmooth_locus

square = sp.And(x >= 0, x <= 1, y >= 0, y <= 1)
strata = region_active_boundary_strata(square, [x, y])
assert sum(s.active_count == 1 for s in strata) == 4
assert sum(s.active_count == 2 for s in strata) == 4
```

The active-boundary result is a defining-constraint stratification, not a Whitney stratification; redundant inequalities can refine it.
