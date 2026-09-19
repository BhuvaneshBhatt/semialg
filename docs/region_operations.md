# Region relations and operations

Use these APIs when you want to ask set-theoretic questions about regions or construct a new region from existing geometry. Canonical regions are handled structurally when semialg has an exact closed form; general semialgebraic inputs fall back to the package's exact formula/QE/CAD machinery.

## Relations

```python
from semialg import Ball, IntervalRegion, is_disjoint, is_equal, is_subset

small = Ball((0, 0), 1)
large = Ball((0, 0), 3)

is_subset(small, large)
# True

is_disjoint(IntervalRegion(0, 1), IntervalRegion(2, 3))
# True

is_equal(IntervalRegion(0, 1), IntervalRegion(0, 1))
# True
```

Canonical geometry objects expose the same relations as methods:

```python
small.subset_of(large)
small.disjoint_from(Ball((4, 0), 1))
small.equals_region(Ball((0, 0), 1))
```

The structural layer currently recognizes common exact cases for points, intervals, boxes, and filled balls. If it cannot certify the answer directly, the relation is lowered to the existing exact implication/equivalence/satisfiability engine.

## Boolean set operations

```python
from semialg import region_union, region_intersection, region_difference, region_complement
```

`region_intersection` preserves canonical structure when a direct representation is known. For example, intersecting two closed intervals yields an `IntervalRegion`, and intersecting compatible boxes yields a `BoxRegion`. General intersections use the unified symbolic-region representation.

`region_union`, `region_difference`, and `region_complement` construct exact Boolean semialgebraic regions.

## Cartesian products

`region_product` forms Cartesian products. Closed intervals and boxes retain a canonical `BoxRegion`; products of individual points retain a canonical `Point`.

```python
from semialg import IntervalRegion, region_product

rectangle = region_product(IntervalRegion(0, 1), IntervalRegion(2, 4))
# BoxRegion(bounds=((0, 1), (2, 4)))
```

For regions without a direct canonical product representation, semialg creates a `SemialgebraicRegion` in concatenated coordinates.

## Images and preimages

`region_image` and `region_preimage` represent exact direct and inverse images under affine or polynomial maps. Affine maps preserve canonical structure whenever the resulting geometry has a supported canonical representation. Nonlinear direct images remain exact transformed regions until formula lowering; nonlinear preimages are exact substitutions.

```python
import sympy as sp
from semialg import IntervalRegion, region_preimage

x = sp.symbols("x", real=True)
preimage = region_preimage(IntervalRegion(0, 1), (x**2,), (x,))

preimage.contains((-1,))
# True
preimage.contains((2,))
# False
```

Geometry objects also provide `.image(...)` and `.preimage(...)` methods.

## Minkowski sums

`minkowski_sum(A, B)` returns `{a + b : a in A, b in B}`. Points, intervals, boxes, and filled balls use direct structural formulas. General cases are computed as exact semialgebraic images.

```python
from semialg import Ball, minkowski_sum

minkowski_sum(Ball((0, 0), 1), Ball((2, 0), 3))
# Ball(center=(2, 0), radius=4)
```

The method form `region.minkowski_sum(other)` uses the same dispatch path.

## Closure, interior, boundary, dimension, components

```python
import sympy as sp
from semialg import region_closure, region_interior, region_boundary, region_dimension

x, y = sp.symbols("x y", real=True)

region_closure(sp.And(x > 0, x < 1), [x])
# (x >= 0) & (x <= 1)

region_boundary(sp.And(x > 0, x < 1), [x])
# Eq(x, 0) | Eq(x, 1)

region_dimension(x**2 + y**2 < 1, [x, y])
# 2
```

These topology operations use exact CAD semantics for composite formulas. `region_dimension` is the maximum Euclidean dimension of a selected cell in a complete adapted CAD.

## Symbolic relation conditions

For parameter-dependent statements, `RegionElement`, `RSubset`, `RDisjoint`, and `REqual` represent symbolic assertions. `region_element_conditions()` and `region_relation_conditions()` expose their formulas; use `eliminate=True` when a quantifier-free condition is required.

## Structural convex intersections

Affine lines and rays are clipped directly against convex canonical regions when an exact structural description is available. Boxes use their coordinate halfspaces directly; simplexes, polygons, parallelepipeds, and polytopes use exact polyhedral representations. Filled balls use the quadratic line parameter, so these cases avoid CAD entirely.

Hyperplane sections are also structural for finite convex polyhedra. The intersection vertices are obtained from vertices lying on the hyperplane and exact crossings of polytope edges. Hyperplane sections of filled balls remain intrinsic balls represented through an affine embedding, and hyperplane sections of ellipsoids remain intrinsic ellipsoids.

Polytope relations exploit convexity before using the general decision engine. To prove `P <= Q`, it is enough to certify every generating vertex of `P` against the supporting halfspaces of `Q`. Disjointness first looks for a separating supporting facet. Minkowski sums use the classical identity

```text
conv(V) + conv(W) = conv({v + w : v in V, w in W})
```

and therefore construct the exact sum from pairwise vertex sums without invoking CAD.
