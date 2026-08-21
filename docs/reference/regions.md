# Regions reference

## Formula-based region operations

- `region_union`
- `region_intersection`
- `region_difference`
- `region_complement`
- `region_closure`
- `region_interior`
- `region_boundary`

Structural and predicate APIs include `region_dimension`, `region_components`, `region_subset`, `region_equal`, `region_disjoint`, `region_bounded`, `region_closed`, and `region_compact`.

## Standard regions

The standard-region hierarchy represents common geometry directly:

- points and intervals;
- boxes;
- simplices and tetrahedra;
- polygons/polyhedra;
- parallelograms/parallelepipeds;
- prisms and pyramids;
- balls, spheres, and spherical shells;
- cylinders and cones;
- stadiums and capsules;
- parametric and transformed regions;
- Boolean combinations of standard regions.

Constructor invariants are part of the public contract. See [Region invariants](../guides/region_invariants.md).

## Boolean standard regions

`RegionUnion`, `RegionIntersection`, `RegionDifference`, and `RegionSymmetricDifference` construct set-theoretic combinations of `StandardRegion` objects.

## Exact endpoint semantics

One-dimensional component merging and region intersections compare exact endpoints symbolically. Close algebraic values are not ordered by fixed-precision decimal conversion.

See [Region operations](../region_operations.md).
