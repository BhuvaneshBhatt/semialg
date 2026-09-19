# Region invariants

Standard-region objects validate basic geometric invariants at construction time. Provably invalid geometry is rejected early rather than being allowed to produce negative measure or obscure downstream matrix errors.

Symbolic geometry is allowed when an ordering cannot yet be decided exactly; most checks reject conditions that are **provably** invalid. `PolygonRegion` is stricter because triangulation requires a certified orientation and simple boundary, so indeterminate symbolic polygon geometry is rejected rather than guessed.

| Region | Important constructor invariants |
|---|---|
| `PointRegion` | coordinates define one ambient dimension |
| `IntervalRegion` | lower endpoint must not be provably greater than upper endpoint |
| `BoxRegion` | every interval satisfies lower ≤ upper |
| `SimplexRegion` / `TetrahedronRegion` | simplices contain at least one vertex; vertices have consistent ambient dimension; tetrahedra use 3-D vertices |
| `PolygonRegion` / `PolyhedronRegion` | polygons are exact simple planar polygons; repeated closing vertices are normalized; self-intersections and zero-area polygons are rejected |
| `ParallelogramRegion` / `ParallelepipedRegion` | origin and spanning vectors have compatible dimensions |
| `PrismRegion` | base and extrusion data have compatible ambient dimensions |
| `PyramidRegion` | base and apex have compatible ambient dimensions |
| `BallRegion` / `SphereRegion` / canonical `Ball` / `Sphere` | radius is not provably negative |
| canonical `Ellipsoid` | center is nonempty; shape matrix is square, symmetric, and positive definite; unresolved symbolic Sylvester conditions are retained |
| `SphericalShellRegion` | $0 \le r_{inner} \le r_{outer}$ when decidable |
| `CylinderRegion` / `ConeRegion` | endpoint dimensions agree and radius is not provably negative |
| `StadiumRegion` / `CapsuleRegion` | endpoint dimensions agree and radius is not provably negative; stadiums are planar while capsules may have arbitrary ambient dimension |
| `ParametricRegion` | each declared parameter has exactly one limit; no undeclared limit variables; multiplicity is provably positive |

## Reversed bounds are invalid input

A region such as `IntervalRegion(2, 1)` is not interpreted as an oriented integral. It is invalid geometric input and is rejected. The same rule applies to each coordinate interval in a box and to explicit integration bounds accepted by semialgebraic integration APIs.

## Degenerate versus invalid

Equal endpoints or zero radius may describe a degenerate region and are not the same as reversed bounds or negative radius. Whether a downstream operation supports that lower-dimensional object depends on the operation and requested measure dimension.

`dimension()` respects certified degeneracy rather than returning the nominal constructor dimension. Examples include a zero-width box coordinate, affinely dependent simplex vertices, a zero-radius ball, and a shell whose inner and outer radii coincide. The empty point set and an open interval with equal endpoints have dimension `-1`.

`PolygonRegion` uses exact ear-clipping triangulation. Concave simple polygons are therefore represented and integrated without the over-counting that a first-vertex fan can introduce. Self-intersecting polygons are invalid because their interior semantics are ambiguous without an explicit winding rule.

## Boolean regions

`RegionUnion`, `RegionIntersection`, and `RegionDifference` preserve set semantics. In particular,

$$
\mu(A\setminus B)=\mu(A)-\mu(A\cap B),
$$

not generally $\mu(A)-\mu(B)$.

Boolean-region constructors also require compatible ambient dimensions and the correct operation arity. Union dimension is structural (`max` of component dimensions); interval intersections are handled exactly. For other Boolean combinations whose dimension depends on geometric incidence, convert to `SemialgebraicRegion` and use `region_dimension()` rather than relying on a structural guess.

## Parametric regions

For a parameterization with parameters $(u_1,\ldots,u_k)$, every parameter must occur exactly once as a limit variable. Unknown, duplicate, or missing limit variables are rejected. Multiplicity must be provably positive because integration divides by it.

## Symbol identity

String variable names in standard/parametric integration are resolved against the actual symbols in the integrand and mapping. See [Symbol handling](symbol_handling.md).
