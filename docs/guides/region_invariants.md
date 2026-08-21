# Region invariants

Standard-region objects validate basic geometric invariants at construction time. Provably invalid geometry is rejected early rather than being allowed to produce negative measure or obscure downstream matrix errors.

Symbolic geometry is allowed when an ordering cannot yet be decided exactly; the checks reject conditions that are **provably** invalid.

| Region | Important constructor invariants |
|---|---|
| `PointRegion` | coordinates define one ambient dimension |
| `IntervalRegion` | lower endpoint must not be provably greater than upper endpoint |
| `BoxRegion` | every interval satisfies lower ≤ upper |
| `SimplexRegion` / `TetrahedronRegion` | vertices have consistent ambient dimension |
| `PolygonRegion` / `PolyhedronRegion` | constituent geometry has consistent dimensions |
| `ParallelogramRegion` / `ParallelepipedRegion` | origin and spanning vectors have compatible dimensions |
| `PrismRegion` | base and extrusion data have compatible ambient dimensions |
| `PyramidRegion` | base and apex have compatible ambient dimensions |
| `BallRegion` / `SphereRegion` | radius is not provably negative |
| `SphericalShellRegion` | $0 \le r_{inner} \le r_{outer}$ when decidable |
| `CylinderRegion` / `ConeRegion` | endpoint dimensions agree and radius is not provably negative |
| `StadiumRegion` / `CapsuleRegion` | endpoint dimensions agree and radius is not provably negative |
| `ParametricRegion` | each declared parameter has exactly one limit; no undeclared limit variables; multiplicity is provably positive |

## Reversed bounds are invalid input

A region such as `IntervalRegion(2, 1)` is not interpreted as an oriented integral. It is invalid geometric input and is rejected. The same rule applies to each coordinate interval in a box and to explicit integration bounds accepted by semialgebraic integration APIs.

## Degenerate versus invalid

Equal endpoints or zero radius may describe a degenerate region and are not the same as reversed bounds or negative radius. Whether a downstream operation supports that lower-dimensional object depends on the operation and requested measure dimension.

## Boolean regions

`RegionUnion`, `RegionIntersection`, and `RegionDifference` preserve set semantics. In particular,

$$
\mu(A\setminus B)=\mu(A)-\mu(A\cap B),
$$

not generally $\mu(A)-\mu(B)$.

## Parametric regions

For a parameterization with parameters $(u_1,\ldots,u_k)$, every parameter must occur exactly once as a limit variable. Unknown, duplicate, or missing limit variables are rejected. Multiplicity must be provably positive because integration divides by it.

## Symbol identity

String variable names in standard/parametric integration are resolved against the actual symbols in the integrand and mapping. See [Symbol handling](symbol_handling.md).
