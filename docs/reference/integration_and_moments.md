# Integration, measure, and moments reference
## Family contract

**Mathematical return.** Integration APIs return exact symbolic integrals over supported semialgebraic regions; measure and moment APIs are derived exact integrals with dimension-aware semantics.

**Exactness and certification.** Bounds/cells used for exact integration are derived from certified semialgebraic decompositions. Numerical quadrature is not silently substituted for an unsupported exact integral.

**Algorithm.** The implementation reduces regions to exact bounds/cells, supports Boolean decomposition and intrinsic regular strata, and can integrate univariate parameter fibers using CAD-controlled algebraic root-function endpoints.

**Complexity and limitations.** Arbitrary singular mixed-dimensional stratification and fully general multidimensional parametric algebraic integration remain incomplete.




## Primary API overview

This table is the substantive coverage target for the primary APIs assigned to this reference page. Each entry states the API's primary role; the family contract and detailed sections below explain shared algorithms, exactness guarantees, and limitations. It is maintained together with `docs/reference/primary_api_manifest.toml`, and documentation tests require every root-level primary API to map here rather than merely appearing in the generated public index.

| API | Kind | Role / return |
|---|---|---|
| `integrate_over_region` | function | Integrate ``integrand`` over a supported semialgebraic region. |
| `reduce_region_integral` | function | Reduce a supported region integral to explicit iterated integrals. |
| `centroid` | function | Return the centroid of a finite-measure semialgebraic region. |
| `covariance_matrix` | function | Return the covariance matrix of the uniform measure on a region. |
| `region_moment` | function | Return a raw moment integral over a semialgebraic region. |
| `region_measure` | function | Return intrinsic or ambient measure for canonical geometry and formula regions. |
| `semialgebraic_measure` | function | Return the exact measure of a supported semialgebraic set. |

## `integrate_over_region`

```text
integrate_over_region(
    integrand, condition, variables, *, bounds=None,
    method="symbolic", precision=50,
    measure_dimension="ambient", return_result=False
)
```

Integrates over a semialgebraic region. Exact symbolic paths use specialized reductions when available and certified CAD-derived cylindrical cells as the general ambient-measure fallback. The integration planner is coordinate-order aware: it can permute integration variables, use Brown-style CAD order scoring, compare resulting bound complexity, and prefer an order with fewer/simpler algebraic boundaries. This avoids introducing radicals merely because the caller supplied an inconvenient coordinate order.

## `reduce_region_integral`

Returns a reduced integral representation without necessarily evaluating it immediately. Useful for inspecting how a region was decomposed. For CAD-derived pieces, diagnostics include the selected `integration_variable_order`, bound-verification status, and whether the order came from explicit cylindrical recognition or CAD order search.

## `region_measure` and geometry methods

`region_measure(region, measure_dimension=None)` is the unified measure front door. With the default `None`, canonical geometry uses intrinsic Hausdorff measure while formula regions use ambient Lebesgue measure. Explicit `"intrinsic"`, `"ambient"`, or integer dimensions override that abstraction-aware default. Canonical `Geometry.measure()` still defaults explicitly to intrinsic measure; `SemialgebraicRegion.measure()` defaults explicitly to ambient measure. Formula-only workflows that need bounds or parameter stratification should use `semialgebraic_measure`.

## `semialgebraic_measure`

Computes ambient measure by default. `measure_dimension="intrinsic"` requests supported intrinsic measure on lower-dimensional regular strata.

## Standard and parametric regions

- `integrate_over_standard_region`
- `integrate_over_parametric_region`
- `reduce_parametric_integral`

Parametric integration preserves the actual ambient SymPy symbol identities and validates parameter limits/multiplicity at region construction.

## Moments

- `region_moment`
- `centroid`
- `covariance_matrix`

These build on the exact region-integration machinery.

## Bounds

Explicit bounds must name declared variables and may not be provably reversed. Malformed, duplicate, or irrelevant bounds are rejected early.

See [Region integration](../region_integration.md), [Moments](../moments.md), and [Region invariants](../guides/region_invariants.md).


### Parameter-dependent formula regions

`integrate_over_region(..., parameters=None, return_stratified=False)` and `semialgebraic_measure(..., parameters=None, return_stratified=False)` accept `parameters=[...]` together with `return_stratified=True` for exact guarded parameter-dependent answers. This API is separate from `integrate_over_parametric_region`, whose `ParametricRegion` argument describes an explicit parametrization.


## Derived moment tensors

`centroid` and `covariance_matrix` are the canonical normalized first- and second-central-moment APIs. `moment_matrix` returns the normalized raw second moment `E[x x.T]`. `inertia_tensor` returns the unit-density tensor

`∫_S (||x||^2 I - x x.T) dx`

about the origin. All four reuse the exact region-integration backend and inherit its supported-region and parameter limitations.
