# Integration, measure, and moments reference

## `integrate_over_region`

```text
integrate_over_region(
    integrand, condition, variables, *, bounds=None,
    method="symbolic", precision=50,
    measure_dimension="ambient", return_result=False
)
```

Integrates over a semialgebraic region. Exact symbolic paths use specialized reductions when available and CAD-derived structure where appropriate.

## `reduce_region_integral`

Returns a reduced integral representation without necessarily evaluating it immediately. Useful for inspecting how a region was decomposed.

## `semialgebraic_measure`

Computes ambient measure by default. `measure_dimension="intrinsic"` requests supported intrinsic measure on lower-dimensional regular strata.

## Standard and parametric regions

- `integrate_over_standard_region`
- `integrate_over_parametric_region`
- `reduce_parametric_integral`

Parametric integration preserves the actual ambient SymPy symbol identities and validates parameter limits/multiplicity at region construction.

## Moments

- `region_moment`
- `region_centroid`
- `region_covariance`

These build on the exact region-integration machinery.

## Bounds

Explicit bounds must name declared variables and may not be provably reversed. Malformed, duplicate, or irrelevant bounds are rejected early.

See [Region integration](../region_integration.md), [Moments](../moments.md), and [Region invariants](../guides/region_invariants.md).
