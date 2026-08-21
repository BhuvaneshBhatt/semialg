# Algebraic roots and exact finite solving reference

## `root_of` and `AlgebraicRootFunction`

`root_of(polynomial, variable, index)` denotes an ordered real algebraic root with binder-aware substitution semantics. `AlgebraicRootFunction` carries the root identity as base parameters vary and supports certified specialization, comparison, derivatives, and regularity checks in supported cases.

When a fully specialized low-degree root has a uniquely certified exact radical presentation, semialg may display that radical while retaining the ordered-root identity as the certification basis.

## `classify_real_roots`

```text
classify_real_roots(polynomial, variable, *, parameters=None)
```

Classifies real-root behavior, including parameter-dependent cases. String variables/parameters resolve to the actual expression symbols.

## `solve_zero_dimensional_system`

```text
solve_zero_dimensional_system(
    equations, inequalities=None, vars=None, *, variables=None,
    backend="rur", real=True, parameter=None,
    max_separating_attempts=64
)
```

Solves supported finite polynomial systems exactly using rational-univariate representation. Optional inequalities are exact filters on algebraic candidate points.

## Exact comparisons

Certified algebraic ordering uses isolating intervals, minimal/defining polynomials, exact sign determination, and appropriate algebraic representations. Fixed-precision sorting is not a proof mechanism.

## Advanced algebraic APIs

The package also exposes RUR, border-basis, subresultant, and related algebraic utilities for advanced users. These are lower-level than the primary decision/solve interfaces and may require stronger preconditions.

See [Root classification](../root_classification.md) and [Exactness and certification](../concepts/exactness_and_certification.md).
