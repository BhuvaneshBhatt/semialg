# Public API overview

This page groups the public API by task.

## CAD and quantifier elimination

- `cad`
- `cad_text`
- `generic_cad`
- `generic_cad_text`
- `qe_by_complete_cad`
- `reduce_formula`
- `resolve_formula`
- `semialgebraicize`
- `root_of`

## Decision and logic

- `is_satisfiable`
- `is_tautology`
- `implies`
- `equivalent`
- `simplify_boole`
- `BooleanSimplificationResult`
- `simplify_system`

Set `return_result=True` on the four decision helpers to receive structured
objects with validated witnesses or counterexamples. `simplify_boole(..., return_result=True)` returns a `BooleanSimplificationResult`; the simplifier canonicalizes polynomial relational atoms, applies exact univariate interval simplification, and removes semantically redundant clauses when the decision layer can prove implication.

### Structured result types

- `SatisfiabilityResult`
- `TautologyResult`
- `ImplicationResult`
- `EquivalenceResult`
- `BooleanSimplificationResult`
- `SemialgebraicSolution`

## Solving and sampling

- `solve_semialgebraic`
- `find_instance`
- `component_instances`
- `sample_point`
- `sample_points`
- `sign_at`
- `sign_vector`

`sample_points` supports `representative`, `rational`, `grid`, `random`, and
`cad_cells` strategies. `sign_at` and `sign_vector` use exact rational,
algebraic, and RUR-backed sign evaluation by default.

## Parameters and roots

- `classify_real_roots`
- `solvability_conditions`
- `root_count_conditions`

## Region construction and operations

- `region_union`
- `region_intersection`
- `region_difference`
- `region_complement`
- `region_closure`
- `region_interior`
- `region_boundary`
- `region_dimension`
- `region_components`

## Region predicates

- `region_subset`
- `region_equal`
- `region_disjoint`
- `region_bounded`
- `region_closed`
- `region_compact`

## Integration and measure

- `reduce_region_integral`
- `integrate_over_region`
- `semialgebraic_measure`

## Moments and geometry

- `region_moment`
- `region_centroid`
- `region_covariance`

## Optimization and range

- `semialgebraic_minimize`
- `semialgebraic_maximize`
- `function_range`

## Simplification and proving

- `simplify_piecewise`
- `simplify_under_assumptions`
- `AssumptionSimplificationResult`
- `prove_positive`
- `prove_nonnegative`
- `prove_negative`
- `prove_nonpositive`

## Implicit-geometry utilities

- `semialgebraic_level_function`
- `decompose_implicit_formula`
- `extract_symbolic_box_bounds`
- `decompose_cylindrical_formula_to_vertical_bounds_2d`

## Exact algebraic backend utilities

- `subresultant_prs`
- `principal_subresultant_coefficients`
- `compute_border_basis`
- `compute_border_basis_linear`

`compute_border_basis` returns a `BorderBasisResult` with the order ideal,
border relations, exact normal forms, quotient coordinates, multiplication
matrices for variables or arbitrary quotient elements, exact commutator
certificates, and construction diagnostics. The default `algorithm="groebner"`
uses the Groebner-derived exact constructor. `algorithm="linear"` or
`compute_border_basis_linear(...)` uses the exact Macaulay exact Macaulay
linear-algebra constructor, which row-reduces bounded polynomial multiples and
extracts border relations directly from the row space. Use `strict=False` to
receive a failed diagnostic result instead of an exception for invalid custom
order ideals, positive-dimensional inputs, or non-stabilization within the chosen
linear-algebra degree bound.


### Piecewise and system simplification

`simplify_piecewise` simplifies branch expressions under their branch assumptions, removes unreachable branches after earlier conditions cover the domain, and can return a `PiecewiseSimplificationResult` with diagnostics. `simplify_system` records simple equality substitutions in `SimplifiedSystem.substitutions`, supports `eliminate_equalities=True`, and can return either a formula, a tuple of constraints, or the structured result via `output=`.

### Structured sign proof results

The sign-proving helpers accept `return_result=True` and return
`SignProofResult`:

```python
result = prove_positive(x**2, [x], return_result=True)
result.proven          # False
result.counterexample  # {x: 0}
```

Default calls continue to return plain booleans. The implementation tries cheap
certificates before the general implication engine: constants, obvious sums of
squares/even-power products, negative sums of squares, and positive
constant-plus-squares forms.
