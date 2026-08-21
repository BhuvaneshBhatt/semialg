# API overview

This page is a map of the public API. It is intentionally navigational; signatures, contracts, result semantics, and examples live in the family reference pages.

## Decision and quantifier elimination

Use these when the primary question is whether a real formula is true, feasible, implied, or equivalent:

- `is_satisfiable`
- `is_tautology`
- `implies`
- `equivalent`
- `qe_by_complete_cad`

→ [Decision and QE reference](reference/decision_and_qe.md)

## Solving, witnesses, and sampling

- `solve_semialgebraic`
- `sample_point`, `sample_points`
- `sign_at`, `sign_vector`
- `find_instance_formula`, `component_instances`

→ [Solving and sampling reference](reference/solving_and_sampling.md)

## CAD and structured geometry

- `cad`, `generic_cad`
- `extract_structured_cad_cells`
- `extract_vertical_bounds_from_cad_2d`
- `build_cad_adjacency_graph`, `extract_cad_connectivity`
- CAD result, cell, bound, and certificate types

→ [CAD reference](reference/cad.md)

## Optimization and function ranges

- `semialgebraic_minimize`, `semialgebraic_maximize`
- `function_range`
- `polynomial_locus_dimension`
- `OptimizationResult` and parametric result types

→ [Optimization and range reference](reference/optimization_and_range.md)

## Regions and geometry

- Boolean operations: `region_union`, `region_intersection`, `region_difference`, `region_complement`
- Topology: `region_closure`, `region_interior`, `region_boundary`
- Predicates and structure: `region_dimension`, `region_components`, `region_subset`, `region_equal`, `region_disjoint`, `region_bounded`, `region_closed`, `region_compact`
- Standard regions: intervals, boxes, balls, spheres, shells, simplices, polytopes, parametric and transformed regions

→ [Regions reference](reference/regions.md)

## Integration, measure, and moments

- `reduce_region_integral`, `integrate_over_region`
- `semialgebraic_measure`
- `integrate_over_standard_region`, `integrate_over_parametric_region`
- `region_moment`, `region_centroid`, `region_covariance`

→ [Integration and moments reference](reference/integration_and_moments.md)

## Algebraic roots and exact finite solving

- `root_of`, `AlgebraicRootFunction`
- `classify_real_roots`
- `solve_zero_dimensional_system`
- rational-univariate-representation and border-basis APIs

→ [Algebraic reference](reference/algebraic.md)

## Parameters and conditional results

- `conditional_result`
- `ParameterStratifiedResult`
- `verify_parameter_stratification`
- `solvability_conditions`, `root_count_conditions`
- parametric optimization/range result types

→ [Parameters reference](reference/parameters.md)


## Applied workflows

- robust parameter/tolerance analysis
- certified symbolic-math validation
- exact numerical-optimization benchmarks
- polynomial control stability regions
- polynomial safety/invariant verification
- polynomial response-surface analysis
- polynomial model comparison
- parameter solvability and root-count regimes
- polynomial/geometric probability

→ [Applications reference](reference/applications.md) — robust design, validation, control/safety verification, response surfaces, model comparison, parameter regimes, probability, Lyapunov/barrier certificates, sensitivity, and constraint diagnostics

## Choosing an API

- Need a Boolean answer? Start with the [decision APIs](reference/decision_and_qe.md).
- Need actual points or a solution representation? Use [solving](reference/solving_and_sampling.md).
- Need the decomposition itself? Use [CAD](reference/cad.md).
- Need an extremum? Use [optimization](reference/optimization_and_range.md).
- Need a geometric set operation? Use [regions](reference/regions.md).
- Need measure or an integral? Use [integration](reference/integration_and_moments.md).
- Need behavior as parameters vary? Use [parameter-stratified APIs](reference/parameters.md).

For guarantees shared across all of these families, read [Exactness and certification](concepts/exactness_and_certification.md).
