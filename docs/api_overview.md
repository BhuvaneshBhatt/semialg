# API overview

For an exhaustive, machine-checked list of every root-level export, see the [Public API index](reference/public_api.md).

This page is a map of the public API. It is navigational; signatures, contracts, result semantics, and examples live in the family reference pages.

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

- `cad`, `parametric_cad`
- `extract_structured_cad_cells`
- `extract_vertical_bounds_from_cad_2d`
- `build_cad_adjacency_graph`, `extract_cad_connectivity`
- CAD result, cell, bound, and certificate types

→ [CAD reference](reference/cad.md)

## Optimization and function ranges

- `semialgebraic_minimize`, `semialgebraic_maximize`
- `function_range`
- `polynomial_locus_dimension`
- `critical_values`
- `OptimizationResult` and parametric result types

→ [Optimization and range reference](reference/optimization_and_range.md)

## Regions and geometry

- Extrema sets: `argmin_set`, `argmax_set`, `extrema_set`
- Level families: `level_set`, `sublevel_set`, `superlevel_set`

- Region relations/operations: `is_subset`, `is_equal`, `is_disjoint`, `intersects`, `is_interior_disjoint`, `region_union`, `region_intersection`, `region_difference`, `region_symmetric_difference`, `region_complement`, `region_product`
- Topology: `region_closure`, `region_interior`, `region_boundary`, `is_path_connected`, `path_between`, `euler_characteristic`
- Predicates and structure: `is_empty`, `is_bounded`, `is_compact`, `is_open`, `is_closed`, `is_subset`, `is_equal`, `is_disjoint`, `intersects`, `is_dense_in`, `contains_point`, `is_convex`, `is_connected`, `is_path_connected`, `is_full_dimensional`, `has_empty_interior`, plus the lower-level `region_*` operations
- Maps and metric queries: `semialgebraic_projection`, `region_image`, `region_preimage`, `fiber`, `translate`, `scale`, `linear_image`, `affine_image`, `minkowski_sum`, `squared_distance_range`, `distance_set`, `bounding_box`, `coordinate_range`, `distance_to_region`, `distance_between_regions`, `nearest_point`, `closest_points`, `diameter`, `support_function`, `width`
- Local algebraic geometry: `singular_locus`, `is_singular`, `is_smooth`, `tangent_space`, `tangent_dimension`, `tangent_cone`
- Standard regions: intervals, boxes, balls, spheres, shells, simplices, polytopes, parametric and transformed regions

→ [Regions reference](reference/regions.md)

## Integration, measure, and moments

- `reduce_region_integral`, `integrate_over_region` (including parameter-stratified formula regions)
- `semialgebraic_measure` (including parameter-stratified measure)
- `integrate_over_standard_region`, `integrate_over_parametric_region`
- `region_moment`, `centroid`, `covariance_matrix`, `moment_matrix`, `inertia_tensor`

→ [Integration and moments reference](reference/integration_and_moments.md)

## Algebraic roots and exact finite solving

- `root_of`, `AlgebraicRootFunction`
- `classify_real_roots`
- `semialg.solve.solve_zero_dimensional_system`
- rational-univariate-representation and border-basis APIs

→ [Algebraic reference](reference/algebraic.md)

## Parameters and conditional results

- `conditional_result`
- `ParameterStratifiedResult`
- `verify_parameter_stratification`
- `solvability_conditions` (root); `semialg.parameters.root_count_conditions` (specialist namespace)
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

For a decision table and comparisons between overlapping abstractions, see [Which function should I use?](guides/choosing_an_api.md).

- Need a Boolean answer? Start with the [decision APIs](reference/decision_and_qe.md).
- Need actual points or a solution representation? Use [solving](reference/solving_and_sampling.md).
- Need the decomposition itself? Use [CAD](reference/cad.md).
- Need an extremum? Use [optimization](reference/optimization_and_range.md).
- Need a geometric set operation? Use [regions](reference/regions.md).
- Need measure or an integral? Use [integration](reference/integration_and_moments.md).
- Need behavior as parameters vary? Use [parameter-stratified APIs](reference/parameters.md).

For guarantees shared across all of these families, read [Exactness and certification](concepts/exactness_and_certification.md).


## Errors and failure modes

Package-specific exceptions derive from `semialg.SemialgError` and the canonical classes in `semialg.errors`, including strategy/unsupported-fragment, backend, normalization, algebraic-solving, QE, reconstruction, certification, exact-evaluation, dimension, and resource-limit failures. `semialg.exceptions` exposes the same classes from a dedicated exception module.
