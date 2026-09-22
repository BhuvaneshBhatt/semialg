# Root-function documentation adequacy

This audit checks documentation *quality dimensions* for every root-level function, not merely whether the name appears in the reference index. It is generated from the live API, docstrings, the primary-reference manifest, and executable examples in fenced code blocks.

The rubric deliberately distinguishes reference coverage from adequacy. Medium-risk functions require a clear purpose, inspectable signature, return semantics, and family-level exactness/limitations. High-risk functions additionally require parameter semantics and an executable-style documentation example. Critical functions additionally require algorithm/backend documentation.

**Current result:** 218/218 functions satisfy their risk-adjusted adequacy requirements; 0 need documentation deepening.

## Gap summary

| Dimension | Functions missing required evidence |
|---|---:|
| `purpose` | 0 |
| `signature` | 0 |
| `parameters` | 0 |
| `return_semantics` | 0 |
| `exactness` | 0 |
| `algorithm` | 0 |
| `limitations` | 0 |
| `example` | 0 |

## Function-by-function audit

| Function | Owner | Risk | Status | Missing dimensions | Reference |
|---|---|---|---|---|---|
| `quantifier_eliminate` | geometry | high | adequate | — | [reference/root_api_usage.md](../reference/root_api_usage.md#quantifier_eliminate) |
| `project_region` | decision | high | adequate | — | [reference/decision_and_qe.md](../reference/decision_and_qe.md#project_regionregion-eliminate-variablesnone-strategyauto-return_resultfalse) |
| `convert_region` | geometry | high | adequate | — | [reference/root_api_usage.md](../reference/root_api_usage.md#convert_region) |
| `polygonal_region_from_paths` | geometry | high | adequate | — | [reference/root_api_usage.md](../reference/root_api_usage.md#polygonal_region_from_paths) |
| `triangulate_polytope` | geometry | high | adequate | — | [reference/root_api_usage.md](../reference/root_api_usage.md#triangulate_polytope) |
| `decompose_polytope` | geometry | high | adequate | — | [reference/root_api_usage.md](../reference/root_api_usage.md#decompose_polytope) |
| `tetrahedralize_cell` | geometry | high | adequate | — | [reference/root_api_usage.md](../reference/root_api_usage.md#tetrahedralize_cell) |
| `tetrahedralize_cells` | geometry | high | adequate | — | [reference/root_api_usage.md](../reference/root_api_usage.md#tetrahedralize_cells) |
| `canonicalize_polygon` | geometry | high | adequate | — | [reference/root_api_usage.md](../reference/root_api_usage.md#canonicalize_polygon) |
| `canonicalize_polyhedron` | geometry | high | adequate | — | [reference/root_api_usage.md](../reference/root_api_usage.md#canonicalize_polyhedron) |
| `canonicalize_region` | geometry | high | adequate | — | [reference/root_api_usage.md](../reference/root_api_usage.md#canonicalize_region) |
| `verify_nonnegative_combination_certificate` | geometry | high | adequate | — | [reference/root_api_usage.md](../reference/root_api_usage.md#verify_nonnegative_combination_certificate) |
| `nonnegative_combination_certificate` | geometry | high | adequate | — | [reference/root_api_usage.md](../reference/root_api_usage.md#nonnegative_combination_certificate) |
| `implied_polynomial_inequality` | geometry | high | adequate | — | [reference/root_api_usage.md](../reference/root_api_usage.md#implied_polynomial_inequality) |
| `redundant_polynomial_inequalities` | geometry | high | adequate | — | [reference/root_api_usage.md](../reference/root_api_usage.md#redundant_polynomial_inequalities) |
| `component_constraint_descriptions` | geometry | high | adequate | — | [reference/root_api_usage.md](../reference/root_api_usage.md#component_constraint_descriptions) |
| `polynomial_constraints` | algebraic | high | adequate | — | [reference/root_api_usage.md](../reference/root_api_usage.md#polynomial_constraints) |
| `active_constraints` | algebraic | high | adequate | — | [reference/root_api_usage.md](../reference/root_api_usage.md#active_constraints) |
| `relative_interior` | algebraic | high | adequate | — | [reference/root_api_usage.md](../reference/root_api_usage.md#relative_interior) |
| `relative_boundary` | algebraic | high | adequate | — | [reference/root_api_usage.md](../reference/root_api_usage.md#relative_boundary) |
| `semialgebraic_tangent_cone` | algebraic | high | adequate | — | [reference/root_api_usage.md](../reference/root_api_usage.md#semialgebraic_tangent_cone) |
| `parameterization_geometry` | geometry | high | adequate | — | [reference/root_api_usage.md](../reference/root_api_usage.md#parameterization_geometry) |
| `parameterization_critical_locus` | geometry | high | adequate | — | [reference/root_api_usage.md](../reference/root_api_usage.md#parameterization_critical_locus) |
| `parameterization_critical_values` | geometry | high | adequate | — | [reference/root_api_usage.md](../reference/root_api_usage.md#parameterization_critical_values) |
| `implicitize_polynomial_map` | algebraic | high | adequate | — | [reference/root_api_usage.md](../reference/root_api_usage.md#implicitize_polynomial_map) |
| `minimal_prime_intersections` | algebraic | high | adequate | — | [reference/root_api_usage.md](../reference/root_api_usage.md#minimal_prime_intersections) |
| `local_dimension_strata` | algebraic | high | adequate | — | [reference/root_api_usage.md](../reference/root_api_usage.md#local_dimension_strata) |
| `local_branch_geometry` | algebraic | high | adequate | — | [reference/root_api_usage.md](../reference/root_api_usage.md#local_branch_geometry) |
| `stratified_singular_geometry` | algebraic | high | adequate | — | [reference/root_api_usage.md](../reference/root_api_usage.md#stratified_singular_geometry) |
| `certified_radicalization` | algebraic | high | adequate | — | [reference/root_api_usage.md](../reference/root_api_usage.md#certified_radicalization) |
| `irreducible_components` | algebraic | high | adequate | — | [reference/root_api_usage.md](../reference/root_api_usage.md#irreducible_components) |
| `reduced_component_singular_loci` | algebraic | high | adequate | — | [reference/root_api_usage.md](../reference/root_api_usage.md#reduced_component_singular_loci) |
| `zariski_closure` | algebraic | high | adequate | — | [reference/root_api_usage.md](../reference/root_api_usage.md#zariski_closure) |
| `is_singular` | algebraic | high | adequate | — | [reference/root_api_usage.md](../reference/root_api_usage.md#is_singular) |
| `is_smooth` | algebraic | high | adequate | — | [reference/root_api_usage.md](../reference/root_api_usage.md#is_smooth) |
| `singular_locus` | algebraic | high | adequate | — | [reference/root_api_usage.md](../reference/root_api_usage.md#singular_locus) |
| `tangent_cone` | geometry | high | adequate | — | [reference/regions.md](../reference/regions.md#primary-api-overview) |
| `tangent_dimension` | algebraic | high | adequate | — | [reference/root_api_usage.md](../reference/root_api_usage.md#tangent_dimension) |
| `tangent_space` | geometry | high | adequate | — | [reference/regions.md](../reference/regions.md#primary-api-overview) |
| `as_cad_region` | geometry | high | adequate | — | [reference/root_api_usage.md](../reference/root_api_usage.md#as_cad_region) |
| `replay_certificate` | geometry | critical | adequate | — | [reference/regions.md](../reference/regions.md#primary-api-overview) |
| `analyze_affine_map` | geometry | high | adequate | — | [reference/root_api_usage.md](../reference/root_api_usage.md#analyze_affine_map) |
| `convexity_certificate` | function-analysis | medium | adequate | — | [reference/root_api_usage.md](../reference/root_api_usage.md#convexity_certificate) |
| `function_convexity` | function-analysis | medium | adequate | — | [reference/convexity_backend.md](../reference/convexity_backend.md#function-convexity) |
| `function_convex_partition` | function-analysis | medium | adequate | — | [reference/convexity_backend.md](../reference/convexity_backend.md#function-convex-partition) |
| `function_monotonicity` | function-analysis | medium | adequate | — | [reference/convexity_backend.md](../reference/convexity_backend.md#function-monotonicity) |
| `function_monotonic_partition` | function-analysis | medium | adequate | — | [reference/convexity_backend.md](../reference/convexity_backend.md#function-monotonic-partition) |
| `function_sign_partition` | function-analysis | medium | adequate | — | [reference/convexity_backend.md](../reference/convexity_backend.md#function-sign-partition) |
| `function_smoothness` | function-analysis | medium | adequate | — | [reference/root_api_usage.md](../reference/root_api_usage.md#function_smoothness) |
| `function_mapping_properties` | function-analysis | medium | adequate | — | [reference/root_api_usage.md](../reference/root_api_usage.md#function_mapping_properties) |
| `is_injective` | function-analysis | medium | adequate | — | [reference/convexity_backend.md](../reference/convexity_backend.md#function-mapping-properties) |
| `is_surjective` | function-analysis | medium | adequate | — | [reference/convexity_backend.md](../reference/convexity_backend.md#function-mapping-properties) |
| `is_bijective` | function-analysis | medium | adequate | — | [reference/convexity_backend.md](../reference/convexity_backend.md#function-mapping-properties) |
| `is_function_continuous` | function-analysis | medium | adequate | — | [reference/convexity_backend.md](../reference/convexity_backend.md#function-smoothness) |
| `is_function_smooth` | function-analysis | medium | adequate | — | [reference/convexity_backend.md](../reference/convexity_backend.md#function-smoothness) |
| `matrix_definiteness` | function-analysis | medium | adequate | — | [reference/convexity_backend.md](../reference/convexity_backend.md#family-contract) |
| `matrix_pd_on` | function-analysis | medium | adequate | — | [reference/convexity_backend.md](../reference/convexity_backend.md#family-contract) |
| `matrix_psd_on` | function-analysis | medium | adequate | — | [reference/convexity_backend.md](../reference/convexity_backend.md#family-contract) |
| `matrix_rank_on` | function-analysis | medium | adequate | — | [reference/convexity_backend.md](../reference/convexity_backend.md#family-contract) |
| `matrix_rank_stratification` | function-analysis | medium | adequate | — | [reference/convexity_backend.md](../reference/convexity_backend.md#family-contract) |
| `is_convex` | function-analysis | medium | adequate | — | [reference/root_api_usage.md](../reference/root_api_usage.md#is_convex) |
| `equivalent` | decision | critical | adequate | — | [reference/root_api_usage.md](../reference/root_api_usage.md#equivalent) |
| `implies` | decision | critical | adequate | — | [reference/decision_and_qe.md](../reference/decision_and_qe.md#primary-api-overview) |
| `is_satisfiable` | decision | critical | adequate | — | [reference/decision_and_qe.md](../reference/decision_and_qe.md#primary-api-overview) |
| `real_algebraic_feasibility` | algebraic | critical | adequate | — | [reference/root_api_usage.md](../reference/root_api_usage.md#real_algebraic_feasibility) |
| `solve_real_algebraic_set` | algebraic | critical | adequate | — | [reference/root_api_usage.md](../reference/root_api_usage.md#solve_real_algebraic_set) |
| `find_negative_point` | geometry | critical | adequate | — | [reference/root_api_usage.md](../reference/root_api_usage.md#find_negative_point) |
| `polynomial_nonnegative` | geometry | critical | adequate | — | [reference/root_api_usage.md](../reference/root_api_usage.md#polynomial_nonnegative) |
| `zeng_negative_point` | geometry | high | adequate | — | [reference/root_api_usage.md](../reference/root_api_usage.md#zeng_negative_point) |
| `find_negative_witness_fast` | geometry | critical | adequate | — | [reference/root_api_usage.md](../reference/root_api_usage.md#find_negative_witness_fast) |
| `is_tautology` | decision | critical | adequate | — | [reference/decision_and_qe.md](../reference/decision_and_qe.md#primary-api-overview) |
| `solve_semialgebraic` | decision | critical | adequate | — | [reference/root_api_usage.md](../reference/root_api_usage.md#solve_semialgebraic) |
| `cad` | geometry | critical | adequate | — | [reference/root_api_usage.md](../reference/root_api_usage.md#cad) |
| `parametric_cad` | geometry | high | adequate | — | [reference/root_api_usage.md](../reference/root_api_usage.md#parametric_cad) |
| `argmax_set` | geometry | high | adequate | — | [reference/root_api_usage.md](../reference/root_api_usage.md#argmax_set) |
| `argmin_set` | geometry | high | adequate | — | [reference/regions.md](../reference/regions.md#primary-api-overview) |
| `centroid` | integration | high | adequate | — | [reference/root_api_usage.md](../reference/root_api_usage.md#centroid) |
| `closest_points` | geometry | high | adequate | — | [reference/root_api_usage.md](../reference/root_api_usage.md#closest_points) |
| `connected_components` | geometry | high | adequate | — | [reference/regions.md](../reference/regions.md#primary-api-overview) |
| `contains_point` | geometry | high | adequate | — | [reference/regions.md](../reference/regions.md#primary-api-overview) |
| `coordinate_range` | geometry | high | adequate | — | [reference/root_api_usage.md](../reference/root_api_usage.md#coordinate_range) |
| `covariance_matrix` | integration | high | adequate | — | [reference/root_api_usage.md](../reference/root_api_usage.md#covariance_matrix) |
| `diameter` | geometry | high | adequate | — | [reference/root_api_usage.md](../reference/root_api_usage.md#diameter) |
| `distance_set` | geometry | high | adequate | — | [reference/root_api_usage.md](../reference/root_api_usage.md#distance_set) |
| `extrema_set` | geometry | high | adequate | — | [reference/root_api_usage.md](../reference/root_api_usage.md#extrema_set) |
| `has_empty_interior` | geometry | high | adequate | — | [reference/root_api_usage.md](../reference/root_api_usage.md#has_empty_interior) |
| `inertia_tensor` | geometry | high | adequate | — | [reference/root_api_usage.md](../reference/root_api_usage.md#inertia_tensor) |
| `intersects` | geometry | high | adequate | — | [reference/root_api_usage.md](../reference/root_api_usage.md#intersects) |
| `is_interior_disjoint` | geometry | high | adequate | — | [reference/root_api_usage.md](../reference/root_api_usage.md#is_interior_disjoint) |
| `is_bounded` | geometry | high | adequate | — | [reference/root_api_usage.md](../reference/root_api_usage.md#is_bounded) |
| `is_closed` | geometry | high | adequate | — | [reference/root_api_usage.md](../reference/root_api_usage.md#is_closed) |
| `is_compact` | geometry | high | adequate | — | [reference/root_api_usage.md](../reference/root_api_usage.md#is_compact) |
| `is_connected` | geometry | high | adequate | — | [reference/root_api_usage.md](../reference/root_api_usage.md#is_connected) |
| `is_dense_in` | geometry | high | adequate | — | [reference/root_api_usage.md](../reference/root_api_usage.md#is_dense_in) |
| `is_disjoint` | geometry | high | adequate | — | [reference/root_api_usage.md](../reference/root_api_usage.md#is_disjoint) |
| `is_empty` | geometry | high | adequate | — | [reference/root_api_usage.md](../reference/root_api_usage.md#is_empty) |
| `is_equal` | geometry | high | adequate | — | [reference/root_api_usage.md](../reference/root_api_usage.md#is_equal) |
| `is_full_dimensional` | geometry | high | adequate | — | [reference/root_api_usage.md](../reference/root_api_usage.md#is_full_dimensional) |
| `is_open` | geometry | high | adequate | — | [reference/root_api_usage.md](../reference/root_api_usage.md#is_open) |
| `is_subset` | geometry | high | adequate | — | [reference/root_api_usage.md](../reference/root_api_usage.md#is_subset) |
| `level_set` | geometry | high | adequate | — | [reference/root_api_usage.md](../reference/root_api_usage.md#level_set) |
| `linear_image` | geometry | high | adequate | — | [reference/root_api_usage.md](../reference/root_api_usage.md#linear_image) |
| `minkowski_sum` | geometry | high | adequate | — | [reference/root_api_usage.md](../reference/root_api_usage.md#minkowski_sum) |
| `moment_matrix` | geometry | high | adequate | — | [reference/root_api_usage.md](../reference/root_api_usage.md#moment_matrix) |
| `nearest_point` | geometry | high | adequate | — | [reference/root_api_usage.md](../reference/root_api_usage.md#nearest_point) |
| `scale` | geometry | high | adequate | — | [reference/root_api_usage.md](../reference/root_api_usage.md#scale) |
| `squared_distance_range` | geometry | high | adequate | — | [reference/root_api_usage.md](../reference/root_api_usage.md#squared_distance_range) |
| `sublevel_set` | geometry | high | adequate | — | [reference/root_api_usage.md](../reference/root_api_usage.md#sublevel_set) |
| `superlevel_set` | geometry | high | adequate | — | [reference/root_api_usage.md](../reference/root_api_usage.md#superlevel_set) |
| `support_function` | geometry | high | adequate | — | [reference/root_api_usage.md](../reference/root_api_usage.md#support_function) |
| `translate` | geometry | high | adequate | — | [reference/root_api_usage.md](../reference/root_api_usage.md#translate) |
| `width` | geometry | high | adequate | — | [reference/root_api_usage.md](../reference/root_api_usage.md#width) |
| `function_domain` | decision | high | adequate | — | [reference/decision_and_qe.md](../reference/decision_and_qe.md#primary-api-overview) |
| `is_real_valued` | solving | medium | adequate | — | [reference/root_api_usage.md](../reference/root_api_usage.md#is_real_valued) |
| `bounding_box` | geometry | high | adequate | — | [reference/root_api_usage.md](../reference/root_api_usage.md#bounding_box) |
| `critical_values` | geometry | high | adequate | — | [reference/root_api_usage.md](../reference/root_api_usage.md#critical_values) |
| `critical_value_image` | geometry | high | adequate | — | [reference/root_api_usage.md](../reference/root_api_usage.md#critical_value_image) |
| `distance_between_regions` | geometry | high | adequate | — | [reference/root_api_usage.md](../reference/root_api_usage.md#distance_between_regions) |
| `distance_to_region` | geometry | high | adequate | — | [reference/root_api_usage.md](../reference/root_api_usage.md#distance_to_region) |
| `euler_characteristic` | geometry | high | adequate | — | [reference/regions.md](../reference/regions.md#primary-api-overview) |
| `fiber` | geometry | high | adequate | — | [reference/root_api_usage.md](../reference/root_api_usage.md#fiber) |
| `is_path_connected` | geometry | high | adequate | — | [reference/root_api_usage.md](../reference/root_api_usage.md#is_path_connected) |
| `path_between` | geometry | high | adequate | — | [reference/root_api_usage.md](../reference/root_api_usage.md#path_between) |
| `semialgebraic_projection` | geometry | high | adequate | — | [reference/root_api_usage.md](../reference/root_api_usage.md#semialgebraic_projection) |
| `semialgebraic_measure` | integration | high | adequate | — | [reference/root_api_usage.md](../reference/root_api_usage.md#semialgebraic_measure) |
| `region_measure` | integration | high | adequate | — | [reference/root_api_usage.md](../reference/root_api_usage.md#region_measure) |
| `region_moment` | integration | high | adequate | — | [reference/root_api_usage.md](../reference/root_api_usage.md#region_moment) |
| `function_range` | optimization | high | adequate | — | [reference/root_api_usage.md](../reference/root_api_usage.md#function_range) |
| `semialgebraic_maximize` | optimization | high | adequate | — | [reference/root_api_usage.md](../reference/root_api_usage.md#semialgebraic_maximize) |
| `semialgebraic_minimize` | optimization | high | adequate | — | [reference/root_api_usage.md](../reference/root_api_usage.md#semialgebraic_minimize) |
| `connected_component_count` | geometry | high | adequate | — | [reference/root_api_usage.md](../reference/root_api_usage.md#connected_component_count) |
| `connected_component_samples` | geometry | high | adequate | — | [reference/root_api_usage.md](../reference/root_api_usage.md#connected_component_samples) |
| `topology_summary` | geometry | high | adequate | — | [reference/root_api_usage.md](../reference/root_api_usage.md#topology_summary) |
| `betti_number` | geometry | high | adequate | — | [reference/root_api_usage.md](../reference/root_api_usage.md#betti_number) |
| `solvability_conditions` | decision | critical | adequate | — | [reference/root_api_usage.md](../reference/root_api_usage.md#solvability_conditions) |
| `clip_affine_subspace_to_box` | geometry | high | adequate | — | [reference/root_api_usage.md](../reference/root_api_usage.md#clip_affine_subspace_to_box) |
| `apply_quantifiers` | geometry | critical | adequate | — | [reference/root_api_usage.md](../reference/root_api_usage.md#apply_quantifiers) |
| `prove_negative` | decision | critical | adequate | — | [reference/root_api_usage.md](../reference/root_api_usage.md#prove_negative) |
| `prove_nonnegative` | decision | critical | adequate | — | [reference/root_api_usage.md](../reference/root_api_usage.md#prove_nonnegative) |
| `prove_nonpositive` | decision | critical | adequate | — | [reference/root_api_usage.md](../reference/root_api_usage.md#prove_nonpositive) |
| `prove_positive` | decision | critical | adequate | — | [reference/root_api_usage.md](../reference/root_api_usage.md#prove_positive) |
| `prove_zero` | decision | critical | adequate | — | [reference/root_api_usage.md](../reference/root_api_usage.md#prove_zero) |
| `prove_nonzero` | decision | critical | adequate | — | [reference/root_api_usage.md](../reference/root_api_usage.md#prove_nonzero) |
| `function_sign` | function-analysis | medium | adequate | — | [reference/convexity_backend.md](../reference/convexity_backend.md#family-contract) |
| `simplify_system` | decision | high | adequate | — | [reference/root_api_usage.md](../reference/root_api_usage.md#simplify_system) |
| `simplify_under_assumptions` | decision | high | adequate | — | [reference/decision_and_qe.md](../reference/decision_and_qe.md#primary-api-overview) |
| `local_dimension` | geometry | high | adequate | — | [reference/root_api_usage.md](../reference/root_api_usage.md#local_dimension) |
| `region_boundary_result` | geometry | high | adequate | — | [reference/root_api_usage.md](../reference/root_api_usage.md#region_boundary_result) |
| `region_active_boundary_strata` | geometry | high | adequate | — | [reference/root_api_usage.md](../reference/root_api_usage.md#region_active_boundary_strata) |
| `region_nonsmooth_locus` | geometry | high | adequate | — | [reference/root_api_usage.md](../reference/root_api_usage.md#region_nonsmooth_locus) |
| `region_regular_locus` | geometry | high | adequate | — | [reference/root_api_usage.md](../reference/root_api_usage.md#region_regular_locus) |
| `region_singular_locus` | geometry | high | adequate | — | [reference/root_api_usage.md](../reference/root_api_usage.md#region_singular_locus) |
| `region_singular_locus_result` | geometry | high | adequate | — | [reference/root_api_usage.md](../reference/root_api_usage.md#region_singular_locus_result) |
| `integrate_over_region` | integration | high | adequate | — | [reference/root_api_usage.md](../reference/root_api_usage.md#integrate_over_region) |
| `reduce_region_integral` | integration | high | adequate | — | [reference/root_api_usage.md](../reference/root_api_usage.md#reduce_region_integral) |
| `region_boundary` | geometry | high | adequate | — | [reference/root_api_usage.md](../reference/root_api_usage.md#region_boundary) |
| `region_closure` | geometry | high | adequate | — | [reference/regions.md](../reference/regions.md#primary-api-overview) |
| `region_complement` | geometry | high | adequate | — | [reference/root_api_usage.md](../reference/root_api_usage.md#region_complement) |
| `region_difference` | geometry | high | adequate | — | [reference/root_api_usage.md](../reference/root_api_usage.md#region_difference) |
| `region_dimension` | geometry | high | adequate | — | [reference/regions.md](../reference/regions.md#primary-api-overview) |
| `region_interior` | geometry | high | adequate | — | [reference/root_api_usage.md](../reference/root_api_usage.md#region_interior) |
| `region_intersection` | geometry | high | adequate | — | [reference/regions.md](../reference/regions.md#primary-api-overview) |
| `region_product` | geometry | high | adequate | — | [reference/regions.md](../reference/regions.md#primary-api-overview) |
| `region_symmetric_difference` | geometry | high | adequate | — | [reference/root_api_usage.md](../reference/root_api_usage.md#region_symmetric_difference) |
| `region_union` | geometry | high | adequate | — | [reference/regions.md](../reference/regions.md#primary-api-overview) |
| `classify_real_roots` | geometry | high | adequate | — | [reference/root_api_usage.md](../reference/root_api_usage.md#classify_real_roots) |
| `sample_point` | solving | medium | adequate | — | [reference/solving_and_sampling.md](../reference/solving_and_sampling.md#primary-api-overview) |
| `random_point` | solving | medium | adequate | — | [reference/solving_and_sampling.md](../reference/solving_and_sampling.md#primary-api-overview) |
| `random_points` | solving | medium | adequate | — | [reference/solving_and_sampling.md](../reference/solving_and_sampling.md#primary-api-overview) |
| `affine_relative_interior_formula` | function-analysis | medium | adequate | — | [reference/convexity_backend.md](../reference/convexity_backend.md#family-contract) |
| `strict_feasible` | function-analysis | medium | adequate | — | [reference/convexity_backend.md](../reference/convexity_backend.md#family-contract) |
| `sample_points` | solving | medium | adequate | — | [reference/solving_and_sampling.md](../reference/solving_and_sampling.md#primary-api-overview) |
| `sign_at` | solving | medium | adequate | — | [reference/solving_and_sampling.md](../reference/solving_and_sampling.md#primary-api-overview) |
| `sign_vector` | solving | medium | adequate | — | [reference/solving_and_sampling.md](../reference/solving_and_sampling.md#primary-api-overview) |
| `discretize_region_geometry` | solving | medium | adequate | — | [reference/solving_and_sampling.md](../reference/solving_and_sampling.md#primary-api-overview) |
| `discretize_solution` | solving | medium | adequate | — | [reference/solving_and_sampling.md](../reference/solving_and_sampling.md#primary-api-overview) |
| `plot_region_geometry` | solving | medium | adequate | — | [reference/solving_and_sampling.md](../reference/solving_and_sampling.md#primary-api-overview) |
| `plot_solution` | solving | medium | adequate | — | [reference/solving_and_sampling.md](../reference/solving_and_sampling.md#primary-api-overview) |
| `find_instance` | solving | medium | adequate | — | [reference/solving_and_sampling.md](../reference/solving_and_sampling.md#primary-api-overview) |
| `is_zero_dimensional` | solving | medium | adequate | — | [reference/solving_and_sampling.md](../reference/solving_and_sampling.md#primary-api-overview) |
| `reduce_formula` | solving | critical | adequate | — | [reference/root_api_usage.md](../reference/root_api_usage.md#reduce_formula) |
| `resolve_formula` | solving | critical | adequate | — | [reference/root_api_usage.md](../reference/root_api_usage.md#resolve_formula) |
| `affine_image` | geometry | high | adequate | — | [reference/root_api_usage.md](../reference/root_api_usage.md#affine_image) |
| `affine_preimage` | geometry | high | adequate | — | [reference/root_api_usage.md](../reference/root_api_usage.md#affine_preimage) |
| `region_image` | geometry | high | adequate | — | [reference/root_api_usage.md](../reference/root_api_usage.md#region_image) |
| `region_preimage` | geometry | high | adequate | — | [reference/root_api_usage.md](../reference/root_api_usage.md#region_preimage) |
| `deduplicate_indexed_vertices` | geometry | high | adequate | — | [reference/root_api_usage.md](../reference/root_api_usage.md#deduplicate_indexed_vertices) |
| `polygon_vertices` | geometry | high | adequate | — | [reference/root_api_usage.md](../reference/root_api_usage.md#polygon_vertices) |
| `outer_polygons` | geometry | high | adequate | — | [reference/root_api_usage.md](../reference/root_api_usage.md#outer_polygons) |
| `inner_polygons` | geometry | high | adequate | — | [reference/root_api_usage.md](../reference/root_api_usage.md#inner_polygons) |
| `polyhedron_vertices` | geometry | high | adequate | — | [reference/root_api_usage.md](../reference/root_api_usage.md#polyhedron_vertices) |
| `polyhedron_face_indices` | geometry | high | adequate | — | [reference/root_api_usage.md](../reference/root_api_usage.md#polyhedron_face_indices) |
| `outer_polyhedra` | geometry | high | adequate | — | [reference/root_api_usage.md](../reference/root_api_usage.md#outer_polyhedra) |
| `inner_polyhedra` | geometry | high | adequate | — | [reference/root_api_usage.md](../reference/root_api_usage.md#inner_polyhedra) |
| `convex_hull` | function-analysis | medium | adequate | — | [reference/root_api_usage.md](../reference/root_api_usage.md#convex_hull) |
| `random_polygon` | geometry | high | adequate | — | [reference/root_api_usage.md](../reference/root_api_usage.md#random_polygon) |
| `random_polytope` | geometry | high | adequate | — | [reference/root_api_usage.md](../reference/root_api_usage.md#random_polytope) |
| `subdivide_triangular_faces` | geometry | high | adequate | — | [reference/root_api_usage.md](../reference/root_api_usage.md#subdivide_triangular_faces) |
| `geodesic_refinement` | geometry | high | adequate | — | [reference/root_api_usage.md](../reference/root_api_usage.md#geodesic_refinement) |
| `polyhedral_intersection` | geometry | high | adequate | — | [reference/root_api_usage.md](../reference/root_api_usage.md#polyhedral_intersection) |
| `polyhedral_boolean` | geometry | high | adequate | — | [reference/root_api_usage.md](../reference/root_api_usage.md#polyhedral_boolean) |
| `RegularPolygon` | geometry | high | adequate | — | [reference/root_api_usage.md](../reference/root_api_usage.md#regularpolygon) |
| `Cube` | geometry | high | adequate | — | [reference/root_api_usage.md](../reference/root_api_usage.md#cube) |
| `Tetrahedron` | geometry | high | adequate | — | [reference/root_api_usage.md](../reference/root_api_usage.md#tetrahedron) |
| `Octahedron` | geometry | high | adequate | — | [reference/root_api_usage.md](../reference/root_api_usage.md#octahedron) |
| `Icosahedron` | geometry | high | adequate | — | [reference/root_api_usage.md](../reference/root_api_usage.md#icosahedron) |
| `Dodecahedron` | geometry | high | adequate | — | [reference/root_api_usage.md](../reference/root_api_usage.md#dodecahedron) |
| `Prism` | geometry | high | adequate | — | [reference/root_api_usage.md](../reference/root_api_usage.md#prism) |
| `Pyramid` | geometry | high | adequate | — | [reference/root_api_usage.md](../reference/root_api_usage.md#pyramid) |
| `as_semialgebraic_region` | geometry | high | adequate | — | [reference/root_api_usage.md](../reference/root_api_usage.md#as_semialgebraic_region) |
| `interior_of_closure` | geometry | high | adequate | — | [reference/root_api_usage.md](../reference/root_api_usage.md#interior_of_closure) |
| `closure_of_interior` | geometry | high | adequate | — | [reference/root_api_usage.md](../reference/root_api_usage.md#closure_of_interior) |
| `region_variables` | geometry | high | adequate | — | [reference/root_api_usage.md](../reference/root_api_usage.md#region_variables) |
| `is_regular_closed_region` | geometry | high | adequate | — | [reference/root_api_usage.md](../reference/root_api_usage.md#is_regular_closed_region) |
| `is_regular_open_region` | geometry | high | adequate | — | [reference/root_api_usage.md](../reference/root_api_usage.md#is_regular_open_region) |
| `simplify_region` | geometry | high | adequate | — | [reference/root_api_usage.md](../reference/root_api_usage.md#simplify_region) |
| `simplify_boole` | geometry | high | adequate | — | [reference/root_api_usage.md](../reference/root_api_usage.md#simplify_boole) |
| `simplify_piecewise` | geometry | high | adequate | — | [reference/root_api_usage.md](../reference/root_api_usage.md#simplify_piecewise) |

## Interpretation

A failing dimension is a concrete documentation task, not a claim that the implementation is defective. In particular, a family reference can accurately describe shared certification and complexity while a function still lacks parameter-specific guidance or a worked call. Conversely, adding prose merely to satisfy a counter is discouraged: examples should demonstrate a representative semantic distinction, and parameter text should explain non-obvious meaning rather than repeat the signature.

The generated TOML registry is `docs/reference/root_function_documentation_adequacy.toml`. CI checks that it remains complete and reproducible as the root API changes.
