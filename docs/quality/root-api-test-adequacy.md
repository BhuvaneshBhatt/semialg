# Root API test adequacy

This generated audit evaluates every root-level public function by **semantic test dimensions**, not line coverage or raw test count. Existing explicitly curated multidimensional evidence is combined with directly-called tests whose names identify a semantic contract.

## Risk-adjusted adequacy floor

- **Critical:** nominal coverage, at least three semantic dimensions, at least one independent mathematical detector, and at least one adverse/boundary/conditional detector.
- **High:** nominal coverage plus a meaningful independent, boundary, failure, assumptions, or parameter-sensitive detector.
- **Medium:** nominal coverage plus at least one second semantic dimension.

> This is an adequacy floor, not a proof of correctness. Test names only count when the test directly calls the root API; explicit multidimensional evidence is treated as authoritative.

## Results

**218/218 adequate; 0 need deepening.**

| API | Owner | Risk | Status | Dimensions | Missing |
|---|---|---|---|---|---|
| `Cube` | geometry | high | adequate | metamorphic, nominal, representation | — |
| `Dodecahedron` | geometry | high | adequate | metamorphic, nominal | — |
| `Icosahedron` | geometry | high | adequate | metamorphic, nominal | — |
| `Octahedron` | geometry | high | adequate | metamorphic, nominal | — |
| `Prism` | geometry | high | adequate | metamorphic, nominal | — |
| `Pyramid` | geometry | high | adequate | metamorphic, negative, nominal | — |
| `RegularPolygon` | geometry | high | adequate | invalid-failure, metamorphic, nominal, representation | — |
| `Tetrahedron` | geometry | high | adequate | metamorphic, nominal | — |
| `active_constraints` | algebraic | high | adequate | boundary-degenerate, nominal | — |
| `affine_image` | geometry | high | adequate | boundary-degenerate, certificate, invalid-failure, metamorphic, nominal, property-generated, representation, round-trip | — |
| `affine_preimage` | geometry | high | adequate | boundary-degenerate, invalid-failure, nominal, representation, round-trip | — |
| `affine_relative_interior_formula` | function-analysis | medium | adequate | boundary-degenerate, negative, nominal | — |
| `analyze_affine_map` | geometry | high | adequate | boundary-degenerate, invalid-failure, metamorphic, negative, nominal, parameter-regime, round-trip | — |
| `apply_quantifiers` | geometry | critical | adequate | nominal, parameter-regime, round-trip | — |
| `argmax_set` | geometry | high | adequate | boundary-degenerate, nominal | — |
| `argmin_set` | geometry | high | adequate | boundary-degenerate, metamorphic, nominal | — |
| `as_cad_region` | geometry | high | adequate | independent-oracle, metamorphic, nominal | — |
| `as_semialgebraic_region` | geometry | high | adequate | metamorphic, negative, nominal, property-generated, representation | — |
| `betti_number` | geometry | high | adequate | boundary-degenerate, nominal | — |
| `bounding_box` | geometry | high | adequate | metamorphic, nominal | — |
| `cad` | geometry | critical | adequate | boundary-degenerate, differential, independent-oracle, invalid-failure, metamorphic, negative, nominal, representation, round-trip | — |
| `canonicalize_polygon` | geometry | high | adequate | metamorphic, nominal, property-generated, representation | — |
| `canonicalize_polyhedron` | geometry | high | adequate | certificate, nominal, representation | — |
| `canonicalize_region` | geometry | high | adequate | metamorphic, nominal, representation | — |
| `centroid` | integration | high | adequate | independent-oracle, metamorphic, nominal, property-generated | — |
| `certified_radicalization` | algebraic | high | adequate | certificate, metamorphic, nominal | — |
| `classify_real_roots` | geometry | high | adequate | boundary-degenerate, differential, independent-oracle, invalid-failure, nominal, parameter-regime, property-generated | — |
| `clip_affine_subspace_to_box` | geometry | high | adequate | independent-oracle, nominal | — |
| `closest_points` | geometry | high | adequate | metamorphic, nominal | — |
| `closure_of_interior` | geometry | high | adequate | independent-oracle, metamorphic, nominal | — |
| `component_constraint_descriptions` | geometry | high | adequate | certificate, independent-semantic, nominal | — |
| `connected_component_count` | geometry | high | adequate | boundary-degenerate, nominal | — |
| `connected_component_samples` | geometry | high | adequate | boundary-degenerate, nominal | — |
| `connected_components` | geometry | high | adequate | boundary-degenerate, metamorphic, negative, nominal, property-generated | — |
| `contains_point` | geometry | high | adequate | boundary-degenerate, independent-oracle, nominal | — |
| `convert_region` | geometry | high | adequate | boundary-degenerate, nominal, round-trip | — |
| `convex_hull` | function-analysis | medium | adequate | nominal, representation | — |
| `convexity_certificate` | function-analysis | medium | adequate | boundary-degenerate, certificate, invalid-failure, nominal, property-generated | — |
| `coordinate_range` | geometry | high | adequate | differential, metamorphic, nominal, property-generated | — |
| `covariance_matrix` | integration | high | adequate | metamorphic, nominal | — |
| `critical_value_image` | geometry | high | adequate | boundary-degenerate, nominal | — |
| `critical_values` | geometry | high | adequate | boundary-degenerate, nominal | — |
| `decompose_polytope` | geometry | high | adequate | independent-semantic, nominal | — |
| `deduplicate_indexed_vertices` | geometry | high | adequate | boundary-degenerate, metamorphic, nominal | — |
| `diameter` | geometry | high | adequate | metamorphic, nominal | — |
| `discretize_region_geometry` | solving | medium | adequate | independent-oracle, metamorphic, negative, nominal, round-trip-presentation | — |
| `discretize_solution` | solving | medium | adequate | boundary-degenerate, nominal | — |
| `distance_between_regions` | geometry | high | adequate | metamorphic, nominal | — |
| `distance_set` | geometry | high | adequate | metamorphic, nominal | — |
| `distance_to_region` | geometry | high | adequate | metamorphic, nominal | — |
| `equivalent` | decision | critical | adequate | assumptions, boundary-degenerate, differential, independent-oracle, invalid-failure, metamorphic, negative, nominal, parameter-regime, representation | — |
| `euler_characteristic` | geometry | high | adequate | independent-semantic, nominal | — |
| `extrema_set` | geometry | high | adequate | boundary-degenerate, independent-semantic, nominal | — |
| `fiber` | geometry | high | adequate | differential, independent-oracle, independent-semantic, invalid-failure, nominal | — |
| `find_instance` | solving | medium | adequate | certificate, differential, independent-oracle, nominal | — |
| `find_negative_point` | geometry | critical | adequate | boundary-degenerate, certificate, negative, nominal | — |
| `find_negative_witness_fast` | geometry | critical | adequate | certificate, invalid-failure, negative, nominal, property-generated | — |
| `function_convex_partition` | function-analysis | medium | adequate | nominal, parameter-regime, property-generated | — |
| `function_convexity` | function-analysis | medium | adequate | boundary-degenerate, certificate, independent-oracle, invalid-failure, negative, nominal, parameter-regime, property-generated | — |
| `function_domain` | decision | high | adequate | metamorphic, nominal, representation | — |
| `function_mapping_properties` | function-analysis | medium | adequate | boundary-degenerate, metamorphic, nominal, parameter-regime, property-generated | — |
| `function_monotonic_partition` | function-analysis | medium | adequate | boundary-degenerate, nominal, parameter-regime, property-generated | — |
| `function_monotonicity` | function-analysis | medium | adequate | boundary-degenerate, certificate, differential, negative, nominal, parameter-regime, property-generated | — |
| `function_range` | optimization | high | adequate | boundary-degenerate, metamorphic, nominal, parameter-regime, round-trip | — |
| `function_sign` | function-analysis | medium | adequate | assumptions, nominal, parameter-regime | — |
| `function_sign_partition` | function-analysis | medium | adequate | independent-oracle, nominal, parameter-regime, property-generated | — |
| `function_smoothness` | function-analysis | medium | adequate | metamorphic, nominal, parameter-regime, property-generated, representation | — |
| `geodesic_refinement` | geometry | high | adequate | independent-semantic, nominal | — |
| `has_empty_interior` | geometry | high | adequate | boundary-degenerate, independent-semantic, nominal, property-generated | — |
| `implicitize_polynomial_map` | algebraic | high | adequate | boundary-degenerate, invalid-failure, nominal, parameter-regime | — |
| `implied_polynomial_inequality` | geometry | high | adequate | certificate, metamorphic, nominal | — |
| `implies` | decision | critical | adequate | assumptions, differential, invalid-failure, metamorphic, negative, nominal | — |
| `inertia_tensor` | geometry | high | adequate | independent-semantic, metamorphic, nominal | — |
| `inner_polygons` | geometry | high | adequate | boundary-degenerate, metamorphic, nominal | — |
| `inner_polyhedra` | geometry | high | adequate | boundary-degenerate, independent-semantic, nominal | — |
| `integrate_over_region` | integration | high | adequate | boundary-degenerate, differential, independent-oracle, invalid-failure, metamorphic, negative, nominal, parameter-regime, representation | — |
| `interior_of_closure` | geometry | high | adequate | independent-oracle, metamorphic, nominal | — |
| `intersects` | geometry | high | adequate | independent-oracle, nominal, property-generated | — |
| `irreducible_components` | algebraic | high | adequate | boundary-degenerate, certificate, metamorphic, nominal | — |
| `is_bijective` | function-analysis | medium | adequate | nominal, property-generated | — |
| `is_bounded` | geometry | high | adequate | nominal, property-generated | — |
| `is_closed` | geometry | high | adequate | independent-oracle, metamorphic, nominal, property-generated | — |
| `is_compact` | geometry | high | adequate | boundary-degenerate, certificate, independent-oracle, metamorphic, nominal, property-generated, round-trip | — |
| `is_connected` | geometry | high | adequate | boundary-degenerate, nominal, property-generated | — |
| `is_convex` | function-analysis | medium | adequate | boundary-degenerate, certificate, differential, metamorphic, nominal, property-generated | — |
| `is_dense_in` | geometry | high | adequate | independent-oracle, nominal | — |
| `is_disjoint` | geometry | high | adequate | boundary-degenerate, metamorphic, nominal, property-generated | — |
| `is_empty` | geometry | high | adequate | metamorphic, nominal, property-generated | — |
| `is_equal` | geometry | high | adequate | boundary-degenerate, certificate, differential, metamorphic, nominal, property-generated, round-trip | — |
| `is_full_dimensional` | geometry | high | adequate | independent-oracle, nominal, property-generated | — |
| `is_function_continuous` | function-analysis | medium | adequate | nominal, property-generated | — |
| `is_function_smooth` | function-analysis | medium | adequate | nominal, property-generated | — |
| `is_injective` | function-analysis | medium | adequate | nominal, property-generated | — |
| `is_interior_disjoint` | geometry | high | adequate | boundary-degenerate, independent-oracle, nominal | — |
| `is_open` | geometry | high | adequate | independent-oracle, metamorphic, nominal, property-generated | — |
| `is_path_connected` | geometry | high | adequate | independent-oracle, nominal | — |
| `is_real_valued` | solving | medium | adequate | assumptions, boundary-degenerate, nominal | — |
| `is_regular_closed_region` | geometry | high | adequate | metamorphic, nominal | — |
| `is_regular_open_region` | geometry | high | adequate | metamorphic, nominal | — |
| `is_satisfiable` | decision | critical | adequate | boundary-degenerate, certificate, differential, invalid-failure, metamorphic, negative, nominal, parameter-regime, representation | — |
| `is_singular` | algebraic | high | adequate | boundary-degenerate, metamorphic, nominal | — |
| `is_smooth` | algebraic | high | adequate | boundary-degenerate, metamorphic, nominal | — |
| `is_subset` | geometry | high | adequate | metamorphic, nominal, property-generated | — |
| `is_surjective` | function-analysis | medium | adequate | nominal, property-generated | — |
| `is_tautology` | decision | critical | adequate | differential, invalid-failure, metamorphic, negative, nominal, parameter-regime | — |
| `is_zero_dimensional` | solving | medium | adequate | boundary-degenerate, invalid-failure, metamorphic, nominal | — |
| `level_set` | geometry | high | adequate | boundary-degenerate, independent-semantic, nominal | — |
| `linear_image` | geometry | high | adequate | metamorphic, nominal | — |
| `local_branch_geometry` | algebraic | high | adequate | negative, nominal, representation | — |
| `local_dimension` | geometry | high | adequate | boundary-degenerate, independent-oracle, nominal | — |
| `local_dimension_strata` | algebraic | high | adequate | boundary-degenerate, metamorphic, nominal | — |
| `matrix_definiteness` | function-analysis | medium | adequate | independent-oracle, metamorphic, nominal | — |
| `matrix_pd_on` | function-analysis | medium | adequate | invalid-failure, nominal, parameter-regime | — |
| `matrix_psd_on` | function-analysis | medium | adequate | invalid-failure, nominal, parameter-regime | — |
| `matrix_rank_on` | function-analysis | medium | adequate | nominal, parameter-regime | — |
| `matrix_rank_stratification` | function-analysis | medium | adequate | boundary-degenerate, nominal, parameter-regime | — |
| `minimal_prime_intersections` | algebraic | high | adequate | boundary-degenerate, independent-semantic, nominal | — |
| `minkowski_sum` | geometry | high | adequate | metamorphic, nominal, representation | — |
| `moment_matrix` | geometry | high | adequate | independent-semantic, metamorphic, nominal | — |
| `nearest_point` | geometry | high | adequate | independent-oracle, nominal | — |
| `nonnegative_combination_certificate` | geometry | high | adequate | certificate, negative, nominal | — |
| `outer_polygons` | geometry | high | adequate | boundary-degenerate, metamorphic, nominal | — |
| `outer_polyhedra` | geometry | high | adequate | boundary-degenerate, independent-semantic, nominal | — |
| `parameterization_critical_locus` | geometry | high | adequate | boundary-degenerate, nominal, parameter-regime | — |
| `parameterization_critical_values` | geometry | high | adequate | boundary-degenerate, nominal, parameter-regime | — |
| `parameterization_geometry` | geometry | high | adequate | boundary-degenerate, nominal, parameter-regime | — |
| `parametric_cad` | geometry | high | adequate | boundary-degenerate, certificate, independent-oracle, negative, nominal, parameter-regime | — |
| `path_between` | geometry | high | adequate | independent-oracle, nominal | — |
| `plot_region_geometry` | solving | medium | adequate | independent-oracle, negative, nominal, round-trip-presentation | — |
| `plot_solution` | solving | medium | adequate | nominal, representation, round-trip-presentation | — |
| `polygon_vertices` | geometry | high | adequate | boundary-degenerate, metamorphic, nominal | — |
| `polygonal_region_from_paths` | geometry | high | adequate | nominal, representation | — |
| `polyhedral_boolean` | geometry | high | adequate | metamorphic, nominal | — |
| `polyhedral_intersection` | geometry | high | adequate | metamorphic, nominal | — |
| `polyhedron_face_indices` | geometry | high | adequate | boundary-degenerate, independent-semantic, metamorphic, nominal | — |
| `polyhedron_vertices` | geometry | high | adequate | boundary-degenerate, independent-semantic, nominal | — |
| `polynomial_constraints` | algebraic | high | adequate | boundary-degenerate, metamorphic, nominal | — |
| `polynomial_nonnegative` | geometry | critical | adequate | certificate, differential, negative, nominal, parameter-regime, property-generated | — |
| `project_region` | decision | high | adequate | boundary-degenerate, nominal, parameter-regime | — |
| `prove_negative` | decision | critical | adequate | certificate, negative, nominal | — |
| `prove_nonnegative` | decision | critical | adequate | certificate, negative, nominal, parameter-regime | — |
| `prove_nonpositive` | decision | critical | adequate | certificate, negative, nominal | — |
| `prove_nonzero` | decision | critical | adequate | assumptions, boundary-degenerate, certificate, negative, nominal, parameter-regime | — |
| `prove_positive` | decision | critical | adequate | certificate, metamorphic, nominal | — |
| `prove_zero` | decision | critical | adequate | assumptions, boundary-degenerate, certificate, negative, nominal, parameter-regime | — |
| `quantifier_eliminate` | geometry | high | adequate | metamorphic, nominal, parameter-regime | — |
| `random_point` | solving | medium | adequate | independent-oracle, nominal, property-generated, representation | — |
| `random_points` | solving | medium | adequate | independent-oracle, nominal, property-generated, representation | — |
| `random_polygon` | geometry | high | adequate | metamorphic, nominal, property-generated, representation | — |
| `random_polytope` | geometry | high | adequate | metamorphic, nominal, property-generated | — |
| `real_algebraic_feasibility` | algebraic | critical | adequate | boundary-degenerate, certificate, invalid-failure, negative, nominal | — |
| `reduce_formula` | solving | critical | adequate | certificate, independent-oracle, negative, nominal | — |
| `reduce_region_integral` | integration | high | adequate | nominal, representation | — |
| `reduced_component_singular_loci` | algebraic | high | adequate | boundary-degenerate, negative, nominal | — |
| `redundant_polynomial_inequalities` | geometry | high | adequate | certificate, metamorphic, nominal | — |
| `region_active_boundary_strata` | geometry | high | adequate | boundary-degenerate, independent-oracle, metamorphic, negative, nominal | — |
| `region_boundary` | geometry | high | adequate | boundary-degenerate, metamorphic, nominal, property-generated | — |
| `region_boundary_result` | geometry | high | adequate | boundary-degenerate, metamorphic, nominal | — |
| `region_closure` | geometry | high | adequate | boundary-degenerate, metamorphic, nominal, property-generated | — |
| `region_complement` | geometry | high | adequate | differential, nominal | — |
| `region_difference` | geometry | high | adequate | differential, nominal | — |
| `region_dimension` | geometry | high | adequate | boundary-degenerate, certificate, independent-oracle, metamorphic, nominal, property-generated, round-trip | — |
| `region_image` | geometry | high | adequate | boundary-degenerate, invalid-failure, nominal, parameter-regime, round-trip | — |
| `region_interior` | geometry | high | adequate | boundary-degenerate, metamorphic, nominal, property-generated | — |
| `region_intersection` | geometry | high | adequate | boundary-degenerate, differential, nominal | — |
| `region_measure` | integration | high | adequate | metamorphic, nominal | — |
| `region_moment` | integration | high | adequate | independent-oracle, metamorphic, nominal | — |
| `region_nonsmooth_locus` | geometry | high | adequate | boundary-degenerate, independent-oracle, negative, nominal | — |
| `region_preimage` | geometry | high | adequate | boundary-degenerate, invalid-failure, metamorphic, nominal, round-trip | — |
| `region_product` | geometry | high | adequate | boundary-degenerate, metamorphic, nominal | — |
| `region_regular_locus` | geometry | high | adequate | boundary-degenerate, independent-oracle, metamorphic, nominal | — |
| `region_singular_locus` | geometry | high | adequate | boundary-degenerate, certificate, independent-oracle, negative, nominal | — |
| `region_singular_locus_result` | geometry | high | adequate | boundary-degenerate, certificate, independent-oracle, metamorphic, nominal | — |
| `region_symmetric_difference` | geometry | high | adequate | independent-oracle, independent-semantic, metamorphic, nominal | — |
| `region_union` | geometry | high | adequate | differential, nominal | — |
| `region_variables` | geometry | high | adequate | independent-oracle, independent-semantic, nominal, parameter-regime | — |
| `relative_boundary` | algebraic | high | adequate | boundary-degenerate, metamorphic, nominal | — |
| `relative_interior` | algebraic | high | adequate | boundary-degenerate, metamorphic, nominal | — |
| `replay_certificate` | geometry | critical | adequate | boundary-degenerate, certificate, invalid-failure, metamorphic, nominal, property-generated, round-trip | — |
| `resolve_formula` | solving | critical | adequate | boundary-degenerate, certificate, independent-oracle, negative, nominal | — |
| `sample_point` | solving | medium | adequate | boundary-degenerate, certificate, nominal | — |
| `sample_points` | solving | medium | adequate | certificate, invalid-failure, negative, nominal, property-generated | — |
| `scale` | geometry | high | adequate | boundary-degenerate, metamorphic, nominal, round-trip | — |
| `semialgebraic_maximize` | optimization | high | adequate | boundary-degenerate, independent-oracle, metamorphic, nominal | — |
| `semialgebraic_measure` | integration | high | adequate | boundary-degenerate, independent-oracle, invalid-failure, metamorphic, negative, nominal, parameter-regime, property-generated | — |
| `semialgebraic_minimize` | optimization | high | adequate | boundary-degenerate, certificate, differential, independent-oracle, invalid-failure, metamorphic, nominal, parameter-regime, round-trip | — |
| `semialgebraic_projection` | geometry | high | adequate | invalid-failure, metamorphic, nominal, round-trip | — |
| `semialgebraic_tangent_cone` | algebraic | high | adequate | independent-semantic, nominal | — |
| `sign_at` | solving | medium | adequate | invalid-failure, nominal, parameter-regime | — |
| `sign_vector` | solving | medium | adequate | metamorphic, nominal | — |
| `simplify_boole` | geometry | high | adequate | invalid-failure, metamorphic, negative, nominal, representation | — |
| `simplify_piecewise` | geometry | high | adequate | assumptions, nominal, parameter-regime | — |
| `simplify_region` | geometry | high | adequate | metamorphic, nominal, representation | — |
| `simplify_system` | decision | high | adequate | nominal, parameter-regime | — |
| `simplify_under_assumptions` | decision | high | adequate | assumptions, boundary-degenerate, nominal, parameter-regime | — |
| `singular_locus` | algebraic | high | adequate | boundary-degenerate, nominal, representation | — |
| `solvability_conditions` | decision | critical | adequate | boundary-degenerate, certificate, independent-oracle, negative, nominal, parameter-regime | — |
| `solve_real_algebraic_set` | algebraic | critical | adequate | certificate, metamorphic, nominal | — |
| `solve_semialgebraic` | decision | critical | adequate | boundary-degenerate, differential, invalid-failure, metamorphic, negative, nominal, parameter-regime, property-generated, representation, round-trip | — |
| `squared_distance_range` | geometry | high | adequate | metamorphic, nominal | — |
| `stratified_singular_geometry` | algebraic | high | adequate | boundary-degenerate, negative, nominal | — |
| `strict_feasible` | function-analysis | medium | adequate | boundary-degenerate, nominal | — |
| `subdivide_triangular_faces` | geometry | high | adequate | independent-semantic, nominal | — |
| `sublevel_set` | geometry | high | adequate | boundary-degenerate, independent-semantic, nominal | — |
| `superlevel_set` | geometry | high | adequate | boundary-degenerate, independent-semantic, nominal | — |
| `support_function` | geometry | high | adequate | metamorphic, nominal | — |
| `tangent_cone` | geometry | high | adequate | independent-oracle, negative, nominal, parameter-regime | — |
| `tangent_dimension` | algebraic | high | adequate | boundary-degenerate, differential, independent-oracle, nominal | — |
| `tangent_space` | geometry | high | adequate | independent-oracle, negative, nominal, parameter-regime | — |
| `tetrahedralize_cell` | geometry | high | adequate | independent-oracle, nominal | — |
| `tetrahedralize_cells` | geometry | high | adequate | independent-semantic, nominal | — |
| `topology_summary` | geometry | high | adequate | independent-oracle, metamorphic, nominal | — |
| `translate` | geometry | high | adequate | boundary-degenerate, invalid-failure, metamorphic, nominal, property-generated, round-trip | — |
| `triangulate_polytope` | geometry | high | adequate | independent-oracle, nominal | — |
| `verify_nonnegative_combination_certificate` | geometry | high | adequate | certificate, invalid-failure, negative, nominal | — |
| `width` | geometry | high | adequate | metamorphic, nominal | — |
| `zariski_closure` | algebraic | high | adequate | differential, independent-oracle, nominal | — |
| `zeng_negative_point` | geometry | high | adequate | certificate, invalid-failure, negative, nominal | — |
