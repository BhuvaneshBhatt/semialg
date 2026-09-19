# Everyday API

`semialg` keeps a curated root namespace, and most applications need only the functions on this page. These names form the recommended discovery surface for ordinary exact semialgebraic work. The complete curated root surface is documented separately in [Public API](public_api.md); expert machinery remains available from its defining submodules.

## Solve and decide

- `solve_semialgebraic` — solve a real semialgebraic formula.
- `is_satisfiable` — decide whether a formula has a real solution.
- `find_instance` — return one or more exact satisfying assignments.
- `equivalent`, `implies` — exact logical comparison of formulas.
- `reduce_formula` — reduce quantified or constrained formulas to a simpler exact result.
- `cad` — construct or query a cylindrical algebraic decomposition when explicit CAD access is needed.
- `polynomial_nonnegative_decision`, `polynomial_nonnegative` — decide global polynomial nonnegativity.
- `prove_positive`, `prove_nonnegative` — prove common sign properties under assumptions.

## Functions and maps

- `function_domain` — exact real domain.
- `function_range` — exact range.
- `function_convexity`, `function_monotonicity`, `function_smoothness` — aggregate certified property analyses.
- `function_mapping_properties` — aggregate injectivity/surjectivity/bijectivity analysis.
- `is_injective`, `is_surjective`, `is_bijective` — direct mapping predicates.
- `is_function_continuous`, `is_function_smooth` — direct regularity predicates.
- `semialgebraic_minimize`, `semialgebraic_maximize` — exact constrained extrema.
- `matrix_definiteness` — exact symmetric-matrix definiteness on a domain.

## Regions and geometry

- `as_semialgebraic_region` — obtain a symbolic region value.
- `contains_point`, `is_empty`, `is_bounded`, `is_connected`, `connected_components` — common region queries.
- `region_union`, `region_intersection`, `region_difference`, `region_complement` — Boolean region operations.
- `region_boundary`, `region_dimension` — basic geometric structure.
- `semialgebraic_image`, `semialgebraic_preimage`, `semialgebraic_projection` — exact map and projection operations.
- `region_measure`, `region_centroid`, `integrate_over_region` — exact geometric measurement and integration.
- `sample_point`, `sample_points` — exact representative points.

## Algebraic geometry and parameters

- `is_singular`, `is_smooth`, `tangent_space`, `tangent_cone` — common local algebraic-geometry queries.
- `solvability_conditions` — common parameter conditions for exact solvability.
- `semialg.parameters.root_count_conditions` — specialist parameter stratification by distinct real-root count.

## Certificates

- `replay_certificate` — independently replay a supported exact certificate or structured certified result.

## Return-value convention

Functions with a single obvious mathematical answer return that answer by default. Rich diagnostics, witnesses, backend traces, or certificates are obtained with `return_result=True` where supported. Aggregate analyses remain structured when several mathematical outputs are inseparable. Tri-state predicates use `None` for an exact question that the selected certified method cannot decide; `None` is never silently converted to `False`.
