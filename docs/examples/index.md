# Worked example gallery

These examples are designed as complete, executable workflows. Each page has a matching script under `examples/gallery/`, and CI executes the scripts so code and documented behavior stay synchronized.

The gallery favors small exact problems that expose the real algorithms—CAD/QE, exact optimization, integration, algebraic geometry, and parameter stratification—without making the default documentation checks prohibitively expensive.

- [Projection as quantifier elimination](01_projection_and_qe.md) — Eliminate a coordinate from a parabolic strip and recover the exact projected interval.
- [Exact function range on a disk](02_exact_function_range.md) — Compute every attainable value of a linear polynomial on the unit disk.
- [Certified global polynomial optimization](03_certified_nonconvex_optimization.md) — Minimize a polynomial over a curved compact region and inspect attainment and certification.
- [An optimizer locus with positive dimension](04_positive_dimensional_argmin.md) — Return the entire minimizer set rather than a single witness.
- [Connected components and Euler characteristic](05_topology_and_components.md) — Use CAD semantics to recover components and basic exact topology.
- [Exact measure, centroid, covariance, and inertia](06_measure_centroid_and_moments.md) — Compute several geometric statistics of the unit disk from exact region integrals.
- [Singular locus, tangent space, and exact tangent cone](07_singular_cusp_geometry.md) — Analyze the cusp y²=x³ at its singular point.
- [Minkowski sums and support functions](08_convex_geometry_operations.md) — Combine exact set operations with convex-geometric queries.
- [Distances and exact bounding boxes](09_distance_and_bounding_geometry.md) — Derive metric and coordinate bounds from semialgebraic sets.
- [One parameter, four exact computations](10_parametric_interval_workflow.md) — Follow the family 0≤x≤a through range, optimization, measure, and integration.
- [Parametric integration with algebraic endpoints](11_parametric_algebraic_endpoints.md) — Integrate over x²≤a while CAD controls the parameter-dependent root branches.
- [Polynomial images and transformed regions](12_images_and_linear_maps.md) — Compute the exact image of a disk under an invertible linear map.
- [Direct mathematical return values](13_direct_mathematical_returns.md) — Use formula/value/witness defaults while structured result objects remain opt-in metadata.
- [Real roots and exact transcendental algebraization](14_real_roots_and_exact_algebraization.md) — Distinguish principal powers from explicit real roots and compute exact trig/exponential ranges.
- [Canonical and stable formula simplification](15_canonical_formula_simplification.md) — Normalize equivalent polynomial atoms and verify deterministic idempotent output.
- [CAD-driven integration variable ordering](16_cad_integration_variable_order.md) — Let CAD choose simpler cylindrical coordinates and exact iterated bounds.

- [Symbolic region conditions and active-boundary strata](17_symbolic_region_conditions.md) — Keep region membership symbolic by default, opt into QE, expose real parameter assumptions, and inspect active boundary strata.
