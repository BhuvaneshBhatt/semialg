# Public API index

This index is the documentation-coverage manifest for every name exported by `semialg.__all__`. Each public export belongs to one primary reference family. The automated documentation test verifies that the index and `semialg.__all__` remain synchronized.

Names may also appear in other guides when they participate in multiple workflows.

## Decision, QE, and formulas

Primary reference: [Decision, QE, and formulas](decision_and_qe.md)

| Public name | Implementation module | Kind | Summary |
|---|---|---|---|
| `apply_quantifiers` | `semialg.quantifiers` | function | Wrap ``formula`` in a prenex quantifier prefix. |
| `computation_context` | `semialg.context` | function | Create/reuse a context for one complete exact solve operation. |
| `equivalent` | `semialg.decision.api` | function | Return whether two semialgebraic formulas define the same real set. |
| `Exists` | `semialg.quantifiers` | class | Existentially quantify one or more variables in a Boolean formula. |
| `ForAll` | `semialg.quantifiers` | class | Universally quantify one or more variables in a Boolean formula. |
| `function_domain` | `semialg.domain_solve` | function | Return exact recognized real-domain constraints for supported expressions. |
| `implies` | `semialg.decision.api` | function | Return whether ``assumptions`` imply ``conclusion`` over the reals. |
| `is_real_valued` | `semialg.domain_solve` | function | Return whether supported domain conditions follow from assumptions. |
| `is_satisfiable` | `semialg.decision.api` | function | Return whether a real semialgebraic formula has a satisfying point. |
| `is_tautology` | `semialg.decision.api` | function | Return whether a real semialgebraic formula is true for all variables. |
| `prove_negative` | `semialg.reasoning` | function | Return whether the expression is certified negative on the stated domain. |
| `prove_nonnegative` | `semialg.reasoning` | function | Return whether the expression is certified nonnegative on the stated domain. |
| `prove_nonpositive` | `semialg.reasoning` | function | Return whether the expression is certified nonpositive on the stated domain. |
| `prove_positive` | `semialg.reasoning` | function | Return whether the expression is certified positive on the stated domain. |
| `reduce_formula` | `semialg.solve.reduce` | function | Reduce a parsed real formula using the selected exact decision strategy. |
| `resolve_formula` | `semialg.solve.resolve` | function | Resolve a parsed formula and return its exact solution representation. |
| `SatisfiabilityResult` | `semialg.decision.solution` | class | Structured result for a real satisfiability query. |
| `SemialgebraicSolution` | `semialg.decision.solution` | class | Structured solution summary for a semialgebraic constraint system. |
| `SemialgOptions` | `semialg.domains` | class | Shared options accepted by high-level semialgebraic APIs. |
| `simplify_boole` | `semialg.symbolic_simplify` | function | Simplify a semialgebraic Boolean formula over the real numbers. |
| `simplify_piecewise` | `semialg.symbolic_simplify` | function | Simplify a Piecewise expression using semialgebraic branch conditions. |
| `simplify_system` | `semialg.reasoning` | function | Simplify a real semialgebraic system with CAD/QE-backed checks. |
| `simplify_under_assumptions` | `semialg.reasoning` | function | Simplify real expressions using provable assumptions. |
| `solve_semialgebraic` | `semialg.decision.api` | function | Reduce, sample, and summarize a semialgebraic system over the reals. |
| `SolveDomain` | `semialg.domains` | class | str(object='') -> str |

## CAD and decomposition

Primary reference: [CAD and decomposition](cad.md)

| Public name | Implementation module | Kind | Summary |
|---|---|---|---|
| `CADOptions` | `semialg.decomposition.cylindrical` | class | Options accepted by :func:`cad`. |
| `CADResult` | `semialg.decomposition.cylindrical` | class | Public result for CAD requests. |
| `cad` | `semialg.decomposition.cylindrical` | function | Compute a cylindrical algebraic decomposition for a real formula. |

## Solving and sampling

Primary reference: [Solving and sampling](solving_and_sampling.md)

| Public name | Implementation module | Kind | Summary |
|---|---|---|---|
| `discretize_region_geometry` | `semialg.solution_geometry` | function | Return lightweight geometry for explicit standard-region objects. |
| `discretize_solution` | `semialg.solution_geometry` | function | Return a small plotting/discretization representation for a solution. |
| `find_instance` | `semialg.solve.find_instance` | function | Find satisfying assignments for a formula. |
| `RealAlgebraicFeasibilityResult` | `semialg.real_algebraic` | class | Certified real-algebraic feasibility result. |
| `real_algebraic_feasibility` | `semialg.real_algebraic` | function | Decide real feasibility of rational polynomial equations by certified ARS reduction. |
| `solve_real_algebraic_set` | `semialg.real_algebraic` | function | Return one exact real point or certified emptiness. |
| `WitnessSearchResult` | `semialg.witness_heuristics` | class | Exact-verified one-sided witness-search result. |
| `find_negative_witness_fast` | `semialg.witness_heuristics` | function | Search odd-degree and rational-line negative witnesses. |
| `is_zero_dimensional` | `semialg.solve.zero_dimensional` | function | Return whether rational polynomial equations define a finite complex set. |
| `plot_region_geometry` | `semialg.solution_geometry` | function | Plot an explicit standard-region object using Matplotlib. |
| `plot_solution` | `semialg.solution_geometry` | function | Plot a 1D/2D solution using Matplotlib when available. |
| `sample_point` | `semialg.sampling` | function | Return one satisfying real sample point for ``formula``, or ``None``. |
| `sample_points` | `semialg.sampling` | function | Return satisfying real sample points for a quantifier-free formula. |
| `sign_at` | `semialg.sampling` | function | Return the sign of a polynomial/expression at a point. |
| `sign_vector` | `semialg.sampling` | function | Return the signs of ``polys`` at ``point`` in input order. |
| `solve_zero_dimensional_system` | `semialg.solve.zero_dimensional` | function | Solve a finite rational polynomial system exactly. |

## Optimization and ranges

Primary reference: [Optimization and ranges](optimization_and_range.md)

| Public name | Implementation module | Kind | Summary |
|---|---|---|---|
| `function_range` | `semialg.optimization` | function | Return a quantifier-free formula describing a real function range. |
| `FunctionRangeResult` | `semialg.optimization_results` | class | Exact range summary for a supported semialgebraic image problem. |
| `OptimizationResult` | `semialg.optimization_results` | class | Exact optimum summary for supported semialgebraic problems. |
| `PolynomialNegativityResult` | `semialg.polynomial_positivity` | class | Certified polynomial-negativity result. |
| `find_negative_point` | `semialg.polynomial_positivity` | function | Return a certified negative point or certified nonnegativity. |
| `polynomial_nonnegative` | `semialg.polynomial_positivity` | function | Decide global nonnegativity on the certified specialized fragment. |
| `zeng_negative_point` | `semialg.polynomial_positivity` | function | Specialized exact critical-value polynomial negativity backend. |
| `root_count_conditions` | `semialg.parameters` | function | Return parameter conditions grouped by distinct real-root count. |
| `semialgebraic_maximize` | `semialg.optimization` | function | Return an exact maximum/supremum for a polynomial semialgebraic problem. |
| `semialgebraic_minimize` | `semialg.optimization` | function | Return an exact minimum/infimum for a polynomial semialgebraic problem. |
| `solvability_conditions` | `semialg.parameters` | function | Return parameter conditions for real solvability of a constraint system. |

## Regions and geometry

Primary reference: [Regions and geometry](regions.md)

| Public name | Implementation module | Kind | Summary |
|---|---|---|---|
| `AffineBoxClip` | `semialg.polyhedral_clipping` | class | Exact vertices of an affine subspace clipped by an axis-aligned box. |
| `AffineMapAnalysis` | `semialg.affine_geometry` | class | Exact algebraic properties of an affine map ``x -> A*x + b``. |
| `analyze_affine_map` | `semialg.affine_geometry` | function | Analyze an affine expression map or an explicit matrix/offset pair. |
| `BoundaryStratum` | `semialg.region_analysis` | class | One exact CAD boundary cell with membership and active-set metadata. |
| `bounded_parametric_cover` | `semialg.parametric_geometry` | function | Return a certified finite bounded parametric cover when one is structural. |
| `clip_affine_subspace_to_box` | `semialg.polyhedral_clipping` | function | Clip a low-dimensional affine subspace to a box without CAD. |
| `ParametricChart` | `semialg.parametric_geometry` | class | One certified map from a bounded parameter domain into a region. |
| `ParametricCover` | `semialg.parametric_geometry` | class | Finite certified cover of a region or bounded region intersection. |
| `ParametricMapDegree` | `semialg.map_degree` | class | Generic complex fiber degree of a rational map in characteristic zero. |
| `parametric_map_degree` | `semialg.map_degree` | function | Return the generic algebraic fiber degree of a rational map. |
| `RegionBoundaryResult` | `semialg.region_analysis` | class | Exact region boundary together with reusable CAD and cell metadata. |
| `SingularLocusResult` | `semialg.region_analysis` | class | Certified singular-locus result with explicit incomplete semantics. |
| `region_boundary_result` | `semialg.region_analysis` | function | Return exact boundary cells, membership status, active residuals, and CAD. |
| `affine_transform` | `semialg.derived_geometry` | function | Return the exact affine image ``A*x + b`` in the original coordinates. |
| `argmax_set` | `semialg.derived_geometry` | function | Return the exact global maximizer set as a semialgebraic formula. |
| `argmin_set` | `semialg.derived_geometry` | function | Return the exact global minimizer set as a semialgebraic formula. |
| `BallRegion` | `semialg.standard_regions` | class | BallRegion(center: 'Sequence[object]', radius: 'object' = 1) |
| `BooleanRegion` | `semialg.standard_regions` | class | BooleanRegion(op: 'str', regions: 'Sequence[StandardRegion]', *, assume_disjoint: 'bool' = False) |
| `bounding_box` | `semialg.geometry_queries` | function | Compute the exact axis-aligned bounding box by coordinate optimization. |
| `BoxRegion` | `semialg.standard_regions` | class | BoxRegion(bounds: 'Sequence[tuple[object, object]]') |
| `CapsuleRegion` | `semialg.standard_regions` | class | CapsuleRegion(start: 'Sequence[object]', end: 'Sequence[object]', radius: 'object' = 1) |
| `centroid` | `semialg.derived_geometry` | function | Return the exact centroid of a measurable semialgebraic region. |
| `closest_points` | `semialg.derived_geometry` | function | Return all exact closest point pairs when the distance is attained. |
| `ConeRegion` | `semialg.standard_regions` | class | ConeRegion(start: 'Sequence[object]', end: 'Sequence[object]', radius: 'object' = 1) |
| `connected_components` | `semialg.derived_geometry` | function | Return exact CAD-connected-component formulas. |
| `convexity_certificate` | `semialg.convexity` | function | Decide convexity through the staged exact certificate hierarchy. |
| `contains_point` | `semialg.derived_geometry` | function | Return whether an exact point belongs to the semialgebraic region. |
| `coordinate_range` | `semialg.derived_geometry` | function | Return the exact range of one coordinate over a region. |
| `covariance_matrix` | `semialg.derived_geometry` | function | Return the exact covariance matrix of the uniform measure on a region. |
| `critical_values` | `semialg.geometry_queries` | function | Return exact objective values from isolated and constant KKT components. |
| `CylinderRegion` | `semialg.standard_regions` | class | CylinderRegion(start: 'Sequence[object]', end: 'Sequence[object]', radius: 'object' = 1) |
| `diameter` | `semialg.derived_geometry` | function | Return the exact Euclidean diameter (supremal pairwise distance). |
| `distance_between_regions` | `semialg.geometry_queries` | function | Compute exact Euclidean distance between two semialgebraic regions. |
| `distance_set` | `semialg.derived_geometry` | function | Return the exact set of Euclidean pairwise distances as a formula. |
| `distance_to_region` | `semialg.geometry_queries` | function | Compute exact Euclidean distance from a point to a semialgebraic region. |
| `euler_characteristic` | `semialg.geometry_queries` | function | Compute the semialgebraic Euler characteristic from selected CAD cells. |
| `extrema_set` | `semialg.derived_geometry` | function | Return the union of the exact global minimum and maximum sets. |
| `fiber` | `semialg.geometry_queries` | function | Specialize a semialgebraic family at fixed parameter/coordinate values. |
| `has_empty_interior` | `semialg.derived_geometry` | function | Return whether the region has empty ambient interior. |
| `inertia_tensor` | `semialg.derived_geometry` | function | Return the unit-density second moment-of-inertia tensor about the origin. |
| `intersects` | `semialg.derived_geometry` | function | Return whether two semialgebraic regions have nonempty intersection. |
| `IntervalRegion` | `semialg.standard_regions` | class | IntervalRegion(lower: 'object', upper: 'object', *, lower_closed: 'bool' = True, upper_closed: 'bool' = True) |
| `is_bounded` | `semialg.derived_geometry` | function | Return whether the semialgebraic region is bounded. |
| `is_closed` | `semialg.derived_geometry` | function | Return whether the semialgebraic region is closed. |
| `is_compact` | `semialg.derived_geometry` | function | Return whether the semialgebraic region is compact. |
| `is_connected` | `semialg.derived_geometry` | function | Decide connectedness; for semialgebraic sets this equals path connectedness. |
| `is_convex` | `semialg.convexity` | function | Decide semialgebraic set convexity through staged exact certificates and complete QE fallback. |
| `is_dense_in` | `semialg.derived_geometry` | function | Return whether ``subset`` is dense in ``ambient`` in the ambient Euclidean topology. |
| `is_disjoint` | `semialg.derived_geometry` | function | Return whether two semialgebraic regions are disjoint. |
| `is_empty` | `semialg.derived_geometry` | function | Return whether the semialgebraic region is empty. |
| `is_equal` | `semialg.derived_geometry` | function | Return whether two semialgebraic regions define the same set. |
| `is_full_dimensional` | `semialg.derived_geometry` | function | Return whether the region has full dimension in its ambient variables. |
| `is_open` | `semialg.derived_geometry` | function | Return whether the semialgebraic region is open. |
| `is_path_connected` | `semialg.geometry_queries` | function | Decide path connectedness via exact CAD connectivity. |
| `is_singular` | `semialg.algebraic_geometry` | function | Return whether ``point`` is singular on the polynomial variety. |
| `is_smooth` | `semialg.algebraic_geometry` | function | Return whether the real polynomial variety has empty singular locus. |
| `is_subset` | `semialg.derived_geometry` | function | Return whether one semialgebraic region is contained in another. |
| `level_set` | `semialg.derived_geometry` | function | Return ``region ∩ {expression = value}``. |
| `linear_image` | `semialg.derived_geometry` | function | Return the exact linear image ``A*x`` in the original coordinates. |
| `minkowski_sum` | `semialg.derived_geometry` | function | Return the exact Minkowski sum of two semialgebraic regions. |
| `moment_matrix` | `semialg.derived_geometry` | function | Return the normalized raw second-moment matrix ``E[x x.T]``. |
| `nearest_point` | `semialg.derived_geometry` | function | Return all exact nearest points when the distance is attained. |
| `ParallelepipedRegion` | `semialg.standard_regions` | class | ParallelepipedRegion(origin: 'Sequence[object]', vectors: 'Sequence[Sequence[object]]') |
| `ParallelogramRegion` | `semialg.standard_regions` | class | ParallelogramRegion(origin: 'Sequence[object]', vectors: 'Sequence[Sequence[object]]') |
| `ParametricRegion` | `semialg.standard_regions` | class | ParametricRegion(parameters: 'Sequence[sp.Symbol \| str]', limits: 'Sequence[tuple[sp.Symbol \| str, object, object]]', mapping: 'Sequence[object]', *, multiplicity: 'object' = 1, assumptions: 'object' = True) |
| `path_between` | `semialg.geometry_queries` | function | Return a certified CAD cell-chain connecting two points in a region. |
| `PointRegion` | `semialg.standard_regions` | class | PointRegion(points: 'Sequence[Sequence[object]] \| Sequence[object]') |
| `PolygonRegion` | `semialg.standard_regions` | class | PolygonRegion(vertices: 'Sequence[Sequence[object]]') |
| `PolyhedronRegion` | `semialg.standard_regions` | class | PolyhedronRegion(tetrahedra: 'Sequence[TetrahedronRegion \| Sequence[Sequence[object]]]') |
| `PrismRegion` | `semialg.standard_regions` | class | PrismRegion(base: 'PolygonRegion \| SimplexRegion \| Sequence[Sequence[object]]', vector: 'Sequence[object]') |
| `PyramidRegion` | `semialg.standard_regions` | class | PyramidRegion(base: 'PolygonRegion \| SimplexRegion \| Sequence[Sequence[object]]', apex: 'Sequence[object]') |
| `region_boundary` | `semialg.regions.operations` | function | Return the Euclidean boundary of a semialgebraic region. |
| `region_closure` | `semialg.regions.operations` | function | Return the Euclidean closure of a semialgebraic region. |
| `region_complement` | `semialg.regions.operations` | function | Return the complement of an implicit semialgebraic region. |
| `region_components` | `semialg.regions.operations` | function | Return connected-component formulas for simple explicit cases. |
| `region_difference` | `semialg.regions.operations` | function | Return ``lhs`` minus ``rhs`` for implicit semialgebraic regions. |
| `region_dimension` | `semialg.regions.operations` | function | Return the exact semialgebraic dimension from a complete adapted CAD. |
| `region_interior` | `semialg.regions.operations` | function | Return the Euclidean interior of a semialgebraic region. |
| `region_intersection` | `semialg.regions.operations` | function | Return the intersection of implicit semialgebraic regions. |
| `region_union` | `semialg.regions.operations` | function | Return the union of implicit semialgebraic regions. |
| `RegionDifference` | `semialg.standard_regions` | function | Return the Boolean difference of two standard regions. |
| `RegionIntersection` | `semialg.standard_regions` | function | Return the Boolean intersection of standard regions. |
| `RegionSymmetricDifference` | `semialg.standard_regions` | function | Return the Boolean symmetric difference of two standard regions. |
| `RegionUnion` | `semialg.standard_regions` | function | Return the Boolean union of standard regions. |
| `scale` | `semialg.derived_geometry` | function | Scale a region about the origin by a scalar factor. |
| `semialgebraic_image` | `semialg.geometry_queries` | function | Return the exact semialgebraic image of a polynomial/rational map. |
| `semialgebraic_preimage` | `semialg.geometry_queries` | function | Return the preimage of a semialgebraic target under a symbolic map. |
| `semialgebraic_projection` | `semialg.geometry_queries` | function | Project ``region`` by existentially eliminating the requested variables. |
| `SimplexRegion` | `semialg.standard_regions` | class | SimplexRegion(vertices: 'Sequence[Sequence[object]]') |
| `singular_locus` | `semialg.algebraic_geometry` | function | Return equations defining the singular locus of an algebraic variety. |
| `SphereRegion` | `semialg.standard_regions` | class | SphereRegion(center: 'Sequence[object]', radius: 'object' = 1) |
| `SphericalShellRegion` | `semialg.standard_regions` | class | SphericalShellRegion(center: 'Sequence[object]', radii: 'tuple[object, object]') |
| `squared_distance_range` | `semialg.derived_geometry` | function | Return the exact range of squared pairwise distances. |
| `StadiumRegion` | `semialg.standard_regions` | class | StadiumRegion(start: 'Sequence[object]', end: 'Sequence[object]', radius: 'object' = 1) |
| `StandardRegion` | `semialg.standard_regions` | class | Base class for explicit region objects supported by semialg. |
| `sublevel_set` | `semialg.derived_geometry` | function | Return ``region ∩ {expression <= value}`` (or strict variant). |
| `superlevel_set` | `semialg.derived_geometry` | function | Return ``region ∩ {expression >= value}`` (or strict variant). |
| `support_function` | `semialg.derived_geometry` | function | Return ``sup(x·direction)`` over the region. |
| `tangent_cone` | `semialg.algebraic_geometry` | function | Return the exact ideal-theoretic Zariski tangent cone at ``point``. |
| `tangent_dimension` | `semialg.algebraic_geometry` | function | Return the exact Zariski tangent-space dimension at ``point``. |
| `tangent_space` | `semialg.algebraic_geometry` | function | Return the Zariski tangent space at a point as the Jacobian nullspace. |
| `TetrahedronRegion` | `semialg.standard_regions` | class | TetrahedronRegion(vertices: 'Sequence[Sequence[object]]') |
| `TransformedRegion` | `semialg.standard_regions` | class | TransformedRegion(base: 'StandardRegion', mapping: 'Sequence[object]', base_variables: 'Sequence[sp.Symbol \| str]') |
| `translate` | `semialg.derived_geometry` | function | Translate a region by ``vector`` while preserving coordinate symbols. |
| `width` | `semialg.derived_geometry` | function | Return exact directional width ``max u·x - min u·x``. |

## Integration, measure, and moments

Primary reference: [Integration, measure, and moments](integration_and_moments.md)

| Public name | Implementation module | Kind | Summary |
|---|---|---|---|
| `integrate_over_region` | `semialg.region_integrate` | function | Integrate ``integrand`` over a supported semialgebraic region. |
| `reduce_region_integral` | `semialg.region_integrate` | function | Reduce a supported region integral to explicit iterated integrals. |
| `region_centroid` | `semialg.moments` | function | Return the centroid of a finite-measure semialgebraic region. |
| `region_covariance` | `semialg.moments` | function | Return the covariance matrix of the uniform measure on a region. |
| `region_moment` | `semialg.moments` | function | Return a raw moment integral over a semialgebraic region. |
| `semialgebraic_measure` | `semialg.measure` | function | Return the exact measure of a supported semialgebraic set. |

## Algebraic computation

Primary reference: [Algebraic computation](algebraic.md)

| Public name | Implementation module | Kind | Summary |
|---|---|---|---|
| `classify_real_roots` | `semialg.root_classification` | function | Classify real roots of a univariate polynomial or polynomial family. |

## Parameters and conditional results

Primary reference: [Parameters and conditional results](parameters.md)

| Public name | Implementation module | Kind | Summary |
|---|---|---|---|

## Errors

Primary reference: [Errors](../guides/errors_and_failure_modes.md)

| Public name | Implementation module | Kind | Summary |
|---|---|---|---|
| `ResourceLimitError` | `semialg.errors` | class | A computation stopped because a configured resource limit was reached. |
| `SemialgError` | `semialg.errors` | class | Base class for semialg-specific failures. |
| `SemialgStrategyFailure` | `semialg.errors` | class | A speculative exact/symbolic strategy could not handle the input. |
| `UnsupportedFragmentError` | `semialg.errors` | class | The input is valid, but outside the symbolic fragment a strategy supports. |

## Package metadata

Primary reference: [Package metadata](package_metadata.md)

| Public name | Implementation module | Kind | Summary |
|---|---|---|---|
| `__version__` | `builtins` | str | str(object='') -> str |



## Unified regions, diagnostics, and scalability APIs

Primary references: [Regions](regions.md), [CAD](cad.md), and [Exactness and certification](../concepts/exactness_and_certification.md)

| Public name | Implementation module | Kind | Summary |
|---|---|---|---|
| `SemialgebraicContext` | `semialg.context` | class | Reusable normalized semialgebraic problem plus exact computation cache. |
| `SemialgebraicRegion` | `semialg.symbolic_regions` | class | A symbolic semialgebraic subset of ````R^n```` with lazy reusable state. |
| `as_semialgebraic_region` | `semialg.symbolic_regions` | function | Coerce a formula or explicit ````StandardRegion```` to ````SemialgebraicRegion````. |
| `local_dimension` | `semialg.region_analysis` | function | Exact local semialgebraic dimension at a point. |
| `region_active_boundary_strata` | `semialg.region_analysis` | function | Return exact pairwise-disjoint active-inequality boundary strata. |
| `region_nonsmooth_locus` | `semialg.region_analysis` | function | Return a conservative exact nonsmooth/corner locus of a region boundary. |
| `region_regular_locus` | `semialg.region_analysis` | function | Return the part of ````region```` outside its algebraic boundary singular locus. |
| `region_singular_locus` | `semialg.region_analysis` | function | Return an exact singular-locus formula or raise when certification is incomplete. |
| `region_singular_locus_result` | `semialg.region_analysis` | function | Return singular-locus geometry with completeness and diagnostic metadata. |
| `replay_certificate` | `semialg.certificates` | function | Replay a supported exact certificate using an independent public path. |
| `result_diagnostics` | `semialg.certificates` | function | Return one stable diagnostic schema for CAD/QE/optimization certificates. |

## Unified symbolic/CAD region facade additions

These exports provide the reusable symbolic-region and CAD-region API.

| Public name | Implementation module | Kind | Summary |
|---|---|---|---|
| `CADRegion` | `semialg.cad_region` | class | A reusable region represented by a public :class:`CADResult`. |
| `as_cad_region` | `semialg.cad_region` | function |  |
| `region_closure_interior` | `semialg.symbolic_regions` | function | Return interior(closure(region)). |
| `region_interior_closure` | `semialg.symbolic_regions` | function | Return closure(interior(region)). |
| `region_variables` | `semialg.symbolic_regions` | function | Return coordinate variables, parameters, or all symbols of a region. |
| `is_regular_closed_region` | `semialg.symbolic_regions` | function |  |
| `is_regular_open_region` | `semialg.symbolic_regions` | function |  |
| `simplify_region` | `semialg.symbolic_regions` | function | Canonicalize a symbolic semialgebraic region formula. |

## Convexity backend primitives

Primary reference: [Convexity backend primitives](convexity_backend.md)

| Public name | Implementation module | Kind | Summary |
|---|---|---|---|
| `affine_relative_interior_formula` | `semialg.strict_feasibility` | function | Construct the relative-interior formula of a conjunctive affine system. |
| `function_convexity` | `semialg.function_analysis` | function | Classify exact convexity/concavity of a supported semialgebraic function. |
| `function_convex_partition` | `semialg.function_analysis` | function | Partition a univariate semialgebraic domain into exact convex/concave regions. |
| `function_monotonicity` | `semialg.function_analysis` | function | Classify exact monotonicity of a supported univariate semialgebraic function. |
| `function_monotonic_partition` | `semialg.function_analysis` | function | Partition a univariate semialgebraic domain into exact monotonicity regions. |
| `matrix_definiteness` | `semialg.matrix_analysis` | function | Decide exact symmetric-matrix definiteness on a semialgebraic domain. |
| `matrix_pd_on` | `semialg.matrix_analysis` | function | Decide positive definiteness on a semialgebraic domain. |
| `matrix_psd_on` | `semialg.matrix_analysis` | function | Decide positive semidefiniteness on a semialgebraic domain. |
| `matrix_rank_on` | `semialg.matrix_analysis` | function | Return the exact constant matrix rank on a region, or `None` when it varies. |
| `matrix_rank_stratification` | `semialg.matrix_analysis` | function | Partition parameter space by exact determinantal rank. |
| `parametric_affine_reduction` | `semialg.affine_reduction` | function | Reduce affine equalities with explicit zero/nonzero parameter-pivot branches. |
| `prove_nonzero` | `semialg.reasoning` | function | Prove that an expression is everywhere nonzero on a semialgebraic region. |
| `prove_zero` | `semialg.reasoning` | function | Prove that an expression is identically zero on a semialgebraic region. |
| `function_sign` | `semialg.reasoning` | function | Return a canonical exact sign classification, optionally stratified by parameters. |
| `function_sign_partition` | `semialg.function_analysis` | function | Partition a univariate semialgebraic function domain into exact positive, zero, and negative regions. |
| `function_smoothness` | `semialg.function_properties` | function | Report continuity, differentiability order, smoothness, and exact exceptional loci. |
| `function_mapping_properties` | `semialg.function_properties` | function | Certify injectivity, surjectivity, bijectivity, image, and exact witnesses. |
| `strict_feasible` | `semialg.strict_feasibility` | function | Decide affine-relative or ambient strict feasibility and return a witness in result mode. |
