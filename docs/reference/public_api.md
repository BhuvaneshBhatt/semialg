# Public API index

This index is the documentation-coverage manifest for every name exported by `semialg.__all__`. Each public export belongs to one primary reference family. The automated documentation test verifies that the index and `semialg.__all__` remain synchronized.

Names may also appear in other guides when they participate in multiple workflows.

## Decision, QE, and formulas

Primary reference: [Decision, QE, and formulas](decision_and_qe.md)

| Public name | Implementation module | Kind | Summary |
|---|---|---|---|
| `apply_quantifiers` | `semialg.quantifiers` | function | Wrap ``formula`` in a prenex quantifier prefix. |
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
| `QuantifierEliminationResult` | `semialg.qe.api` | class | Certified real quantifier-elimination result with dispatcher provenance. |
| `quantifier_eliminate` | `semialg.qe.api` | function | Eliminate real quantifiers with a certified specialist-first dispatcher. |
| `project_region` | `semialg.qe.api` | function | Return the exact existential projection of a real semialgebraic region. |
| `resolve_formula` | `semialg.solve.resolve` | function | Resolve a parsed formula and return its exact solution representation. |
| `simplify_boole` | `semialg.symbolic_simplify` | function | Simplify a semialgebraic Boolean formula over the real numbers. |
| `simplify_piecewise` | `semialg.symbolic_simplify` | function | Simplify a Piecewise expression using semialgebraic branch conditions. |
| `simplify_system` | `semialg.reasoning` | function | Simplify a real semialgebraic system with CAD/QE-backed checks. |
| `simplify_under_assumptions` | `semialg.reasoning` | function | Simplify real expressions using provable assumptions. |
| `solve_semialgebraic` | `semialg.decision.api` | function | Reduce, sample, and summarize a semialgebraic system over the reals. |

## CAD and decomposition

Primary reference: [CAD and decomposition](cad.md)

| Public name | Implementation module | Kind | Summary |
|---|---|---|---|
| `cad` | `semialg.decomposition.cylindrical` | function | Compute a cylindrical algebraic decomposition for a real formula. |
| `parametric_cad` | `semialg.decomposition.parametric` | function | Compute a parameter-aware cylindrical decomposition with generic, exceptional, and specialized-fiber strata. |

## Solving and sampling

Primary reference: [Solving and sampling](solving_and_sampling.md)

| Public name | Implementation module | Kind | Summary |
|---|---|---|---|
| `discretize_region_geometry` | `semialg.solution_geometry` | function | Return lightweight geometry for explicit standard-region objects. |
| `discretize_solution` | `semialg.solution_geometry` | function | Return a small plotting/discretization representation for a solution. |
| `find_instance` | `semialg.solve.find_instance` | function | Find satisfying assignments for a formula. |
| `real_algebraic_feasibility` | `semialg.real_algebraic` | function | Decide real feasibility of rational polynomial equations by certified ARS reduction. |
| `solve_real_algebraic_set` | `semialg.real_algebraic` | function | Return one exact real point or certified emptiness. |
| `find_negative_witness_fast` | `semialg.witness_heuristics` | function | Search odd-degree and rational-line negative witnesses. |
| `is_zero_dimensional` | `semialg.solve.zero_dimensional` | function | Return whether rational polynomial equations define a finite complex set. |
| `plot_region_geometry` | `semialg.solution_geometry` | function | Plot an explicit standard-region object using Matplotlib. |
| `plot_solution` | `semialg.solution_geometry` | function | Plot a 1D/2D solution using Matplotlib when available. |
| `sample_point` | `semialg.sampling` | function | Return one satisfying real sample point for ``formula``, or ``None``. |
| `sample_points` | `semialg.sampling` | function | Return satisfying real sample points for a quantifier-free formula. |
| `random_point` | `semialg.region_sampling` | function | Return one random point, intrinsically uniform for supported canonical regions. |
| `random_points` | `semialg.region_sampling` | function | Return seeded random points, intrinsically uniform for supported canonical regions. |
| `sign_at` | `semialg.sampling` | function | Return the sign of a polynomial/expression at a point. |
| `sign_vector` | `semialg.sampling` | function | Return the signs of ``polys`` at ``point`` in input order. |

## Optimization and ranges

Primary reference: [Optimization and ranges](optimization_and_range.md)

| Public name | Implementation module | Kind | Summary |
|---|---|---|---|
| `function_range` | `semialg.optimization` | function | Return a quantifier-free formula describing a real function range. |
| `find_negative_point` | `semialg.polynomial_positivity` | function | Return a certified negative point or certified nonnegativity. |
| `polynomial_nonnegative` | `semialg.polynomial_positivity` | function | Decide certified global nonnegativity; use `return_result=True` for the SOS/Zeng/ARS/CAD portfolio trace. |
| `zariski_closure` | `semialg.algebraic_geometry` | function | Return the Zariski closure of a polynomial-map image. |
| `zeng_negative_point` | `semialg.polynomial_positivity` | function | Specialized exact critical-value polynomial negativity backend. |
| `semialgebraic_maximize` | `semialg.optimization` | function | Return an exact maximum/supremum for a polynomial semialgebraic problem. |
| `semialgebraic_minimize` | `semialg.optimization` | function | Return an exact minimum/infimum for a polynomial semialgebraic problem. |
| `solvability_conditions` | `semialg.parameters` | function | Return parameter conditions for real solvability of a constraint system. |

## Regions and geometry

Primary reference: [Regions and geometry](regions.md)

| Public name | Implementation module | Kind | Summary |
|---|---|---|---|
| `AffineBoxClip` | `semialg.polyhedral_clipping` | class | Exact vertices of an affine subspace clipped by an axis-aligned box. |
| `analyze_affine_map` | `semialg.affine_geometry` | function | Analyze an affine expression map or an explicit matrix/offset pair. |
| `clip_affine_subspace_to_box` | `semialg.polyhedral_clipping` | function | Clip a low-dimensional affine subspace to a box without CAD. |
| `region_boundary_result` | `semialg.region_analysis` | function | Return exact boundary cells, membership status, active residuals, and CAD. |
| `argmax_set` | `semialg.derived_geometry` | function | Return the exact global maximizer set as a semialgebraic formula. |
| `argmin_set` | `semialg.derived_geometry` | function | Return the exact global minimizer set as a semialgebraic formula. |
| `AffineHalfSpace` | `semialg.standard_regions` | class | Half-space intrinsic to an affine subspace. |
| `AffineSpace` | `semialg.standard_regions` | class | Affine space represented by one point and independent directions. |
| `Geometry` | `semialg.standard_regions` | class | Common structural geometry protocol with semialgebraic formula lowering. |
| `HalfSpace` | `semialg.standard_regions` | class | Closed half-space defined by an oriented hyperplane. |
| `Hyperplane` | `semialg.standard_regions` | class | Codimension-one affine space defined by a normal and point. |
| `Line` | `semialg.standard_regions` | class | Infinite affine line represented by a point and direction. |
| `HRepresentation` | `semialg.polyhedral` | class | Exact closed halfspace representation `A x <= b` with bounded H→V conversion and redundancy certification. |
| `affine_image` | `semialg.region_transformations` | function | Return the exact affine image, dispatching between structured regions and formula-level sets. |
| `affine_preimage` | `semialg.region_transformations` | function | Return the exact affine preimage, preserving structured representations when possible. |
| `region_image` | `semialg.region_transformations` | function | Return the exact image under a symbolic map with representation-aware dispatch. |
| `region_preimage` | `semialg.region_transformations` | function | Return the exact preimage under a symbolic map with representation-aware dispatch. |
| `Point` | `semialg.standard_regions` | class | Canonical single-point geometry. |
| `Ray` | `semialg.standard_regions` | class | Closed affine ray represented by a point and direction. |
| `Polytope` | `semialg.standard_regions` | class | Canonical convex polytope represented by vertices. |
| `Simplex` | `semialg.standard_regions` | class | Canonical simplex with affinely independent vertices. |
| `Triangle` | `semialg.standard_regions` | class | Convenience constructor namespace returning canonical simplices. |
| `Polygon` | `semialg.standard_regions` | class | Canonical simple polygon in two dimensions. |
| `RegularPolygon` | `semialg.standard_regions` | function | Construct a regular polygon as a canonical polygon. |
| `PolyhedralCone` | `semialg.standard_regions` | class | Affine polyhedral cone with lineality directions and nonnegative generating rays. |
| `Parallelepiped` | `semialg.standard_regions` | class | Canonical affine image of a unit box. |
| `Cube` | `semialg.standard_regions` | function | Construct an axis-aligned cube as a canonical parallelepiped. |
| `Tetrahedron` | `semialg.standard_regions` | function | Construct an explicit or regular tetrahedron as a canonical simplex. |
| `Octahedron` | `semialg.standard_regions` | function | Construct a regular octahedron as a canonical polytope. |
| `Icosahedron` | `semialg.standard_regions` | function | Construct a regular icosahedron as a canonical polytope. |
| `Dodecahedron` | `semialg.standard_regions` | function | Construct a regular dodecahedron as a canonical polytope. |
| `Prism` | `semialg.standard_regions` | function | Extrude a vertex-defined base into a canonical polytope. |
| `Pyramid` | `semialg.standard_regions` | function | Join a base to an apex as a canonical polytope. |
| `Hexahedron` | `semialg.standard_regions` | class | Convex three-dimensional polytope with eight vertices and six quadrilateral facets. |
| `PolygonalComponent` | `semialg.boundary_topology` | class | One connected polygonal component with an outer boundary and explicit holes. |
| `PolygonalSet` | `semialg.boundary_topology` | class | Canonical boundary topology for disconnected planar polygonal sets with holes. |
| `PolyhedralShell` | `semialg.boundary_topology` | class | Connected closed consistently oriented polygonal boundary shell. |
| `PolyhedralComponent` | `semialg.boundary_topology` | class | One solid component with an outer shell and explicit cavity shells. |
| `Polyhedron` | `semialg.boundary_topology` | class | Boundary representation of one or more three-dimensional solid components. |
| `deduplicate_indexed_vertices` | `semialg.boundary_topology` | function | Deduplicate exact indexed vertices and remap cell indices. |
| `polygon_vertices` | `semialg.boundary_topology` | function | Return unique vertices from polygonal boundary topology. |
| `outer_polygons` | `semialg.boundary_topology` | function | Return outer polygon boundaries. |
| `inner_polygons` | `semialg.boundary_topology` | function | Return polygon hole boundaries. |
| `polyhedron_vertices` | `semialg.boundary_topology` | function | Return unique vertices from polyhedral boundary topology. |
| `polyhedron_face_indices` | `semialg.boundary_topology` | function | Return face indices grouped by polyhedral shell. |
| `outer_polyhedra` | `semialg.boundary_topology` | function | Return outer polyhedral boundary shells. |
| `inner_polyhedra` | `semialg.boundary_topology` | function | Return cavity boundary shells. |
| `Ball` | `semialg.standard_regions` | class | Canonical closed Euclidean ball in arbitrary dimension. |
| `Sphere` | `semialg.standard_regions` | class | Canonical sphere boundary in arbitrary dimension. |
| `Circle` | `semialg.standard_regions` | constructor namespace | Construct planar canonical spheres. |
| `Ellipsoid` | `semialg.standard_regions` | class | Canonical filled ellipsoid with a positive-definite shape matrix. |
| `EllipsoidBoundary` | `semialg.standard_regions` | class | Canonical ellipsoid boundary with a positive-definite shape matrix. |
| `Cylinder` | `semialg.standard_regions` | class | Canonical flat-ended right circular cylinder. |
| `Cone` | `semialg.standard_regions` | class | Canonical right circular cone with base at `start` and apex at `end`. |
| `Torus` | `semialg.standard_regions` | class | Canonical ring or horn torus surface in three-dimensional space. |
| `FilledTorus` | `semialg.standard_regions` | class | Canonical solid ring or horn torus in three-dimensional space. |
| `BooleanRegion` | `semialg.standard_regions` | class | BooleanRegion(op: 'str', regions: 'Sequence[StandardRegion]', *, assume_disjoint: 'bool' = False) |
| `bounding_box` | `semialg.geometry_queries` | function | Compute the exact axis-aligned bounding box by coordinate optimization. |
| `Box` | `semialg.standard_regions` | class | Box(bounds: 'Sequence[tuple[object, object]]') |
| `Capsule` | `semialg.standard_regions` | class | Capsule(start: 'Sequence[object]', end: 'Sequence[object]', radius: 'object' = 1) |
| `centroid` | `semialg.moments` | function | Return the centroid of a finite-measure semialgebraic region. |
| `closest_points` | `semialg.derived_geometry` | function | Return all exact closest point pairs when the distance is attained. |
| `connected_components` | `semialg.derived_geometry` | function | Return exact CAD-connected-component formulas. |
| `convexity_certificate` | `semialg.convexity` | function | Decide convexity through the staged exact certificate hierarchy. |
| `contains_point` | `semialg.derived_geometry` | function | Return whether an exact point belongs to the semialgebraic region. |
| `RegionElement` | `semialg.symbolic_regions` | class | Symbolic assertion that a point belongs to a region. |
| `RegionNotElement` | `semialg.symbolic_regions` | class | Symbolic assertion that a point does not belong to a region. |
| `coordinate_range` | `semialg.derived_geometry` | function | Return the exact range of one coordinate over a region. |
| `covariance_matrix` | `semialg.moments` | function | Return the covariance matrix of the uniform measure on a region. |
| `critical_values` | `semialg.geometry_queries` | function | Return exact objective values from isolated and constant KKT components. |
| `critical_value_image` | `semialg.geometry_queries` | function | Return the exact supported value-set image of critical loci. |
| `CriticalValueImage` | `semialg.geometry_queries` | class | Structured exact image of isolated and positive-dimensional critical values. |
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
| `is_interior_disjoint` | `semialg.derived_geometry` | function | Return whether two regions have disjoint ambient interiors. |
| `Interval` | `semialg.standard_regions` | class | Interval(lower: 'object', upper: 'object', *, lower_closed: 'bool' = True, upper_closed: 'bool' = True) |
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
| `Parallelogram` | `semialg.standard_regions` | class | Parallelogram(origin: 'Sequence[object]', vectors: 'Sequence[Sequence[object]]') |
| `ParametricRegion` | `semialg.standard_regions` | class | ParametricRegion(parameters: 'Sequence[sp.Symbol \| str]', limits: 'Sequence[tuple[sp.Symbol \| str, object, object]]', mapping: 'Sequence[object]', *, multiplicity: 'object' = 1, assumptions: 'object' = True) |
| `path_between` | `semialg.geometry_queries` | function | Return a certified CAD cell-chain connecting two points in a region. |
| `FinitePointSet` | `semialg.standard_regions` | class | FinitePointSet(points: 'Sequence[Sequence[object]] \| Sequence[object]') |
| `TetrahedralComplex` | `semialg.standard_regions` | class | TetrahedralComplex(tetrahedra: 'Sequence[Simplex \| Sequence[Sequence[object]]]') |
| `region_boundary` | `semialg.regions.operations` | function | Return the Euclidean boundary of a semialgebraic region. |
| `region_closure` | `semialg.regions.operations` | function | Return the Euclidean closure of a semialgebraic region. |
| `region_complement` | `semialg.regions.operations` | function | Return the complement of an implicit semialgebraic region. |
| `region_difference` | `semialg.regions.operations` | function | Return ``lhs`` minus ``rhs`` for implicit semialgebraic regions. |
| `region_dimension` | `semialg.regions.operations` | function | Return the exact semialgebraic dimension from a complete adapted CAD. |
| `region_interior` | `semialg.regions.operations` | function | Return the Euclidean interior of a semialgebraic region. |
| `region_intersection` | `semialg.regions.operations` | function | Return the intersection of implicit semialgebraic regions. |
| `region_product` | `semialg.regions.operations` | function | Return the exact Cartesian product of semialgebraic regions. |
| `region_symmetric_difference` | `semialg.regions.operations` | function | Return the exact symmetric difference of two semialgebraic regions. |
| `region_union` | `semialg.regions.operations` | function | Return the union of implicit semialgebraic regions. |
| `scale` | `semialg.derived_geometry` | function | Scale a region about the origin by a scalar factor. |
| `semialgebraic_projection` | `semialg.geometry_queries` | function | Project ``region`` by existentially eliminating the requested variables. |
| `implicitize_polynomial_map` | `semialg.algebraic_geometry` | function | Implicitize a polynomial map by exact Groebner elimination. |
| `certified_radicalization` | `semialg.algebraic_geometry` | function | Return a certificate-verified reduced presentation of an affine variety. |
| `irreducible_components` | `semialg.algebraic_geometry` | function | Return a complete certified minimal-prime decomposition when available. |
| `reduced_component_singular_loci` | `semialg.algebraic_geometry` | function | Compute singular loci independently on certified reduced irreducible components. |
| `ReducedAlgebraicVariety` | `semialg.algebraic_geometry` | type | Certified reduced affine-variety presentation. |
| `IrreducibleAlgebraicComponent` | `semialg.algebraic_geometry` | type | Certified minimal-prime component. |
| `ReducedComponentSingularLocus` | `semialg.algebraic_geometry` | type | Component together with its reduced singular-locus formula. |
| `singular_locus` | `semialg.algebraic_geometry` | function | Return equations defining the singular locus of an algebraic variety. |
| `SphericalShell` | `semialg.standard_regions` | class | SphericalShell(center: 'Sequence[object]', radii: 'tuple[object, object]') |
| `squared_distance_range` | `semialg.derived_geometry` | function | Return the exact range of squared pairwise distances. |
| `Stadium` | `semialg.standard_regions` | class | Stadium(start: 'Sequence[object]', end: 'Sequence[object]', radius: 'object' = 1) |
| `StandardRegion` | `semialg.standard_regions` | class | Base class for explicit region objects supported by semialg. |
| `sublevel_set` | `semialg.derived_geometry` | function | Return ``region ∩ {expression <= value}`` (or strict variant). |
| `superlevel_set` | `semialg.derived_geometry` | function | Return ``region ∩ {expression >= value}`` (or strict variant). |
| `support_function` | `semialg.derived_geometry` | function | Return ``sup(x·direction)`` over the region. |
| `tangent_cone` | `semialg.algebraic_geometry` | function | Return the exact ideal-theoretic Zariski tangent cone at ``point``. |
| `tangent_dimension` | `semialg.algebraic_geometry` | function | Return the exact Zariski tangent-space dimension at ``point``. |
| `tangent_space` | `semialg.algebraic_geometry` | function | Return the Zariski tangent space at a point as the Jacobian nullspace. |
| `TransformedRegion` | `semialg.standard_regions` | class | TransformedRegion(base: 'StandardRegion', mapping: 'Sequence[object]', base_variables: 'Sequence[sp.Symbol \| str]') |
| `translate` | `semialg.derived_geometry` | function | Translate a region by ``vector`` while preserving coordinate symbols. |
| `width` | `semialg.derived_geometry` | function | Return exact directional width ``max u·x - min u·x``. |

## Integration, measure, and moments

Primary reference: [Integration, measure, and moments](integration_and_moments.md)

| Public name | Implementation module | Kind | Summary |
|---|---|---|---|
| `integrate_over_region` | `semialg.region_integrate` | function | Integrate ``integrand`` over a supported semialgebraic region. |
| `reduce_region_integral` | `semialg.region_integrate` | function | Reduce a supported region integral to explicit iterated integrals. |
| `region_moment` | `semialg.moments` | function | Return a raw moment integral over a semialgebraic region. |
| `region_measure` | `semialg.region_measurement` | function | Return intrinsic or ambient measure for canonical geometry and formula regions. |
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
| `SemialgebraicRegion` | `semialg.symbolic_regions` | class | A symbolic semialgebraic subset of ````R^n```` with lazy reusable state. |
| `as_semialgebraic_region` | `semialg.symbolic_regions` | function | Coerce a formula or explicit ````StandardRegion```` to ````SemialgebraicRegion````. |
| `local_dimension` | `semialg.region_analysis` | function | Exact local semialgebraic dimension at a point. |
| `region_active_boundary_strata` | `semialg.region_analysis` | function | Return exact pairwise-disjoint active-inequality boundary strata. |
| `region_nonsmooth_locus` | `semialg.region_analysis` | function | Return a conservative exact nonsmooth/corner locus of a region boundary. |
| `region_regular_locus` | `semialg.region_analysis` | function | Return the part of ````region```` outside its algebraic boundary singular locus. |
| `region_singular_locus` | `semialg.region_analysis` | function | Return an exact singular-locus formula or raise when certification is incomplete. |
| `region_singular_locus_result` | `semialg.region_analysis` | function | Return singular-locus geometry with completeness and diagnostic metadata. |
| `replay_certificate` | `semialg.certificates` | function | Replay a supported exact certificate using an independent public path. |

## Unified symbolic/CAD region facade additions

These exports provide the reusable symbolic-region and CAD-region API.

| Public name | Implementation module | Kind | Summary |
|---|---|---|---|
| `CADRegion` | `semialg.cad_region` | class | A reusable region represented by a public :class:`CADResult`. |
| `as_cad_region` | `semialg.cad_region` | function |  |
| `interior_of_closure` | `semialg.symbolic_regions` | function | Return interior(closure(region)). |
| `closure_of_interior` | `semialg.symbolic_regions` | function | Return closure(interior(region)). |
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
| `prove_nonzero` | `semialg.reasoning` | function | Prove that an expression is everywhere nonzero on a semialgebraic region. |
| `prove_zero` | `semialg.reasoning` | function | Prove that an expression is identically zero on a semialgebraic region. |
| `function_sign` | `semialg.reasoning` | function | Return a canonical exact sign classification, optionally stratified by parameters. |
| `function_sign_partition` | `semialg.function_analysis` | function | Partition a univariate semialgebraic function domain into exact positive, zero, and negative regions. |
| `function_smoothness` | `semialg.function_properties` | function | Report continuity, differentiability order, smoothness, and exact exceptional loci. |
| `function_mapping_properties` | `semialg.function_properties` | function | Certify injectivity, surjectivity, bijectivity, image, and exact witnesses. |
| `is_injective` | `semialg.function_properties` | function | Return the certified injectivity projection for a semialgebraic map. |
| `is_surjective` | `semialg.function_properties` | function | Return the certified surjectivity projection for a semialgebraic map. |
| `is_bijective` | `semialg.function_properties` | function | Return the certified bijectivity projection for a semialgebraic map. |
| `is_function_continuous` | `semialg.function_properties` | function | Return the certified continuity projection for a semialgebraic function. |
| `is_function_smooth` | `semialg.function_properties` | function | Return the certified smoothness projection for a semialgebraic function. |
| `strict_feasible` | `semialg.strict_feasibility` | function | Decide affine-relative or ambient strict feasibility and return a witness in result mode. |

## Algebraic decomposition and GTZ infrastructure

Primary reference: [Algebraic roots and exact finite solving](algebraic.md)

| Public name | Implementation module | Kind | Summary |
|---|---|---|---|

## Package lifecycle and cache controls

Primary reference: [Performance](../guides/performance.md#cache-lifecycle-and-teardown)

| Public name | Implementation module | Kind | Summary |
|---|---|---|---|

## Exact real-root and topology core

| Public name | Implementation module | Kind | Summary |
|---|---|---|---|
| `connected_component_count` | `semialg.roadmaps` | function | Count exact semialgebraic connected components. |
| `connected_component_samples` | `semialg.roadmaps` | function | Return one exact sample per connected component. |
| `topology_summary` | `semialg.topological_invariants` | function | Compute component, Euler, and supported Betti data. |
| `betti_number` | `semialg.topological_invariants` | function | Return a certified supported Betti number. |

## Semialgebraic topology and families

Primary reference: [Semialgebraic topology and families](regions.md#semialgebraic-topology-and-families)

| Public name | Implementation module | Kind | Summary |
|---|---|---|---|

| `minimal_prime_intersections` | `semialg.algebraic_geometry` | function | Return exact intersections of certified minimal-prime components. |
| `local_dimension_strata` | `semialg.algebraic_geometry` | function | Partition a reduced variety into constructible constant-local-dimension loci. |
| `stratified_singular_geometry` | `semialg.algebraic_geometry` | function | Separate component intersections, intrinsic singularities, and local-dimension strata. |
| `MinimalPrimeIntersection` | `semialg.algebraic_geometry` | type | Exact intersection of certified minimal-prime components. |
| `LocalDimensionStratum` | `semialg.algebraic_geometry` | type | Constructible constant-local-dimension locus. |
| `SingularGeometryStratum` | `semialg.algebraic_geometry` | type | Disjoint constructible singular-geometry stratum. |
| `StratifiedSingularGeometry` | `semialg.algebraic_geometry` | type | Certified stratified singular-geometry result. |

### Local branch geometry

| Public name | Implementation module | Kind | Summary |
|---|---|---|---|
| `local_branch_geometry` | `semialg.algebraic_geometry` | function | Compute certified tangent-cone, tangent-excess, transversality, and local multiplicity data at a point. |
| `polynomial_constraints` | `semialg.semialgebraic_local_geometry` | function | Extract a DNF-preserving polynomial constraint system. |
| `active_constraints` | `semialg.semialgebraic_local_geometry` | function | Identify exact active constraints at a feasible point. |
| `relative_interior` | `semialg.semialgebraic_local_geometry` | function | Compute interior relative to a certified affine hull. |
| `relative_boundary` | `semialg.semialgebraic_local_geometry` | function | Compute boundary relative to a certified affine hull. |
| `semialgebraic_tangent_cone` | `semialg.semialgebraic_local_geometry` | function | Compute the exact Bouligand tangent cone. |
| `parameterization_geometry` | `semialg.parameterization_geometry` | function | Analyze exact rank, image, fiber, and critical geometry of a parameterization. |
| `parameterization_critical_locus` | `semialg.parameterization_geometry` | function | Return the exact generic-rank-drop locus. |
| `parameterization_critical_values` | `semialg.parameterization_geometry` | function | Return the exact image of the critical locus. |
| `LocalBranchGeometry` | `semialg.algebraic_geometry` | type | Complete local branch-geometry result at a point. |
| `BranchTangentGeometry` | `semialg.algebraic_geometry` | type | Tangent cone, tangent space, excess, and multiplicity for one irreducible branch. |
| `ComponentIntersectionGeometry` | `semialg.algebraic_geometry` | type | Tangent and transversality data for the intersection of branches through a point. |

### Constraint certificates

| Public name | Implementation module | Kind | Summary |
|---|---|---|---|
| `implied_polynomial_inequality` | `semialg.constraint_certificates` | function | Certify a polynomial inequality by exact implication, nonnegative combination, or SOS certificate. |
| `redundant_polynomial_inequalities` | `semialg.constraint_certificates` | function | Identify polynomial inequalities implied by the remaining constraints. |
| `nonnegative_combination_certificate` | `semialg.constraint_certificates` | function | Construct an exact nonnegative-combination/SOS certificate. |
| `verify_nonnegative_combination_certificate` | `semialg.constraint_certificates` | function | Replay an exact nonnegative-combination/SOS certificate. |
| `component_constraint_descriptions` | `semialg.constraint_certificates` | function | Restrict a real model to each certified irreducible algebraic component. |

### Algebraic and constraint result types

| Public name | Implementation module | Kind | Summary |
|---|---|---|---|
| `PolynomialMapImplicitizationResult` | `semialg.algebraic_geometry` | type | Exact polynomial-map implicitization result and elimination metadata. |
| `ZariskiClosureResult` | `semialg.algebraic_geometry` | type | Certified Zariski-closure description. |
| `PolynomialConstraint` | `semialg.semialgebraic_local_geometry` | type | One normalized polynomial relation. |
| `PolynomialConstraintClause` | `semialg.semialgebraic_local_geometry` | type | One conjunctive clause of polynomial constraints. |
| `PolynomialConstraintSystem` | `semialg.semialgebraic_local_geometry` | type | Boolean-preserving polynomial constraint system. |
| `ActiveConstraintResult` | `semialg.semialgebraic_local_geometry` | type | Exact active-constraint classification at a point. |
| `ParameterizationGeometry` | `semialg.parameterization_geometry` | type | Exact parameterization rank, fiber, and critical-locus geometry. |
| `NonnegativeCombinationCertificate` | `semialg.constraint_certificates` | type | Replayable exact nonnegative-combination certificate. |
| `ImpliedPolynomialInequality` | `semialg.constraint_certificates` | type | Certified implied-inequality result. |
| `RedundantPolynomialInequality` | `semialg.constraint_certificates` | type | Certified redundant-inequality result. |
| `ComponentConstraintDescription` | `semialg.constraint_certificates` | type | Real constraint description attached to an irreducible component. |

### Polytope canonicalization and decomposition

Primary reference: [Regions](regions.md)

| Public name | Implementation module | Kind | Summary |
|---|---|---|---|
| `canonicalize_polygon` | `semialg.canonicalization` | function | Deterministically normalize supported polygonal geometry. |
| `canonicalize_polyhedron` | `semialg.canonicalization` | function | Canonicalize supported polyhedral boundary geometry without convexifying uncertified nonconvex shells. |
| `canonicalize_region` | `semialg.canonicalization` | function | Dispatch canonicalization for supported region families. |
| `PolytopeDecomposition` | `semialg.polytope_decomposition` | type | Certified exact simplicial decomposition of a convex polytope. |
| `triangulate_polytope` | `semialg.polytope_decomposition` | function | Exact deterministic pulling, placing, or barycentric triangulation. |
| `decompose_polytope` | `semialg.polytope_decomposition` | function | Produce an exact requested convex-polytope decomposition. |
| `MixedCellTetrahedralization` | `semialg.mixed_tetrahedralization` | type | Certified tetrahedralization of one convex 3-cell. |
| `MixedMeshTetrahedralization` | `semialg.mixed_tetrahedralization` | type | Conforming tetrahedralization of a mixed convex-cell mesh. |
| `tetrahedralize_cell` | `semialg.mixed_tetrahedralization` | function | Tetrahedralize a convex 3-cell using global vertex identifiers. |
| `tetrahedralize_cells` | `semialg.mixed_tetrahedralization` | function | Tetrahedralize adjacent cells with shared-face conformity. |

### Exact hulls, generation, refinement, and polyhedral Boolean operations

Primary reference: [Regions](regions.md#exact-hulls-generation-refinement-and-polyhedral-boolean-fast-paths)

| Public name | Implementation module | Kind | Summary |
|---|---|---|---|
| `convex_hull` | `semialg.convex_hull` | function | Construct an exact convex hull and optionally canonicalize it. |
| `random_polygon` | `semialg.generation_refinement` | function | Generate a seeded exact convex lattice polygon. |
| `random_polytope` | `semialg.generation_refinement` | function | Generate a seeded exact full-dimensional lattice polytope. |
| `subdivide_triangular_faces` | `semialg.generation_refinement` | function | Conformingly refine triangular faces with shared exact midpoints. |
| `geodesic_refinement` | `semialg.generation_refinement` | function | Refine triangular faces and project vertices exactly to a sphere. |
| `polyhedral_intersection` | `semialg.polyhedral_boolean` | function | Exact structural intersection of full-dimensional convex polytopes. |
| `polyhedral_boolean` | `semialg.polyhedral_boolean` | function | Structure-preserving convex-polyhedral Boolean fast path. |

### Canonical region and conversion API

Primary reference: [Regions](regions.md#canonicalization-polygonal-paths-conversion-and-polytope-decomposition)

| Public name | Implementation module | Kind | Summary |
|---|---|---|---|
| `Zonotope` | `semialg.standard_regions` | type | Minkowski sum of finitely many generator segments. |
| `RegionConversion` | `semialg.region_conversion` | type | Exactness/certification metadata for a region conversion. |
| `convert_region` | `semialg.region_conversion` | function | Convert among canonical, formula, CAD, parametric, simplicial, mesh, and boundary representations. |
| `polygonal_region_from_paths` | `semialg.polygonal_paths` | function | Build exact polygonal regions from crossing/winding paths. |
