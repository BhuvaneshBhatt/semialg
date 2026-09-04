# Regions reference
## Family contract

**Mathematical return.** Region APIs construct exact semialgebraic sets, decide exact set/geometric properties, or return exact images, projections, components, distances, transforms, and local algebraic-geometric objects.

**Exactness and certification.** Boolean/set relations are semantic rather than syntactic. Geometry queries reuse exact QE, CAD, optimization, and algebraic-geometry machinery; unsupported exact cases do not silently become numerical approximations.

**Algorithm.** Cheap substitutions/derived identities are used where possible; general projections/images and global geometric predicates may reduce to QE/CAD or exact optimization.

**Complexity and limitations.** Simple-looking geometric operations such as image, convexity, or connected components can be as hard as general QE. See [Limitations](../limitations.md).




## Primary API overview

This table is the substantive coverage target for the primary APIs assigned to this reference page. Each entry states the API's primary role; the family contract and detailed sections below explain shared algorithms, exactness guarantees, and limitations. It is maintained together with `docs/reference/primary_api_manifest.toml`, and documentation tests require every root-level primary API to map here rather than merely appearing in the generated public index.

| API | Kind | Role / return |
|---|---|---|
| `AffineMapAnalysis` | class | Exact rank, inverse, determinant, and scaled-isometry properties of an affine map. |
| `analyze_affine_map` | function | Analyze an affine expression map or matrix/offset pair without QE. |
| `AffineBoxClip` | class | Exact low-dimensional affine-subspace/box clipping result. |
| `clip_affine_subspace_to_box` | function | Clip a one- or two-dimensional affine subspace to a box by exact linear algebra. |
| `affine_transform` | function | Return the exact affine image ``A*x + b`` in the original coordinates. |
| `argmax_set` | function | Return the exact global maximizer set as a semialgebraic formula. |
| `argmin_set` | function | Return the exact global minimizer set as a semialgebraic formula. |
| `BallRegion` | class | BallRegion(center: 'Sequence[object]', radius: 'object' = 1) |
| `BooleanRegion` | class | BooleanRegion(op: 'str', regions: 'Sequence[StandardRegion]', *, assume_disjoint: 'bool' = False) |
| `bounding_box` | function | Compute the exact axis-aligned bounding box by coordinate optimization. |
| `BoxRegion` | class | BoxRegion(bounds: 'Sequence[tuple[object, object]]') |
| `CapsuleRegion` | class | CapsuleRegion(start: 'Sequence[object]', end: 'Sequence[object]', radius: 'object' = 1) |
| `centroid` | function | Return the exact centroid of a measurable semialgebraic region. |
| `closest_points` | function | Return all exact closest point pairs when the distance is attained. |
| `ConeRegion` | class | ConeRegion(start: 'Sequence[object]', end: 'Sequence[object]', radius: 'object' = 1) |
| `connected_components` | function | Return exact CAD-connected-component formulas. |
| `convexity_certificate` | function | Decide semialgebraic set convexity through an exact staged hierarchy. |
| `contains_point` | function | Return whether an exact point belongs to the semialgebraic region. |
| `coordinate_range` | function | Return the exact range of one coordinate over a region. |
| `covariance_matrix` | function | Return the exact covariance matrix of the uniform measure on a region. |
| `critical_values` | function | Return exact objective values from isolated and constant positive-dimensional KKT components. |
| `CylinderRegion` | class | CylinderRegion(start: 'Sequence[object]', end: 'Sequence[object]', radius: 'object' = 1) |
| `diameter` | function | Return the exact Euclidean diameter (supremal pairwise distance). |
| `distance_between_regions` | function | Compute exact Euclidean distance between two semialgebraic regions. |
| `distance_set` | function | Return the exact set of Euclidean pairwise distances as a formula. |
| `distance_to_region` | function | Compute exact Euclidean distance from a point to a semialgebraic region. |
| `euler_characteristic` | function | Compute the semialgebraic Euler characteristic from selected CAD cells. |
| `extrema_set` | function | Return the union of the exact global minimum and maximum sets. |
| `fiber` | function | Specialize a semialgebraic family at fixed parameter/coordinate values. |
| `has_empty_interior` | function | Return whether the region has empty ambient interior. |
| `inertia_tensor` | function | Return the unit-density second moment-of-inertia tensor about the origin. |
| `intersects` | function | Return whether two semialgebraic regions have nonempty intersection. |
| `IntervalRegion` | class | IntervalRegion(lower: 'object', upper: 'object', *, lower_closed: 'bool' = True, upper_closed: 'bool' = True) |
| `is_bounded` | function | Return whether the semialgebraic region is bounded. |
| `is_closed` | function | Return whether the semialgebraic region is closed. |
| `is_compact` | function | Return whether the semialgebraic region is compact. |
| `is_connected` | function | Decide connectedness; for semialgebraic sets this equals path connectedness. |
| `is_convex` | function | Return whether a semialgebraic region is convex, exactly. |
| `is_dense_in` | function | Return whether ``subset`` is dense in ``ambient`` in the ambient Euclidean topology. |
| `is_disjoint` | function | Return whether two semialgebraic regions are disjoint. |
| `is_empty` | function | Return whether the semialgebraic region is empty. |
| `is_equal` | function | Return whether two semialgebraic regions define the same set. |
| `is_full_dimensional` | function | Return whether the region has full dimension in its ambient variables. |
| `is_open` | function | Return whether the semialgebraic region is open. |
| `is_path_connected` | function | Decide path connectedness via exact CAD connectivity. |
| `is_singular` | function | Return whether ``point`` is singular on the polynomial variety. |
| `is_smooth` | function | Return whether the real polynomial variety has empty singular locus. |
| `is_subset` | function | Return whether one semialgebraic region is contained in another. |
| `level_set` | function | Return ``region ∩ {expression = value}``. |
| `linear_image` | function | Return the exact linear image ``A*x`` in the original coordinates. |
| `minkowski_sum` | function | Return the exact Minkowski sum of two semialgebraic regions. |
| `moment_matrix` | function | Return the normalized raw second-moment matrix ``E[x x.T]``. |
| `nearest_point` | function | Return all exact nearest points when the distance is attained. |
| `ParallelepipedRegion` | class | ParallelepipedRegion(origin: 'Sequence[object]', vectors: 'Sequence[Sequence[object]]') |
| `ParallelogramRegion` | class | ParallelogramRegion(origin: 'Sequence[object]', vectors: 'Sequence[Sequence[object]]') |
| `ParametricRegion` | class | ParametricRegion(parameters: 'Sequence[sp.Symbol \| str]', limits: 'Sequence[tuple[sp.Symbol \| str, object, object]]', mapping: 'Sequence[object]', *, multiplicity: 'object' = 1, assumptions: 'object' = True) |
| `ParametricChart` | class | One certified bounded parameter-domain map into a semialgebraic region. |
| `ParametricCover` | class | Finite certified collection of parametric charts covering a region or bounded intersection. |
| `bounded_parametric_cover` | function | Build a structural bounded cover without CAD when recognized geometry is available. |
| `ParametricMapDegree` | class | Generic algebraic fiber-degree result for a rational parametrization. |
| `parametric_map_degree` | function | Compute generic complex fiber multiplicity from a quotient-algebra degree. |
| `path_between` | function | Return a certified CAD cell-chain connecting two points in a region. |
| `PointRegion` | class | PointRegion(points: 'Sequence[Sequence[object]] \| Sequence[object]') |
| `PolygonRegion` | class | PolygonRegion(vertices: 'Sequence[Sequence[object]]') |
| `PolyhedronRegion` | class | PolyhedronRegion(tetrahedra: 'Sequence[TetrahedronRegion \| Sequence[Sequence[object]]]') |
| `PrismRegion` | class | PrismRegion(base: 'PolygonRegion \| SimplexRegion \| Sequence[Sequence[object]]', vector: 'Sequence[object]') |
| `PyramidRegion` | class | PyramidRegion(base: 'PolygonRegion \| SimplexRegion \| Sequence[Sequence[object]]', apex: 'Sequence[object]') |
| `region_boundary` | function | Return the Euclidean boundary of a semialgebraic region. |
| `region_closure` | function | Return the Euclidean closure of a semialgebraic region. |
| `region_complement` | function | Return the complement of an implicit or unified semialgebraic region. |
| `region_components` | function | Return connected-component formulas for simple explicit cases. |
| `region_difference` | function | Return ``lhs`` minus ``rhs`` for implicit or unified regions. |
| `region_dimension` | function | Return the exact semialgebraic dimension from a complete adapted CAD. |
| `region_interior` | function | Return the Euclidean interior of a semialgebraic region. |
| `region_intersection` | function | Return the intersection of implicit or unified semialgebraic regions. |
| `region_union` | function | Return the union of implicit or unified semialgebraic regions. |
| `RegionDifference` | function | Return the Boolean difference of two standard regions. |
| `RegionIntersection` | function | Return the Boolean intersection of standard regions. |
| `RegionSymmetricDifference` | function | Return the Boolean symmetric difference of two standard regions. |
| `RegionUnion` | function | Return the Boolean union of standard regions. |
| `scale` | function | Scale a region about the origin by a scalar factor. |
| `semialgebraic_image` | function | Return the exact semialgebraic image of a polynomial/rational map. |
| `semialgebraic_preimage` | function | Return the preimage of a semialgebraic target under a symbolic map. |
| `semialgebraic_projection` | function | Project ``region`` by existentially eliminating the requested variables. |
| `SimplexRegion` | class | SimplexRegion(vertices: 'Sequence[Sequence[object]]') |
| `singular_locus` | function | Return equations defining the singular locus of an algebraic variety. |
| `SphereRegion` | class | SphereRegion(center: 'Sequence[object]', radius: 'object' = 1) |
| `SphericalShellRegion` | class | SphericalShellRegion(center: 'Sequence[object]', radii: 'tuple[object, object]') |
| `squared_distance_range` | function | Return the exact range of squared pairwise distances. |
| `StadiumRegion` | class | StadiumRegion(start: 'Sequence[object]', end: 'Sequence[object]', radius: 'object' = 1) |
| `StandardRegion` | class | Base class for explicit region objects supported by semialg. |
| `sublevel_set` | function | Return ``region ∩ {expression <= value}`` (or strict variant). |
| `superlevel_set` | function | Return ``region ∩ {expression >= value}`` (or strict variant). |
| `support_function` | function | Return ``sup(x·direction)`` over the region. |
| `tangent_cone` | function | Return the exact ideal-theoretic Zariski tangent cone at ``point``. |
| `tangent_dimension` | function | Return the exact Zariski tangent-space dimension at ``point``. |
| `tangent_space` | function | Return the Zariski tangent space at a point as the Jacobian nullspace. |
| `TetrahedronRegion` | class | TetrahedronRegion(vertices: 'Sequence[Sequence[object]]') |
| `TransformedRegion` | class | TransformedRegion(base: 'StandardRegion', mapping: 'Sequence[object]', base_variables: 'Sequence[sp.Symbol \| str]') |
| `translate` | function | Translate a region by ``vector`` while preserving coordinate symbols. |
| `width` | function | Return exact directional width ``max u·x - min u·x``. |
| `SemialgebraicContext` | class | Reusable normalized semialgebraic problem plus exact computation cache. |
| `SemialgebraicRegion` | class | A symbolic semialgebraic subset of ``R^n`` with lazy reusable state. |
| `as_semialgebraic_region` | function | Coerce a formula or explicit ``StandardRegion`` to ``SemialgebraicRegion``. |
| `local_dimension` | function | Exact local semialgebraic dimension at a point. |
| `region_active_boundary_strata` | function | Return pairwise-disjoint exact strata classified by active inequality boundaries. |
| `BoundaryStratum` | class | One exact boundary CAD cell with inclusion and active-constraint metadata. |
| `RegionBoundaryResult` | class | Exact boundary formula plus reusable CAD and boundary strata. |
| `SingularLocusResult` | class | Certified singular-locus result with explicit incomplete-result metadata. |
| `region_boundary_result` | function | Compute a rich exact boundary description while retaining decomposition metadata. |
| `region_nonsmooth_locus` | function | Return algebraic singularities and transverse inequality-boundary corners/ridges. |
| `region_regular_locus` | function | Return the part of ``region`` outside its algebraic boundary singular locus. |
| `region_singular_locus` | function | Return an exact singular-locus formula or raise when certification is incomplete. |
| `region_singular_locus_result` | function | Return singular-locus geometry with completeness and diagnostic metadata. |
| `replay_certificate` | function | Replay a supported exact certificate using an independent public path. |
| `result_diagnostics` | function | Return one stable diagnostic schema for CAD/QE/optimization certificates. |
| `CADRegion` | class | A reusable region represented by a public :class:`CADResult`. |
| `as_cad_region` | function | Coerce a region/formula/CAD result to a reusable :class:`CADRegion`. |
| `region_closure_interior` | function | Return interior(closure(region)). |
| `region_interior_closure` | function | Return closure(interior(region)). |
| `region_variables` | function | Return coordinate variables, parameters, or all symbols of a region. |
| `is_regular_closed_region` | function | Return whether a region equals the closure of its interior. |
| `is_regular_open_region` | function | Return whether a region equals the interior of its closure. |
| `simplify_region` | function | Canonicalize a symbolic semialgebraic region formula. |


## Parametric covers and structural dimension

`bounded_parametric_cover` exposes bounded geometry as a finite collection of
`ParametricChart` objects. Recognized boxes, simplices, polygons, polyhedra,
parallelepipeds, points, intervals, and bounded `ParametricRegion` objects use
their structural parametrizations. Formula regions can be clipped by explicit
finite bounds and represented by an exact identity chart without invoking CAD.

A chart carries explicit parameter bounds, an additional parameter condition,
and the coordinate mapping. When the parameter-domain dimension and generic
Jacobian rank are certified, `region_dimension` reuses that information before
falling back to a complete adapted CAD. The fallback preserves the existing
exact dimension contract for cases whose chart rank or domain dimension cannot
be certified structurally.

## Boundary metadata

`region_boundary_result` computes one exact boundary CAD and retains it on the
result. Each `BoundaryStratum` records its exact cell formula, Euclidean
dimension, whether that boundary cell belongs to the original region, and the
recognized active inequality residuals. This explicitly distinguishes, for
example, the included endpoint of `[0, 1)` from its excluded endpoint and lets
plotting, meshing, topology, and singularity routines reuse the same boundary
decomposition.

## Distance critical points and affine geometry

For a smooth pure equality variety, `distance_to_region` first constructs the
exact Lagrange critical-point system for squared distance. When that system is
zero-dimensional, the algebraic solver returns all real critical points and the
nearest one is selected by exact comparison; otherwise the general exact
optimization backend remains the fallback.

`analyze_affine_map` centralizes affine rank, determinant, inverse, and scaled
isometry recognition. `affine_transform` reuses this analysis for its exact
invertible-map shortcut. `parametric_map_degree` computes the generic algebraic
fiber degree of polynomial/rational maps; real-domain restrictions are kept
separate from this generic complex degree so later integration code can combine
the two deliberately.

`clip_affine_subspace_to_box` is intentionally narrow: it handles one- and
two-dimensional affine subspaces by exact facet intersection and does not
replace general region intersection. It exists for geometry/meshing paths where
the specialized linear-algebra computation is measurably cheaper than CAD.

## Formula-based region operations

- `region_union`
- `region_intersection`
- `region_difference`
- `region_complement`
- `region_closure`
- `region_interior`
- `region_boundary`

Structural and predicate APIs include `region_dimension`, `region_components`, `region_subset`, `region_equal`, `region_disjoint`, `region_bounded`, `region_closed`, and `region_compact`.

## Standard regions

The standard-region hierarchy represents common geometry directly:

- points and intervals;
- boxes;
- simplices and tetrahedra;
- polygons/polyhedra;
- parallelograms/parallelepipeds;
- prisms and pyramids;
- balls, spheres, and spherical shells;
- cylinders and cones;
- stadiums and capsules;
- parametric and transformed regions;
- Boolean combinations of standard regions.

Constructor invariants are part of the public contract. See [Region invariants](../guides/region_invariants.md).

## Boolean standard regions

`RegionUnion`, `RegionIntersection`, `RegionDifference`, and `RegionSymmetricDifference` construct set-theoretic combinations of `StandardRegion` objects.


## Unified symbolic and reusable CAD regions

`SemialgebraicRegion(formula, variables)` is the canonical symbolic region value. It exposes membership, Boolean operations, exact projection/image/preimage, topology, local dimension, algebraic singular/regular loci, simplification, measure, and integration while lazily caching `SemialgebraicContext` and an optional CAD.

Use `region.as_cad_region()` when repeated decomposition-level queries justify constructing a complete CAD. `CADRegion` exposes `locate_point`, exact projection-tower `sign_vector`, `cell_complex`, Euler characteristic, CAD signatures, CAD Boolean combination/extension, and direct CAD-cell integration.

`CADCellComplex` groups selected cells by dimension and records exact codimension-one closure incidence. Incidence prefers the original CAD closure relation, which avoids reconstructing radical boundary expressions and relaunching QE.

See [Region representations](../guides/region_representations.md), [Reuse one CAD](../guides/cad_reuse.md), and [Exact versus numerical region results](../guides/exact_vs_numerical_regions.md).

## Numerical CAD geometry

`StructuredCADCell` exposes nested typed bounds and algebraic boundary descriptors. `triangulate_cad_cell`, `triangulate_cad_cells`, and `triangulate_cad_region` produce numerical simplicial meshes with source-cell provenance. Multi-cell meshing can require sampled face conformity and shares coincident boundary vertices. Numerical curve/surface evaluators preserve the CAD root index used for branch identity.

These meshes and sampled coordinates are not exact topology certificates. See the [region capability matrix](../quality/region_capability_matrix.md).

## Exact endpoint semantics

One-dimensional component merging and region intersections compare exact endpoints symbolically. Close algebraic values are not ordered by fixed-precision decimal conversion.

See [Region operations](../region_operations.md).


## High-level semialgebraic maps and metric queries

`semialgebraic_projection(region, eliminate, variables)` computes a Tarski-Seidenberg projection by existential QE. `semialgebraic_image(mapping, domain, ...)` constructs the graph of a symbolic map and projects away the source variables; `semialgebraic_preimage` substitutes a map into a target formula; `fiber` specializes a family at fixed parameter or coordinate values.

`bounding_box` optimizes every coordinate exactly. `distance_to_region` and `distance_between_regions` minimize squared Euclidean distance through the exact optimization backend and optionally return the underlying optimization certificate/result. `critical_values` exposes exact values at isolated KKT/singular candidates and additionally certifies values on positive-dimensional connected KKT components when the objective is constant there; attained global extrema are also included.

## Convexity and path connectivity

`is_convex` first recognizes affine/polyhedral intersections exactly. For basic polynomial sets it certifies convex sublevel and concave superlevel constraints from Hessian semidefiniteness: constant Hessians are discharged directly, while variable polynomial Hessians are certified exactly through the principal-minor characterization of positive semidefiniteness, using complete CAD only for unresolved sign conditions in the original variables. Nonlinear equalities are not accepted by this fast path. Any unrecognized case falls back to the defining first-order convexity sentence and complete QE.

`is_path_connected` uses CAD connectivity. For semialgebraic sets, connectedness implies semialgebraic path connectedness, so this gives an exact decision when the CAD connectivity extraction completes. `path_between` returns an explicit piecewise-linear path for the current one-dimensional certified case. In higher dimensions it returns a certified chain of CAD cells and connector witnesses; it deliberately does not claim to construct a full algebraic roadmap parameterization.

## Euler characteristic

`euler_characteristic(region, variables, compact_support=True)` computes the additive compactly-supported semialgebraic Euler characteristic from selected CAD cells as `sum((-1)**cell.dimension)`. On compact regions this agrees with the ordinary Euler characteristic. For example a closed interval has value `1`, while an open interval has compactly-supported value `-1`.

## Local algebraic geometry

`singular_locus` applies the Jacobian criterion using an inferred or explicitly supplied codimension. `tangent_space` returns the evaluated Jacobian and its nullspace. `tangent_cone` computes the exact ideal-theoretic tangent-cone ideal `in_m(I)` at the requested point. It uses the m-adic deformation `t**(-ord(f))*f(t*d)`, saturation by `t`, and specialization at `t=0`; this captures initial forms created by cancellations between generators and is therefore independent of the supplied generating set. The `initial_forms` field contains a reduced Gröbner basis of that exact ideal, also available as `ideal_generators`.


## Derived exact geometry API

High-level derived operations are available directly from `semialg`:

- extrema loci: `argmin_set`, `argmax_set`, `extrema_set`;
- exact set predicates: `is_empty`, `is_bounded`, `is_compact`, `is_open`, `is_closed`, `is_subset`, `is_equal`, `is_disjoint`, `intersects`, `is_dense_in`, and `contains_point`;
- topology/structure: `connected_components`, `is_connected`, `is_full_dimensional`, and `has_empty_interior`;
- metric and convex geometry: `coordinate_range`, `diameter`, `nearest_point`, `closest_points`, `support_function`, `width`, `squared_distance_range`, and `distance_set`;
- transforms: `translate`, `scale`, `linear_image`, `affine_transform`, and `minkowski_sum`;
- level constructions: `level_set`, `sublevel_set`, and `superlevel_set`.

These are compositions of the certified lower-level algorithms. Direct substitutions such as translation are cheap; general images, Minkowski sums, and some metric operations may invoke optimization or complete CAD/QE and retain its worst-case complexity.

## Region wrappers, regularization, and display geometry

`BooleanRegion` represents named Boolean combinations of standard regions, and
`TransformedRegion` represents an exact transformed standard region before it is
lowered to a semialgebraic formula.

`region_closure_interior` computes the closure of a region's interior and
`region_interior_closure` computes the interior of its closure. The predicates
`is_regular_closed_region` and `is_regular_open_region` test the corresponding
fixed-point identities. `simplify_region` removes Boolean redundancy while
preserving exact region semantics.

`region_singular_locus` uses reduced-real component-relative algebraic regularity
and analyzes an inequality boundary only on exact CAD cells where that residual is
actually active. Boundary sections are tested relative to the corresponding equality
component, and the result is restricted to the actual topological boundary. The
formula-only API raises `NotImplementedError` if those component-sensitive proof
obligations cannot be certified; it never falls back to one global maximum-dimension
Jacobian threshold. `region_singular_locus_result` exposes the same computation with
explicit `complete`, `formula`, `known_singular_formula`, and diagnostic metadata.
`region_nonsmooth_locus` additionally detects transverse intersections on realized
multi-active inequality strata (polygonal corners and higher-dimensional ridges).
`region_active_boundary_strata` reuses the same rich boundary metadata and returns
exact pairwise-disjoint formulas classified by the realized active residual set; it
is a defining-constraint stratification and does not claim Whitney regularity.
`region_regular_locus` remains the complement of the algebraic singular locus;
`local_dimension` complements them by measuring the largest incident CAD-cell
dimension at a point.

`discretize_solution` and `discretize_region_geometry` produce numerical display
data from exact solution or region objects. `plot_solution` and
`plot_region_geometry` are convenience plotting functions built on those
representations; plotting does not change the exact mathematical object.
