# semialg

`semialg` is a Python package for exact computational semialgebraic  and real algebraic geometry. It combines cylindrical algebraic decomposition (CAD), quantifier elimination (QE), exact real-algebraic arithmetic, polynomial-system and algebraic-decomposition methods, and specialized algorithms for solving, decision problems, geometric analysis, optimization, and integration.

The package has certified paths that prefer an exact answer (or an explicit failure).

## Why use semialg?

Use semialg when the distinction between a plausible numerical answer and a mathematically certified answer matters. It is especially useful when a problem mixes polynomial equations and inequalities, Boolean conditions, parameters, (exact) algebraic numbers, or global questions such as feasibility, equivalence, projection, optimization, topology, or measure. Specialized algebraic and geometric methods handle inexpensive cases first; CAD and QE provide fallback machinery for supported real-polynomial formulas.

A result being symbolic is not by itself a certificate. semialg distinguishes between  exact representations, certified conclusions, candidate/heuristic information, and explicitly numerical approximations. See [Exactness and certification](https://github.com/BhuvaneshBhatt/semialg/blob/main/docs/concepts/exactness_and_certification.md).

If you primarily need floating-point nonlinear optimization, general numerical transcendental solving, or large numerical geometry (instead of exact semialgebraic computation), established numerical libraries will usually be a better fit. Consider SciPy for numerical optimization and nonlinear systems, SymPy/mpmath for numerical transcendental root finding, Shapely for planar computational geometry, and Trimesh for numerical 3-D mesh geometry. See the [capability matrix](https://github.com/BhuvaneshBhatt/semialg/blob/main/docs/feature_matrix.md) and [limitations](https://github.com/BhuvaneshBhatt/semialg/blob/main/docs/limitations.md) before choosing an algorithm.

## Install

```bash
python -m pip install semialg
```

For development:

```bash
python -m pip install -e .[dev]
```

## Quick start

```python
import sympy as sp
from semialg import equivalent, implies, is_satisfiable

x, y = sp.symbols("x y", real=True)

is_satisfiable((x**2 + y**2 <= 1) & (x > 0) & (y > 0), [x, y])
# True

implies(x > 1, x**2 > 1, [x])
# True

equivalent(x**2 <= 1, (x >= -1) & (x <= 1), [x])
# True
```

For a decision query, semialg first tries inexpensive structure before falling back to general elimination machinery:

![Decision computation flow](https://raw.githubusercontent.com/BhuvaneshBhatt/semialg/main/docs/assets/decision-flow.svg)

The [computation-flow guide](https://github.com/BhuvaneshBhatt/semialg/blob/main/docs/concepts/computation_flows.md) shows the corresponding CAD, optimization, and integration pipelines.

First-order formulas can be built directly with semialg's symbolic quantifiers:

```python
from semialg import Exists, ForAll
from semialg.solve import reduce_complete_expr

statement = ForAll(x, Exists(y, sp.Eq(x + y, 0)))
reduce_complete_expr(statement)
# True
```

Optimization and integration:

```python
from semialg import semialgebraic_minimize, semialgebraic_measure

opt = semialgebraic_minimize(x**2 + y**2, x + y >= 1, [x, y], return_result=True)
opt.value
# 1/2

semialgebraic_measure(x**2 + y**2 <= 1, [x, y])
# pi
```

## Core workflows

### Regions and reusable CAD decompositions

Composite Boolean regions are interpreted geometrically through CAD, not atom-by-atom rewriting. The same disjoint CAD-cell decomposition is used by multidimensional measure and integration, so overlapping unions are not double-counted and internal Boolean seams are not reported as boundaries.

The canonical symbolic region API is `SemialgebraicRegion`; build a `CADRegion` only when repeated point-location, topology, cell-complex, or CAD-integration queries justify reusing one decomposition:

```python
from semialg import SemialgebraicRegion

disk_region = SemialgebraicRegion(x**2 + y**2 <= 1, (x, y))
cad_region = disk_region.as_cad_region()

location = cad_region.locate_point((0, 0))
complex_ = cad_region.cell_complex()
location.selected
# True
complex_.euler_characteristic()
# 1
```

See the [region representation guide](https://github.com/BhuvaneshBhatt/semialg/blob/main/docs/guides/region_representations.md), [CAD reuse guide](https://github.com/BhuvaneshBhatt/semialg/blob/main/docs/guides/cad_reuse.md), and [exact-versus-numerical guide](https://github.com/BhuvaneshBhatt/semialg/blob/main/docs/guides/exact_vs_numerical_regions.md). The [API surface policy](https://github.com/BhuvaneshBhatt/semialg/blob/main/docs/architecture/api_surface.md) distinguishes primary root-level names from expert-level/specialized submodule APIs and internal implementation details.

### Parameter-dependent computation

Parameter-dependent region integrals return certified guarded answers:

```python
a = sp.Symbol("a", real=True)
from semialg import integrate_over_region

parametric = integrate_over_region(
    x,
    (x >= 0) & (x <= a),
    [x],
    parameters=[a],
    return_stratified=True,
)
parametric.select({a: 2})
# 2
parametric.select({a: -1})
# 0
```

Parametric optimization and range APIs keep quantified relations by default. Pass `eliminate_quantifiers=True` together with `return_stratified=True` when a quantifier-free relation is worth the additional CAD/QE cost. Common affine one-dimensional parameter families are reconstructed directly before a second CAD/QE is attempted.

### Quantifier elimination and CAD planning

The complete-QE path also has a structural presolver and automatic CAD ordering. Safe affine equalities and linear innermost existential blocks are eliminated before CAD; Brown-style ordering is applied only within semantically interchangeable free/quantifier blocks. Projection-set scoring is available explicitly through `suggest_cad_variable_order(..., strategy="projection")`.

```python
from semialg.heuristics import suggest_cad_variable_order
from semialg.presolve import presolve_semialgebraic

presolved = presolve_semialgebraic(sp.And(sp.Eq(y, x + 1), y > 0), [x, y], eliminate=[y])
presolved.formula
# x > -1

score = suggest_cad_variable_order([x**2 + y**4, x * y + 1], [x, y])
set(score.order)
# {x, y}
```

Univariate parameter-dependent integrals can use algebraic CAD root functions as moving endpoints. For example, the measure of `x**2 <= a` is stratfied exactly as zero for `a <= 0` and the distance between the two ordered roots for `a > 0`; specialziation gives `4` at `a=4` and `6` at `a=9`.

### Geometric queries

High-level geometric queries compose the same certified backends:

```python
from semialg import (
    argmin_set,
    bounding_box,
    centroid,
    connected_components,
    diameter,
    distance_to_region,
    is_compact,
    minkowski_sum,
    region_image,
    singular_locus,
    tangent_cone,
    width,
)

u = sp.Symbol("u", real=True)
region_image((x >= -1) & (x <= 2), x, [x], image_variables=[u])
# (u >= -1) & (u <= 2)

bounding_box((x >= -1) & (x <= 2), [x])
# {x: (-1, 2)}

distance_to_region((2,), (x >= 0) & (x <= 1), [x])
# 1

singular_locus(y**2 - x**3, [x, y]).subs({x: 0, y: 0})
tangent_cone(y**2 - x**3, {x: 0, y: 0}, [x, y]).certified
# Exact ideal-theoretic cone for multiple generators; cancellation terms are retained.
tangent_cone((x**2 + y**3, x**2 - y**3), {x: 0, y: 0}, [x, y]).ideal_generators
# Eq(d_y**2, 0)
```

Derived geometry APIs expose common constructions without duplicating algorithms:

```python
interval = (x >= 0) & (x <= 2)

argmin_set(x**2, interval, [x])
# (x >= 0) & (x <= 2) & Eq(x**2, 0)

diameter(interval, [x])
# 2

width(interval, (1,), [x])
# 2

is_compact(interval, [x])
# True

connected_components((x < 0) | (x > 0), [x])
# (x < 0, x > 0)

centroid(interval, [x])
# {x: 1}
```

The same layer includes set relations (`is_subset`, `is_equal`, `is_disjoint`, `is_interior_disjoint`, `intersects`), level/sublevel/superlevel sets, nearest/closest points, moment/covariance/inertia matrices, affine transforms and Minkowski sums, support functions, directional widths, and algebraic `is_smooth` / `is_singular` / `tangent_dimension` conveniences. Operations such as Minkowski sums and general affine images remain exact existential-QE problems and can therefore inherit CAD complexity.

Region topology uses CAD semantics even for atomic polynomial inequalities; compact direct formulas are used only after cell-wise equivalence verification. `region_dimension` is likewise exact and is read from selected CAD-cell dimensions. Free parameters are preserved by `region_image` when source variables are explicit, and all geometry point/string inputs use contextual symbol resolution.

`is_path_connected` is an exact decision through CAD connectivity. `path_between` returns an explicit piecewise-linear path in the certified 1-D case and a certified CAD cell/connector chain in higher dimensions; full Canny/Basu-Pollack-Roy roadmap parameterization is not implemented yet.

## Primary and specialized APIs

The package root is the everyday mathematical interface. Specialized algorithms and certificate-building functions live in their owning namespaces instead of being duplicated at `semialg.*`. For example, use `semialg.parameters.root_count_conditions`, `semialg.roadmaps.roadmap`, `semialg.topology.semialgebraic.triangulate_region`, `semialg.parametric_geometry.bounded_parametric_cover`, and `semialg.map_degree.parametric_map_degree` when those are needed. See [API namespaces](https://github.com/BhuvaneshBhatt/semialg/blob/main/docs/api_namespaces.md).

Semialgebraic function analysis is composition-aware for the supported graph fragment. Nested radicals/rational powers, `Abs`, `sign`, `Heaviside`, `Min`/`Max`, finite `Piecewise`, and real `re`/`im`/`conjugate` wrappers can be algebraized when their real-domain side conditions are certifiable; unsupported graphs are declined.

## What semialg provides

- **Decision and QE:** satisfiability, tautology, implication, equivalence, CAD and quantifier elimination.
- **Specialized real algebraic decisions:** certified positive-dimensional equality feasibility, replayable algebraic decomposition, polynomial nonnegativity, and one-sided exact witness search.
- **Solving and sampling:** an exact solver with affine presolve, Boolean/incidence decomposition, univariate and linear fast paths, zero-dimensional RUR, Groebner presolve, and CAD/QE fallback, plus structured witnesses and samples.
- **Algebraic roots and parameters:** exact root isolation, certified algebraic root functions, root classification, and parameter-stratified results.
- **Toric and lattice algebra:** saturated integer kernels, lattice and toric ideals, Laurent-monomial elimination, and exact Markov bases.
- **Algebraic statistics:** exact semialgebraic local geometry for polynomial statistical models, active/redundant constraints, tangent cones and local dimension, parameterization geometry, toric/lattice constructions, and certificate-aware hypothesis-testing geometry.
- **Regions:** Boolean region operations, exact set relations/properties, topology, CAD connected components, transforms, standard geometric regions, and CAD-derived geometry.
- **Geometric queries:** exact projection/image/preimage/fibers, extrema sets, bounding boxes, distances/diameter, support/width, transforms/Minkowski sums, staged exact convexity certificates, connectivity/property predicates, Euler characteristic, moments, singular loci, tangent spaces, and tangent cones.
- **Algebraic image geometry:** exact polynomial-map implicitization, Zariski closures by Groebner elimination, and Jacobian-criterion singular loci.
- **Optimization and ranges:** exact polynomial optimization, KKT/active-set analysis, global certification, and semialgebraic image/range computation.
- **Integration and moments:** ambient and intrinsic measure, region integrals, moments, centroids, and covariance.

See the [feature matrix](https://github.com/BhuvaneshBhatt/semialg/blob/main/docs/feature_matrix.md) for a more detailed capability summary and [limitations](https://github.com/BhuvaneshBhatt/semialg/blob/main/docs/limitations.md) for important scope boundaries.

## Applied workflows

`semialg.applications` provides thin workflows built on the certified core:

- **Robust parameter analysis** — derive exact parameter regions for existential feasibility or universal satisfaction.
- **Symbolic-math validation** — certify identities, formula equivalence, and proposed function ranges.
- **Optimization benchmark oracle** — produce certified exact optima for checking numerical optimizers.
- **Polynomial control stability** — derive exact strict Hurwitz-stability regions from characteristic polynomials.
- **Polynomial safety invariants** — certify initiation, inductiveness, and unsafe-state exclusion for discrete polynomial systems.
- **Polynomial response surfaces** — compute exact extrema, ranges, gradients, and threshold regions for polynomial surrogate models.
- **Polynomial model comparison** — certify worst-case discrepancy, dominance, and equivalence over a domain.
- **Parameter regime analysis** — partition parameter space by solvability or real-root count.
- **Polynomial probability** — integrate exact polynomial densities over semialgebraic events and supports.

See [`docs/applications.md`](https://github.com/BhuvaneshBhatt/semialg/blob/main/docs/applications.md). Core operations such as `function_range`, `semialgebraic_measure`, integration, optimization, CAD, and QE remain in the core namespace; they are not duplicated under `applications`.

## Selected advanced capabilities

The sections below highlight capabilities that go beyond the introductory decision, region, optimization, and integration workflows. For the complete inventory, see the [feature matrix](https://github.com/BhuvaneshBhatt/semialg/blob/main/docs/feature_matrix.md).

### Region sampling, measure, and centroid

Canonical geometry objects expose `sample_point()` / `sample_points()` for certified representative witnesses and `random_point()` / `random_points()` for seeded distributional sampling. Supported bounded canonical regions are sampled uniformly with respect to intrinsic Euclidean/Hausdorff measure. `Geometry.measure()` defaults to intrinsic measure, while `measure(measure_dimension="ambient")` returns ambient Lebesgue measure. The root `region_measure(..., measure_dimension=None)` follows the same intrinsic default for canonical geometry but uses ambient measure for formula regions; pass `"intrinsic"`, `"ambient"`, or an integer dimension to override it explicitly. `Geometry.centroid()` uses the same uniform intrinsic measure.

### Canonical boundary topology

The explicit-geometry layer includes `PolygonalSet` for disconnected polygonal sets with holes and `Polyhedron` for general 3-D boundary topology with outer and cavity shells. `PolyhedralShell` validates closed oriented manifold incidence, while `deduplicate_indexed_vertices()` provides exact coordinate
normalization.

### Reduced algebraic varieties

Algebraic geometry operations use certified reduced ideals as their semantic boundary. `certified_radicalization()` reconstructs and verifies `sqrt(I)`, `irreducible_components()` returns components only after a complete certified minimal-prime decomposition, and `reduced_component_singular_loci()` applies the Jacobian criterion separately to those reduced components. `singular_locus()` now radicalizes by default, so a presentation such as `x**2 = 0` has the same smooth
line geometry as `x = 0`.

### Stratified singular geometry

`stratified_singular_geometry()` refines the reduced-variety layer into certified minimal-prime intersections, per-component intrinsic singular loci, constructible local-dimension strata, and disjoint singular strata. This distinguishes smooth component crossings from singularities intrinsic to an irreducible component and handles mixed-dimensional varieties without conflating component dimension with local dimension.

### Local branch geometry

`local_branch_geometry(...)` refines stratified singular geometry at an exact point. It reports the irreducible branches through the point, each branch's exact ideal-theoretic tangent cone and Zariski tangent space, tangent-dimension excess, tangent-cone degree/local multiplicity, and—when several branches meet—the intersection tangent cone and whether the smooth branches meet transversely.

### Algebraic-statistics geometry

`polynomial_constraints`, `active_constraints`, `relative_interior`, `relative_boundary`, `semialgebraic_tangent_cone`, and `local_dimension` expose the local geometry needed by inequality-constrained statistical models. `parameterization_geometry` then analyzes full-dimensional polynomial/rational parameter domains, including generic rank, exact image and generic-fiber dimensions, critical loci, and critical-value images.

- Added certified implied/redundant polynomial inequalities, nonnegative-combination/SOS certificate replay, and component-aware semialgebraic constraint descriptions for hypothesis-testing geometry.

Polyhedral geometry also includes canonical polygon/polyhedron normalization, deterministic pulling/placing/barycentric polytope triangulations, and conforming mixed-cell tetrahedralization based on shared global vertex identifiers.

The polyhedral stack includes canonicalization, winding-rule polygonal path filling, explicit region conversion, pulling/placing/barycentric polytope decomposition, shell-based nonconvex polyhedron semantics, and conforming global-ID tetrahedralization of mixed convex 3-cells.

### Generation, refinement, and convex hulls
`convex_hull` constructs exact finite-point hulls (including embedded lower-dimensional hulls); `random_polygon` and `random_polytope` provide seeded exact lattice examples; and `subdivide_triangular_faces` / `geodesic_refinement` refine triangular polyhedral shells without floating-point geometric decisions. Convex-polytope intersections use the structural Boolean fast path when applicable.

### Certified real quantifier elimination

`quantifier_eliminate(...)` is the primary real-QE interface.  It uses specialized methods first and falls back to certified CAD; `project_region(...)` is the corresponding existential projection operation. Use
`return_result=True` to retain method and backend provenance.

The certified QE dispatcher exploits independent variable-incidence blocks, equality-ideal preprocessing, RUR/virtual substitution, and CAD without confusing Zariski elimination with real existential projection. Reusable parameterized CAD objects support specialization without rebuilding the decomposition.

## Worked examples

The documentation includes a [worked example gallery](https://github.com/BhuvaneshBhatt/semialg/blob/main/docs/examples/index.md) with executable examples for projection/QE, exact ranges and optimization, topology, moments, singular geometry, convex operations, distances, and parameterized algebraic integration.

## Documentation

- [Benchmark and conformance suite](https://github.com/BhuvaneshBhatt/semialg/blob/main/docs/benchmarks.md) — published Wilson CAD and TTICAD corpora plus parameterized Gröbner families, includes commands for reproducible conformance/performance runs.

New users should start with **[Getting started](https://github.com/BhuvaneshBhatt/semialg/blob/main/docs/getting_started.md)** and then follow the task-oriented links in the **[documentation index](https://github.com/BhuvaneshBhatt/semialg/blob/main/docs/index.md)**.

Here is some conceptual material:

- [Introduction to semialgebraic geometry](https://github.com/BhuvaneshBhatt/semialg/blob/main/docs/semialgebraic_geometry/introduction.md)
- [Exactness and certification](https://github.com/BhuvaneshBhatt/semialg/blob/main/docs/concepts/exactness_and_certification.md)
- [Performance guide](https://github.com/BhuvaneshBhatt/semialg/blob/main/docs/guides/performance.md)
- [Errors and failure modes](https://github.com/BhuvaneshBhatt/semialg/blob/main/docs/guides/errors_and_failure_modes.md)
- [Symbol handling](https://github.com/BhuvaneshBhatt/semialg/blob/main/docs/guides/symbol_handling.md)
- [Region invariants](https://github.com/BhuvaneshBhatt/semialg/blob/main/docs/guides/region_invariants.md)
- [API overview](https://github.com/BhuvaneshBhatt/semialg/blob/main/docs/api_overview.md)
- [Toric and lattice algebra](https://github.com/BhuvaneshBhatt/semialg/blob/main/docs/toric_algebra.md)

A progressive executable tutorial is available in [`notebooks/semialg_demo.ipynb`](https://github.com/BhuvaneshBhatt/semialg/blob/main/notebooks/semialg_demo.ipynb).

## Development

```bash
ruff format .
ruff check .
pytest -m "not slow"
pytest -m slow --durations=20
pytest -m performance --durations=20
python scripts/verify_source_quality.py
mkdocs build --strict
pytest -q -m slow tests/test_documentation_example_gallery.py
```

See [architecture guide](https://github.com/BhuvaneshBhatt/semialg/blob/main/docs/architecture/design.md) for implementation structure.

## References

- Michel Coste, *An Introduction to Semialgebraic Geometry*.

## License

See [LICENSE](https://github.com/BhuvaneshBhatt/semialg/blob/main/LICENSE).
