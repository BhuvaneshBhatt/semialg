# semialg

`semialg` is a Python package for exact computational semialgebraic  and real algebraic geometry. It combines cylindrical algebraic decomposition (CAD), quantifier elimination (QE), exact real-algebraic arithmetic, polynomial-system and algebraic-decomposition methods, and specialized algorithms for solving, decision problems, geometric analysis, optimization, and integration.

## Why use `semialg`?

`semialg` is especially useful when a problem mixes polynomial equations and inequalities, Boolean conditions, parameters, (exact) algebraic numbers, or global questions such as feasibility, equivalence, projection, optimization, topology, or measure. Specialized algebraic and geometric methods handle inexpensive cases first; CAD and QE provide fallback machinery for supported real-polynomial formulas. Cheap Boolean certificates, affine presolve, exact linear feasibility, univariate reduction, finite algebraic solving, and local geometric structure are used when their applicability conditions can be certified.

A result being symbolic is not by itself a certificate. `semialg` distinguishes between  exact representations, certified conclusions, candidate/heuristic information, and explicitly numerical approximations. See [Exactness and certification](https://github.com/BhuvaneshBhatt/semialg/blob/main/docs/concepts/exactness_and_certification.md).

If you primarily need floating-point nonlinear optimization, general numerical transcendental solving, or large numerical geometry (instead of exact semialgebraic computation), established numerical libraries will usually be a better fit. Consider SciPy for numerical optimization and nonlinear systems, SymPy/mpmath for numerical transcendental root finding, Shapely for planar computational geometry, and Trimesh for numerical 3-D mesh geometry. See the [capability matrix](https://github.com/BhuvaneshBhatt/semialg/blob/main/docs/feature_matrix.md) and [limitations](https://github.com/BhuvaneshBhatt/semialg/blob/main/docs/limitations.md) before choosing an algorithm.

## Install

```bash
python -m pip install semialg
```

Python 3.11+ is required. Dependencies include SymPy `>=1.14,<2`, python-flint `>=0.9,<1`, algroots `>=0.4.1`, and exprtest `>=0.1.0`; `pip` installs these automatically.

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

For a decision query, `semialg` first tries inexpensive structure before falling back to general elimination machinery:

![Decision computation flow](https://raw.githubusercontent.com/BhuvaneshBhatt/semialg/main/docs/assets/decision-flow.svg)

The [computation-flow guide](https://github.com/BhuvaneshBhatt/semialg/blob/main/docs/concepts/computation_flows.md) shows the corresponding CAD, optimization, and integration pipelines.

First-order formulas can be built directly with `semialg`'s symbolic quantifiers:

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

## Local and approach geometry

`semialg` can certify whether a point is approachable through a semialgebraic region, preserve correlations in polynomial/rational map images, and resolve local algebraic branches and incident strata near singular points. These operations are useful when downstream symbolic analysis must distinguish genuine approach geometry from independent coordinate bounds.

```python
from semialg import correlated_map_image, point_in_closure

t, u, v = sp.symbols("t u v", real=True)
image = correlated_map_image((t, t**2), sp.And(t >= -1, t <= 1), (t,), image_variables=(u, v))
assert sp.simplify(image.formula.subs({u: sp.Rational(1, 2), v: sp.Rational(1, 4)})) is sp.true
assert sp.simplify(image.formula.subs({u: sp.Rational(1, 2), v: sp.Rational(1, 3)})) is sp.false

x, y = sp.symbols("x y", real=True)
punctured = sp.And(x**2 + y**2 < 1, sp.Ne(x**2 + y**2, 0))
assert point_in_closure(punctured, (0, 0), (x, y))
```

Local germs expose branch and sign information directly:

```python
from semialg import local_components, local_sign_strata

crossing = sp.Eq(y**2, x**2)
branches = local_components(crossing, (0, 0), (x, y))
assert len(branches) == 4

signs = local_sign_strata((x, y), crossing, (0, 0), (x, y))
assert all(piece.signs in {(-1, -1), (-1, 1), (0, 0), (1, -1), (1, 1)} for piece in signs)
```

Exact real witnesses returned by `find_instance` also respect sign assumptions carried by supplied SymPy symbols.

For correlated directional images, closure semantics, local germs, parameter strata, and proof diagnostics, see [Approach and local geometry](https://github.com/BhuvaneshBhatt/semialg/blob/main/docs/guides/approach_and_local_geometry.md).

For certified approach geometry used in asymptotic analysis, including germs, bounds, and weighted charts, see [Local geometry for asymptotic analysis](https://github.com/BhuvaneshBhatt/semialg/blob/main/docs/guides/local-asymptotic-geometry.md).

Reusable [chart geometry certificates](https://github.com/BhuvaneshBhatt/semialg/blob/main/docs/chart_certificates.md) check scalar/vector level maps, inverses, Gram rank, coverage, overlap and null sets, and identify bounded parameter-event candidates. Integration callers are responsibile for convergence and exceptional-fiber evaluation.

## Core workflows

### Cylindrical algebraic decomposition

`cad(...)` computes an exact cylindrical algebraic decomposition and reconstructs the selected cells as a semialgebraic formula by default. Variable order is explicit and can change the shape of the decomposition without changing the represented set.

```python
from semialg import cad

x, y = sp.symbols("x y", real=True)

cad(x**2 + y**2 < 1, (x, y))
# (x > -1) & (x < 1) & (y + sqrt(1 - x**2) > 0) & (y - sqrt(1 - x**2) < 0)

# Intervals for univariate polynomials
cad((x - 1)*(x - 2)**2*(x - 3)*(x - 4)*(x - 5) > 0, (x,))
# (x > 5) | (x < 1) | ((x > 3) & (x < 4))

upper_half = (x**2 + y**2 < 1) & (y > x)
cad(upper_half, (x, y))
cad(upper_half, (y, x)) # Note the different variable order

# Mixed-dimensional Boolean geometry: a closed disk together with a line.
cad((x**2 + y**2 <= 1) | sp.Eq(x, 1), (x, y))
```

Quantified projections combine QE with CAD. For example, projecting out `y` and then decomposing the remaining `x`-set gives an exact interval:

```python
from semialg import Exists, quantifier_eliminate

projected = quantifier_eliminate(
    Exists(y, (x**2 + y**2 < 1) & (y > x))
)
cad(projected, (x,))
# (x > -1) & (x < sqrt(2)/2)
```

Full Lazard projection and valuation-aware lifting handle degree drops and supported nullification cases. Native CAD proof replay checks the projection and lifting geometry independently of stored completion flags (for supported certificate formats). The expert interfaces are in `semialg.cad_algorithms`:

```python
from semialg.cad_algorithms import decomp_lazard_complete, verify_cad

polynomial = y**2 - x
native = decomp_lazard_complete((polynomial,), (x, y))
proof = verify_cad(native, polynomials=(polynomial,), return_result=True)
assert proof.valid
```

For recursive equational-constraint projection and the limits of certificate replay, see [Native CAD proofs](https://github.com/BhuvaneshBhatt/semialg/blob/main/docs/native_cad_proofs.md).

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

For parameter-dependent geometry, `region_family` identifies exact parameter strata, including empty fibers and changes in fiber dimension or (for scalar fibers) component count. `bounding_cover` gives a potentially loose exact outer box without optimization, while `optimization_cover` retains a structural parametrization including boundary and lower-dimensional pieces:

```python
from semialg import Ball, bounding_cover, optimization_cover, region_family

a = sp.Symbol("a", real=True)
ball = Ball((1, 2), 3)
assert bounding_cover(ball, (x, y)) == {x: (-2, 4), y: (-1, 5)}
assert optimization_cover(ball, (x, y)).charts
family = region_family(sp.Eq(x**2, a), (x,), (a,))
assert sp.simplify(family.empty_condition.subs(a, -1)) is sp.true
```

`bounding_cover` computes exact, potentially loose outer coordinate bounds for standard geometric regions, supported Boolean combinations, and transformed/parametrized regions whose coordinate maps admit interval bounds. `optimization_cover` constructs finite parametrization charts for supported bounded geometry, including balls, spheres, ellipsoids, simplices, polygons, finite point sets, unions, and transformed regions. It retains boundaries and lower-dimensional pieces. General formula-defined regions require a different approach, such as `bounded_parametric_cover` with explicit finite clipping bounds.

`region_family` analyzes polynomial semialgebraic families using certified CAD. It reports parameter strata, empty fibers, and fiber dimensions. One-dimensional fibers additionally have exact component counts and graph/band topology; polynomial ball and sphere families have structural topology invariants. General higher-dimensional fibers have cylindrical decomposition data, but arbitrary component counts, Betti numbers, and Hardt trivializations are not currently provided.

See [Optimization geometry and parameter families](https://github.com/BhuvaneshBhatt/semialg/blob/main/docs/optimization_geometry.md).

### Quantifier elimination and CAD planning

`quantifier_eliminate(...)` is the operation for eliminating quantified real variables while preserving free parameters.

Automatic quantifier elimination can use quantifier-guided partial CAD, pruning existential branches after success and universal branches after failure. Reduced CAD exploits recursively justified equational constraints where its projection hypotheses can be certified, with safe fallback otherwise. `cad(...)` and `quantifier_eliminate(...)` accept `max_cells` and cooperative `timeout` limits; resource exhaustion never certifies a partial result as complete. See [CAD and QE](https://github.com/BhuvaneshBhatt/semialg/blob/main/docs/cad_qe.md) and [Parameter-cell QE and budgets](https://github.com/BhuvaneshBhatt/semialg/blob/main/docs/tutorials/parameter_cell_qe.md).

```python
from semialg import Exists, ForAll, quantifier_eliminate

a, b, c = sp.symbols("a b c", real=True)

# Project x**2 = a onto the real parameter axis.
condition = quantifier_eliminate(Exists(x, sp.Eq(x**2, a)))
# a >= 0

quantifier_eliminate(Exists((x, y), x**2 + y**2 < 1))
# True

quantifier_eliminate(Exists(x, sp.Eq(a*x + b, 0)))
# Ne(a, 0) | (Eq(a, 0) & Eq(b, 0))

quantifier_eliminate(Exists(x, sp.Eq(a*x**2 + b*x + c, 0)))
# (Eq(a, 0) & Eq(b, 0) & Eq(c, 0)) | \
# (Eq(a, 0) & Ne(b, 0) & Eq(a*c, 0)) | \
# (Ne(a, 0) & (4*a*c - b**2 <= 0))

# A nested universal/existential problem: for which a,b does every
# quadratic constant term c admit a real root?
quantifier_eliminate(
    ForAll(c, Exists(x, sp.Eq(a*x**2 + b*x + c, 0)))
)
# (Eq(a, 0) & (b > 0)) | (Eq(a, 0) & (b < 0))

# Pure existential feasibility can collapse all the way to a truth value.
quantifier_eliminate(Exists(x, (x**2 < 2) & (x > 1)))
# True
```

The complete-QE path also has a structural presolver and automatic CAD ordering. Safe affine equalities and linear innermost existential blocks are eliminated before CAD; automatic Brown/chordal structural ordering is applied only within semantically interchangeable free/quantifier blocks.  Projection-set scoring is available explicitly through `suggest_cad_variable_order(..., strategy="projection")`.

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

Univariate parameter-dependent integrals can use algebraic CAD root functions as moving endpoints. For example, the measure of `x**2 <= a` is stratfied exactly as zero for `a <= 0` and the distance between the two ordered roots for `a > 0`; specialization gives `4` at `a=4` and `6` at `a=9`.

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

`is_path_connected` is an exact decision through CAD connectivity. `path_between` returns an explicit piecewise-linear path in the certified 1-D case and a certified CAD cell/connector chain in higher dimensions; full Canny/Basu-Pollack-Roy roadmap parametrization is not implemented yet.

For supported quadratic inequalities and Cartesian-product regions, structural topology rules can determine closures, interiors, boundaries, and some connectivity properties before general CAD. See [Exact rules before CAD](https://github.com/BhuvaneshBhatt/semialg/blob/main/docs/structural_rules.md).

### Exact optimization specializations

Affine objectives on supported ellipsoids and univariate polynomial objectives can use exact structural methods before general CAD. Optimization results preserve whether an infimum or supremum is actually attained:

```python
from semialg import semialgebraic_minimize

ellipse = (x - 2)**2 + 4*(y + 1)**2 <= 9
minimum = semialgebraic_minimize(3*x + 4*y + 5, ellipse, (x, y), return_result=True)
assert minimum.attained and minimum.certified
punctured = semialgebraic_minimize(x**2, sp.Ne(x, 0), (x,), return_result=True)
assert punctured.value == 0 and not punctured.attained
```

These shortcuts require supported exact polynomial inputs, otherwise the ordinary exact optimization backends apply. See [Structural rules](https://github.com/BhuvaneshBhatt/semialg/blob/main/docs/structural_rules.md).

## Primary and specialized APIs

The package root is the everyday mathematical interface. Specialized algorithms and certificate-building functions live in their owning namespaces such as  `semialg.parameters.root_count_conditions`, `semialg.roadmaps.roadmap`, `semialg.topology.semialgebraic.triangulate_region`, `semialg.parametric_geometry.bounded_parametric_cover`, and `semialg.map_degree.parametric_map_degree`. See [API namespaces](https://github.com/BhuvaneshBhatt/semialg/blob/main/docs/api_namespaces.md).

Semialgebraic function analysis is composition-aware for the supported graph fragment. Nested radicals/rational powers, `Abs`, `sign`, `Heaviside`, `Min`/`Max`, finite `Piecewise`, and real `re`/`im`/`conjugate` wrappers can be "algebraized" when their real-domain side conditions are certifiable; unsupported graphs are declined.

## What semialg provides

- **Decision and QE:** satisfiability, tautology, implication, equivalence, CAD and quantifier elimination.
- **Specialized real algebraic decisions:** certified positive-dimensional equality feasibility, replayable algebraic decomposition, polynomial nonnegativity, and one-sided exact witness search.
- **Solving and sampling:** an exact solver with affine presolve, Boolean/incidence decomposition, univariate and linear fast paths, zero-dimensional RUR, Groebner presolve, and CAD/QE fallback, plus structured witnesses and samples.
- **Algebraic roots and parameters:** exact root isolation, certified algebraic root functions, root classification, and parameter-stratified results.
- **Toric and lattice algebra:** saturated integer kernels, lattice and toric ideals, Laurent-monomial elimination, and exact Markov bases.
- **Algebraic statistics:** exact semialgebraic local geometry for polynomial statistical models, active/redundant constraints, tangent cones and local dimension, parametrisation geometry, toric/lattice constructions, and certificate-aware hypothesis-testing geometry.
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

The explicit-geometry layer includes `PolygonalSet` for disconnected polygonal sets with holes and `Polyhedron` for general 3-D boundary topology with outer and cavity shells. `PolyhedralShell` validates closed oriented manifold incidence, while `deduplicate_indexed_vertices()` provides exact coordinate normalization.

### Reduced algebraic varieties

Algebraic geometry operations use certified reduced ideals as their semantic boundary. `certified_radicalization()` reconstructs and verifies `sqrt(I)`, `irreducible_components()` returns components only after a complete certified minimal-prime decomposition, and `reduced_component_singular_loci()` applies the Jacobian criterion separately to those reduced components. `singular_locus()` radicalizes by default, so a presentation such as `x**2 = 0` has the same smooth
line geometry as `x = 0`.

### Stratified singular geometry

`stratified_singular_geometry()` refines the reduced-variety layer into certified minimal-prime intersections, per-component intrinsic singular loci, constructible local-dimension strata, and disjoint singular strata. This distinguishes smooth component crossings from singularities intrinsic to an irreducible component and handles mixed-dimensional varieties without conflating component dimension with local dimension.

### Local branch geometry

`local_branch_geometry(...)` refines stratified singular geometry at an exact point. It reports the irreducible branches through the point, each branch's exact ideal-theoretic tangent cone and Zariski tangent space, tangent-dimension excess, tangent-cone degree/local multiplicity, and—when several branches meet—the intersection tangent cone and whether the smooth branches meet transversely.

### Algebraic-statistics geometry

`polynomial_constraints`, `active_constraints`, `relative_interior`, `relative_boundary`, `semialgebraic_tangent_cone`, and `local_dimension` expose the local geometry needed by inequality-constrained statistical models. `parameterization_geometry` then analyzes full-dimensional polynomial/rational parameter domains, including generic rank, exact image and generic-fiber dimensions, critical loci, and critical-value images.

Semialgebraic constraint analysis identifies implied and redundant polynomial inequalities, supports independently verifiable nonnegative-combination and sum-of-squares (SOS) certificates, and provides component-aware constraint descriptions for hypothesis-testing geometry.

Polyhedral geometry also includes canonical polygon/polyhedron normalization, deterministic pulling/placing/barycentric polytope triangulations, and conforming mixed-cell tetrahedralization based on shared global vertex identifiers.

The polyhedral stack includes canonicalization, winding-rule polygonal path filling, explicit region conversion, pulling/placing/barycentric polytope decomposition, shell-based nonconvex polyhedron semantics, and conforming global-ID tetrahedralization of mixed convex 3-cells.

### Generation, refinement, and convex hulls
`convex_hull` constructs exact finite-point hulls (including embedded lower-dimensional hulls), `random_polygon` and `random_polytope` provide seeded exact lattice examples, and `subdivide_triangular_faces` / `geodesic_refinement` refine triangular polyhedral shells without floating-point geometric decisions. Convex-polytope intersections use the structural Boolean fast path when applicable.

### Certified real quantifier elimination

`quantifier_eliminate(...)` is the primary real-QE interface.  It uses specialized methods first and falls back to certified CAD; `project_region(...)` is the corresponding existential projection operation. Use `return_result=True` to retain method and backend provenance.

The certified QE dispatcher exploits independent variable-incidence blocks, equality-ideal preprocessing, RUR/virtual substitution, and CAD. Reusable parametrized CAD objects support specialization without rebuilding the decomposition.

`solvability_conditions(..., groebner_backend="native", groebner_limits=NativeF4Limits(...))` can preprocess outer conjunctive equalities using `algroots`'s verified number/parameter-field backends. Generic bases carry exact specialization guards, while original equalities and inequalities remain on exceptional fibers and the full formula goes through complete real elimination. Native limits cover equality preprocessing, not the subsequent CAD calculation. This option does not change the default dispatch or replace real arithmetic with a prime-field answer. See [algroots coefficient fields and specialization](https://github.com/BhuvaneshBhatt/algroots/blob/main/docs/coefficient-fields.md).

## Worked examples

The documentation includes a [worked example gallery](https://github.com/BhuvaneshBhatt/semialg/blob/main/docs/examples/index.md) with executable examples for projection/QE, exact ranges and optimization, topology, moments, singular geometry, convex operations, distances, and parametrized algebraic integration.

## Documentation

New users should start with [Getting started](https://github.com/BhuvaneshBhatt/semialg/blob/main/docs/getting_started.md) and then follow the task-oriented links in the [documentation index](https://github.com/BhuvaneshBhatt/semialg/blob/main/docs/index.md).

Here is some conceptual material:

- [Introduction to semialgebraic geometry](https://github.com/BhuvaneshBhatt/semialg/blob/main/docs/semialgebraic_geometry/introduction.md)
- [Comprehensive semialgebraic geometry overview](https://github.com/BhuvaneshBhatt/semialg/blob/main/docs/semialgebraic_geometry/index.md)
- [Region invariants](https://github.com/BhuvaneshBhatt/semialg/blob/main/docs/guides/region_invariants.md)
- [Toric and lattice algebra](https://github.com/BhuvaneshBhatt/semialg/blob/main/docs/toric_algebra.md)
- [API overview](https://github.com/BhuvaneshBhatt/semialg/blob/main/docs/api_overview.md)
- [Exactness and certification](https://github.com/BhuvaneshBhatt/semialg/blob/main/docs/concepts/exactness_and_certification.md)
- [Performance guide](https://github.com/BhuvaneshBhatt/semialg/blob/main/docs/guides/performance.md)
- [Errors and failure modes](https://github.com/BhuvaneshBhatt/semialg/blob/main/docs/guides/errors_and_failure_modes.md)
- [Symbol handling](https://github.com/BhuvaneshBhatt/semialg/blob/main/docs/guides/symbol_handling.md)
- [Benchmark and conformance suite](https://github.com/BhuvaneshBhatt/semialg/blob/main/docs/benchmarks.md) — published Wilson CAD and TTICAD corpora plus parametrized Gröbner families, includes commands for reproducible conformance/performance runs.

An executable tutorial is available in [`notebooks/semialg_demo.ipynb`](https://github.com/BhuvaneshBhatt/semialg/blob/main/notebooks/semialg_demo.ipynb).

## Development

Contributions – whether they are bug reports, suggestions, or code – are welcome!

```bash
ruff format .
ruff check .
pytest
mkdocs build --strict
```

The default suite includes computational and performance regression tests. Use `pytest -m slow --durations=20` or `pytest -m performance --durations=20` to inspect those subsets separately. Large [stress workloads](https://github.com/BhuvaneshBhatt/semialg/blob/main/docs/quality/stress_workloads.md) require `pytest -m stress tests/stress --durations=0`.

See [architecture guide](https://github.com/BhuvaneshBhatt/semialg/blob/main/docs/architecture/design.md) for implementation structure.

## References

- Michel Coste, *An Introduction to Semialgebraic Geometry*.

## License

`semialg` is licensed under the GNU General Public License version 3 only (GPL-3.0-only). See [LICENSE](https://github.com/BhuvaneshBhatt/semialg/blob/main/LICENSE).
