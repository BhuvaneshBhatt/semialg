# semialg

`semialg` is a Python package for exact symbolic computation with real polynomial and semialgebraic conditions. It combines cylindrical algebraic decomposition (CAD), quantifier elimination (QE), exact algebraic-number methods, and specialized solvers for decision problems, solving, regions, optimization, and integration.

The package has certified paths that prefer an exact answer (or an explicit conservative failure) over silently treating a numerical approximation as proof.

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

First-order formulas can also be built directly with semialg's symbolic quantifiers:

```python
from semialg import Exists, ForAll
from semialg.solve import reduce_complete_expr

statement = ForAll(x, Exists(y, sp.Eq(x + y, 0)))
reduce_complete_expr(statement)
# True
```

Optimization and integration use the same exact-first model:

```python
from semialg import semialgebraic_minimize, semialgebraic_measure

opt = semialgebraic_minimize(x**2 + y**2, x + y >= 1, [x, y], return_result=True)
opt.value
# 1/2

semialgebraic_measure(x**2 + y**2 <= 1, [x, y])
# pi
```

Composite Boolean regions are interpreted geometrically through CAD rather than atom-by-atom rewriting. The same disjoint CAD-cell decomposition is used by multidimensional measure and integration, so overlapping unions are not double-counted and internal Boolean seams are not reported as boundaries.

The canonical symbolic region API is `SemialgebraicRegion`; you should build a `CADRegion` only when repeated point-location, topology, cell-complex, or CAD-integration queries justify reusing one decomposition:

```python
from semialg import SemialgebraicRegion

disk_region = SemialgebraicRegion(x**2 + y**2 <= 1, (x, y))
cad_region = disk_region.as_cad_region()

location = cad_region.locate_point((0, 0))
complex_ = cad_region.cell_complex()
assert location.selected
assert complex_.euler_characteristic() == 1
```

See the [region representation guide](https://github.com/BhuvaneshBhatt/semialg/blob/main/docs/guides/region_representations.md), [CAD reuse guide](https://github.com/BhuvaneshBhatt/semialg/blob/main/docs/guides/cad_reuse.md), and [exact-versus-numerical guide](https://github.com/BhuvaneshBhatt/semialg/blob/main/docs/guides/exact_vs_numerical_regions.md). The [API surface policy](https://github.com/BhuvaneshBhatt/semialg/blob/main/docs/architecture/api_surface.md) distinguishes primary root-level names from expert submodule APIs and internal implementation details.

Parameter-dependent region integrals can return certified guarded answers:

```python
a = sp.Symbol("a", real=True)
from semialg import integrate_over_region

parametric = integrate_over_region(
    x, (x >= 0) & (x <= a), [x],
    parameters=[a], return_stratified=True,
)
parametric.select({a: 2})
# 2
parametric.select({a: -1})
# 0
```

Parametric optimization and range APIs keep exact quantified relations by default. Pass `eliminate_quantifiers=True` together with `return_stratified=True` when a quantifier-free relation is worth the additional CAD/QE cost. Common affine one-dimensional parameter families are reconstructed directly before a second CAD/QE is attempted.

The complete-QE path also has a conservative structural presolver and automatic CAD ordering. Safe affine equalities and linear innermost existential blocks are eliminated before CAD; Brown-style ordering is applied only within semantically interchangeable free/quantifier blocks. Projection-set scoring is available explicitly through `suggest_cad_variable_order(..., strategy="projection")`.

```python
from semialg.heuristics import suggest_cad_variable_order
from semialg.presolve import presolve_semialgebraic

presolved = presolve_semialgebraic(
    sp.And(sp.Eq(y, x + 1), y > 0), [x, y], eliminate=[y]
)
# presolved.formula == (x > -1)

score = suggest_cad_variable_order([x**2 + y**4, x*y + 1], [x, y])
score.order
```

Univariate parameter-dependent integrals can use algebraic CAD root functions as moving endpoints. For example, the measure of `x**2 <= a` is stratfied exactly as zero for `a <= 0` and the distance between the two ordered roots for `a > 0`; specialization gives `4` at `a=4` and `6` at `a=9`.

High-level geometric queries compose the same certified backends:

```python
from semialg import (
    argmin_set, bounding_box, centroid, connected_components,
    diameter, distance_to_region, is_compact, minkowski_sum,
    semialgebraic_image, singular_locus, tangent_cone, width,
)

u = sp.Symbol("u", real=True)
semialgebraic_image(x, (x >= -1) & (x <= 2), [x], image_variables=[u])
# (u >= -1) & (u <= 2)

bounding_box((x >= -1) & (x <= 2), [x])
# {x: (-1, 2)}

distance_to_region((2,), (x >= 0) & (x <= 1), [x])
# 1

singular_locus(y**2 - x**3, [x, y])
tangent_cone(y**2 - x**3, {x: 0, y: 0}, [x, y]).formula
# Exact ideal-theoretic cone for multiple generators; cancellation terms are retained.
tangent_cone((x**2 + y**3, x**2 - y**3), {x: 0, y: 0}, [x, y]).ideal_generators
# Eq(d_y**2, 0)
```

Derived geometry APIs expose common exact constructions without duplicating algorithms:

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

The same layer includes set relations (`is_subset`, `is_equal`, `is_disjoint`, `intersects`), level/sublevel/superlevel sets, nearest/closest points, moment/covariance/inertia matrices, affine transforms and Minkowski sums, support functions, directional widths, and algebraic `is_smooth` / `is_singular` / `tangent_dimension` conveniences. Operations such as Minkowski sums and general affine images remain exact existential-QE problems and can therefore inherit CAD complexity.

Region topology uses CAD semantics even for atomic polynomial inequalities; compact direct formulas are used only after cell-wise equivalence verification. `region_dimension` is likewise exact and is read from selected CAD-cell dimensions. Free parameters are preserved by `semialgebraic_image` when source variables are explicit, and all geometry point/string inputs use contextual symbol resolution.

`is_path_connected` is an exact decision through CAD connectivity. `path_between` returns an explicit piecewise-linear path in the certified 1-D case and a certified CAD cell/connector chain in higher dimensions; full Canny/Basu-Pollack-Roy roadmap parameterization is outside this API.

## What semialg provides

- **Decision and QE:** satisfiability, tautology, implication, equivalence, CAD and quantifier elimination.
- **Specialized real algebraic decisions:** certified positive-dimensional equality feasibility, replayable algebraic decomposition, exact polynomial nonnegativity, and one-sided exact witness search.
- **Solving and sampling:** an exact solver with affine presolve, Boolean/incidence decomposition, univariate and linear fast paths, zero-dimensional RUR, Groebner presolve, and CAD/QE fallback, plus structured witnesses and samples.
- **Algebraic roots and parameters:** exact root isolation, certified algebraic root functions, root classification, and parameter-stratified results.
- **Regions:** Boolean region operations, exact set relations/properties, topology, CAD connected components, transforms, standard geometric regions, and CAD-derived geometry.
- **Geometric queries:** exact projection/image/preimage/fibers, extrema sets, bounding boxes, distances/diameter, support/width, transforms/Minkowski sums, staged exact convexity certificates, connectivity/property predicates, Euler characteristic, moments, singular loci, tangent spaces, and tangent cones.
- **Optimization and ranges:** exact polynomial optimization, KKT/active-set analysis, global certification, and semialgebraic image/range computation.
- **Integration and moments:** ambient and intrinsic measure, region integrals, moments, centroids, and covariance.

See [feature matrix](https://github.com/BhuvaneshBhatt/semialg/blob/main/docs/feature_matrix.md) for a more detailed capability summary and [limitations](https://github.com/BhuvaneshBhatt/semialg/blob/main/docs/limitations.md) for important scope boundaries.

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

See [`docs/applications.md`](https://github.com/BhuvaneshBhatt/semialg/blob/main/docs/applications.md). Core operations such as `function_range`, `semialgebraic_measure`, integration, optimization, CAD, and QE remain in the core namespace rather than being duplicated under `applications`.

## Exact vs certified

A result being symbolic is not by itself a certificate. semialg distinguishes exact representations, certified conclusions, candidate/heuristic information, and explicitly numerical approximations. See [Exactness and certification](https://github.com/BhuvaneshBhatt/semialg/blob/main/docs/concepts/exactness_and_certification.md).

## Documentation

New users should start with [Getting started](https://github.com/BhuvaneshBhatt/semialg/blob/main/docs/getting_started.md) and then follow the task-oriented links in the [documentation index](https://github.com/BhuvaneshBhatt/semialg/blob/main/docs/index.md).

Important conceptual material:

- [Introduction to semialgebraic geometry](https://github.com/BhuvaneshBhatt/semialg/blob/main/docs/semialgebraic_geometry/introduction.md)
- [Exactness and certification](https://github.com/BhuvaneshBhatt/semialg/blob/main/docs/concepts/exactness_and_certification.md)
- [Performance guide](https://github.com/BhuvaneshBhatt/semialg/blob/main/docs/guides/performance.md)
- [Errors and failure modes](https://github.com/BhuvaneshBhatt/semialg/blob/main/docs/guides/errors_and_failure_modes.md)
- [Symbol handling](https://github.com/BhuvaneshBhatt/semialg/blob/main/docs/guides/symbol_handling.md)
- [Region invariants](https://github.com/BhuvaneshBhatt/semialg/blob/main/docs/guides/region_invariants.md)
- [API overview](https://github.com/BhuvaneshBhatt/semialg/blob/main/docs/api_overview.md)

A progressive executable tutorial is available in [`notebooks/semialg_demo.ipynb`](https://github.com/BhuvaneshBhatt/semialg/blob/main/notebooks/semialg_demo.ipynb).

## Worked examples

The documentation includes a [worked example gallery](https://github.com/BhuvaneshBhatt/semialg/blob/main/docs/examples/index.md) with executable examples for projection/QE, exact ranges and optimization, topology, moments, singular geometry, convex operations, distances, and parameterized algebraic integration.

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

See the [architecture guide](https://github.com/BhuvaneshBhatt/semialg/blob/main/docs/architecture/design.md) for implementation structure.

## References

- Michel Coste, *An Introduction to Semialgebraic Geometry*.

## License

See [LICENSE](https://github.com/BhuvaneshBhatt/semialg/blob/main/LICENSE).
