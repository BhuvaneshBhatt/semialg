# Which function should I use?

Several `semialg` operations can express the same mathematical problem. Prefer the most specific API: it communicates intent and may use a cheaper exact algorithm than manually encoding everything as generic quantifier elimination.

| Goal | Preferred API |
|---|---|
| Decide whether a polynomial formula has a real solution | `is_satisfiable` |
| Decide a quantified first-order statement | complete QE / `qe_by_complete_cad` or the high-level solve/reduce interface |
| Prove implication or equivalence | `implies`, `equivalent` |
| Find an exact witness | `find_instance_formula` / solving APIs |
| Eliminate coordinates from a set | `semialgebraic_projection` |
| Compute values attained by a scalar function | `function_range` |
| Find a global minimum/maximum | `semialgebraic_minimize`, `semialgebraic_maximize` |
| Return the entire optimizer locus | `argmin_set`, `argmax_set`, `extrema_set` |
| Compute an image under a polynomial map | `region_image` |
| Pull a set back under a map | `region_preimage` |
| Compute volume, area, length, or intrinsic measure | `semialgebraic_measure` |
| Integrate a function over a region | `integrate_over_region` |
| Test a set relation | `is_subset`, `is_equal`, `is_disjoint`, `intersects` |
| Find connected pieces | `connected_components` |
| Compute distance or nearest points | `distance_to_region`, `distance_between_regions`, `nearest_point`, `closest_points` |
| Analyze a variety locally | `singular_locus`, `tangent_space`, `tangent_cone` |
| Inspect a CAD itself | `cad`, `extract_structured_cad_cells` |

## QE versus a specialized operation

For a projection, one can always write an existential formula. If
\(S(x,y)\) is a region, its projection onto \(x\) is

\[
\{x:\exists y\,S(x,y)\}.
\]

Use `semialgebraic_projection` rather than manually constructing the quantifier unless the quantified formula itself is the object you want to study.

Similarly,

\[
\operatorname{range}_S f
=
\{t:\exists x\,(x\in S\land t=f(x))\},
\]

but `function_range` is the preferred API because it can exploit range-specific structure and returns a range-oriented result.

## Range versus optimization

Use `function_range(f, ...)` when you need **all attainable values**.

Use `semialgebraic_minimize` or `semialgebraic_maximize` when you need an extremum, attainment information, or optimizer witnesses. Do not compute a complete range merely to obtain one endpoint unless the range itself is also useful.

Use `argmin_set`/`argmax_set` when the optimizer can be positive-dimensional. For example, minimizing \(x^2\) on a rectangle may have an entire line segment of minimizers.

## Projection versus image

`semialgebraic_projection` forgets coordinates.

`region_image` applies a map. A projection is a special linear image, but the projection API is clearer and can avoid unnecessary graph construction.

## Measure versus integration

`semialgebraic_measure(S, variables)` computes the measure of \(S\).

`integrate_over_region(f, S, variables)` computes

\[
\int_S f.
\]

Measure is conceptually the special case \(f=1\), but use the dedicated measure API when that is what you mean.

## Decision predicates versus constructing sets

If you only need to know whether \(A\subseteq B\), call `is_subset(A, B, ...)`. Constructing `region_difference(A, B)` and then separately deciding emptiness is mathematically equivalent but less direct.

## When should I call CAD directly?

Usually only when you need:

- the cylindrical cells themselves;
- sample points or cell metadata;
- adjacency/connectivity data;
- projection/lifting diagnostics;
- direct control over CAD strategy.

For ordinary decision, optimization, geometry, or integration tasks, use the corresponding high-level API and let the planner choose an exact backend.

## Virtual substitution, RUR, or CAD?

You normally do **not** choose these manually.

- Structure-aware presolve first removes safe affine/linear structure.
- **Quadratic virtual substitution** is available as a low-degree QE and witness backend.
- Zero-dimensional algebraic systems can use rational univariate representation (RUR).
- Complete CAD remains the general polynomial fallback.

See [How semialg chooses an algorithm](../concepts/algorithm_selection.md).
