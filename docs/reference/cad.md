# CAD and structured geometry reference
## Family contract

**Mathematical return.** CAD APIs return cylindrical cells/decompositions, exact cell bounds and samples, or structured connectivity/integration data derived from them.

**Exactness and certification.** Root order, sign/truth invariance, and cell reconstruction are exact on certified paths. Reduced projection is used only when its side conditions are established; otherwise the implementation falls back conservatively.

**Algorithm.** Projection constructs lower-dimensional control polynomials; lifting decomposes fibers into sections and sectors. Variable ordering can materially change cost.

**Complexity and limitations.** CAD is a general but potentially expensive backend. Prefer task-specific APIs unless the decomposition itself is needed. See [CAD concepts](../concepts/cad.md).




## Primary API overview

This table is the substantive coverage target for the primary APIs assigned to this reference page. Each entry states the API's primary role; the family contract and detailed sections below explain shared algorithms, exactness guarantees, and limitations. It is maintained together with `docs/reference/primary_api_manifest.toml`, and documentation tests require every root-level primary API to map here rather than merely appearing in the generated public index.

| API | Kind | Role / return |
|---|---|---|
| `CADOptions` | class | Options accepted by :func:`cad`. |
| `CADResult` | class | Public result for CAD requests. |
| `cad` | function | Compute a cylindrical algebraic decomposition for a real formula. |

## `cad(...)`

The primary entry point is imported directly from the package root:

```python
from semialg import cad
```

CAD implementation modules are available to expert users under
`semialg.cad_algorithms`; that package name is deliberately distinct from the
root `cad` function.

Constructs or queries a cylindrical algebraic decomposition for supported semialgebraic input. Depending on options/output, CAD can support decision procedures, topology, structured cells, bounds, and integration adapters.

### Understanding direct CAD output and `CADResult`

By default, `cad()` returns the representation selected by `output`. The default
`output="formula"` therefore returns the reconstructed exact SymPy Boolean formula:

```python
import sympy as sp
from semialg import cad

x, y = sp.symbols("x y", real=True)
condition = x**2 + y**2 <= 1

formula = cad(condition, [x, y])
cells = cad(condition, [x, y], output="cells")
membership = cad(condition, [x, y], output="function")
tree = cad(condition, [x, y], output="tree")
```

Use `return_result=True` when the decomposition, diagnostics, status, and all output
representations are needed together:

```python
result = cad(condition, [x, y], return_result=True)
assert result.status == "complete"
formula = result.formula
cells = result.cell_set
membership = result.as_function()
```

A `CADResult` has the following principal fields and views:

- `result.formula` is the exact Boolean formula reconstructed from the selected
  final CAD cells (or from the requested topological operation). It represents the
  selected set when `result.status == "complete"`.
- `result.cad` is the cylindrical algebraic decomposition used internally. It
  contains cells at every lifting level, projection data, and certified samples.
- `result.cell_set` is the subset of final CAD cells on which the input formula (or
  requested operation) is true. `result.cells` is a convenience alias for these
  selected cells.
- `result.function` is an optional `CADFunction` membership representation.
  `result.as_function()` constructs it from a complete result when necessary.
- `result.variables` records the CAD variable order. CAD is cylindrical with
  respect to this order: every level projects onto a cell at the previous level,
  and changing the order can dramatically change both the decomposition and cost.
- `result.status` reports whether the requested decomposition is complete. An
  `"unknown"` result is diagnostic and must not be interpreted as emptiness or
  universality.
- `result.diagnostics` records strategy, preprocessing, resource-limit, and related
  metadata when available.

The central CAD invariant is **truth invariance**: after projection and lifting, the
input formula has a constant truth value on each final cell. Consequently, for a
complete result, `cell_set` is exactly the union of final cells on which the input
formula (or requested operation) is true. Cell indices encode the cylindrical path
through successive lifting levels; they are decomposition identifiers rather than
Cartesian coordinates.

### Variable ordering is a first-class performance decision

CAD is unusually sensitive to variable order. Two orders describe the same
mathematical set but can produce projection towers with very different polynomial
degrees, coefficient sizes, numbers of projection factors, algebraic root
complexity, and final cell counts. On nontrivial problems these differences can
translate into orders-of-magnitude differences in runtime and memory. A simple
answer therefore does **not** imply that every CAD order is cheap.

`result.variables` records the actual cylindrical/lifting order. In that order,
the first variable is the base coordinate and the final variable is eliminated
first during Collins-style projection. If the order is supplied explicitly,
`semialg` respects it where the operation's logical structure permits. Automatic
planning uses structural incidence/degree information and, on selected candidates,
projection and bounded pilot-lifting estimates. Expert users can inspect candidate
orders with `semialg.heuristics.suggest_variable_order(...)` or
`suggest_cad_variable_order(...)`.

A good order is problem-dependent. Useful signals include:

- put variables that yield low-degree, sparse projection polynomials in favorable
  elimination positions;
- exploit equations that permit dimension reduction or equational-constraint CAD;
- prefer orders that avoid introducing unnecessary algebraic root functions in
  downstream tasks such as integration;
- compare candidate projection/cell estimates for expensive problem families rather
  than relying on variable names or the written order of a formula.

For example, the region

\[
0\le y\le1,\qquad y^2\le x\le y
\]

has especially simple cylindrical bounds in the order `(y, x)`: first
`0 <= y <= 1`, then `y**2 <= x <= y`. Using `(x, y)` can instead force a
description involving `sqrt(x)`. The CAD-driven region-integration planner therefore
searches alternative legal orders and scores the resulting certified cells/bounds
rather than always inheriting the caller's coordinate order.

Quantifiers restrict what may be reordered. Variables inside the same quantifier
block can often be permuted without changing meaning, but moving variables across
alternating `Exists`/`ForAll` blocks changes the logical formula and is not a legal
performance optimization. Likewise, free-variable order can affect the desired
shape of reconstructed output. Automatic planning must preserve those semantic
constraints.

When diagnosing an unexpectedly expensive CAD, inspect variable order early. See
[CAD concepts](../concepts/cad.md) for the projection/lifting reason and the
[Performance guide](../guides/performance.md) for practical tuning guidance.

The `output` option controls the direct return value whenever `return_result=False`:

- `output="formula"` returns a SymPy Boolean formula;
- `output="cells"` returns a `CellSet`;
- `output="function"` returns a `CADFunction`;
- `output="tree"` returns the CAD tree.

Setting `return_result=True` instead returns the full `CADResult`; its `output` field
records the requested representation, while the other representations remain
available from the structured result.

### Preprocessing resource guard

Supported semialgebraic syntax such as `Abs` and rational powers may be converted to
polynomial constraints using existential auxiliary variables.  `max_preprocess_aux_vars`
(default `3`) bounds how many such auxiliaries `cad()` will eliminate automatically;
set it to `None` to disable this guard.  If the guard is exceeded, a structured call
returns `status="unknown"` with diagnostic details.  Calls with `return_result=False`
or `strict=True` raise `ResourceLimitError` instead of fabricating a mathematical
formula.

## `generic_cad(...)`

Builds a parameter-aware/generic CAD. By default the direct return follows `output` (`"formula"`, `"cases"`, `"cells"`, or `"function"`), matching `cad()`. Set `return_result=True` for `GenericCADResult`, including exceptional-set cases and CAD diagnostics. Specialization preserves the defining polynomial identity of algebraic root sections.

## Structured cell APIs

- `extract_structured_cad_cells`
- `structured_cad_cells_to_vertical_bounds_2d`
- `extract_vertical_bounds_from_cad_2d`
- `cylindrical_solution_from_structured`

These expose typed cylindrical structure rather than requiring callers to parse arbitrary Boolean formulas.

## Connectivity

- `build_cad_adjacency_graph`
- `extract_cad_connectivity`

These derive adjacency/components from CAD cell information.

## Bounds and certificates

Public CAD types include `CADBound`, `CADCellBoundsCertificate`, `StructuredCADLevel`, `StructuredCADCell`, `StructuredCADCellDecomposition`, `CADResult`, and related result/certificate objects.

Algebraic section bounds can be represented by certified root functions rather than approximate decimal endpoints.

## Exactness

Certified lifting, root ordering, sector sampling, and topology decisions use exact algebraic comparisons. Numerical approximations may be used for presentation/exploration only where explicitly allowed.

See [CAD and QE](../cad_qe.md) and [Performance](../guides/performance.md).

## Primary CAD configuration and context

`CADOptions` collects the user-facing controls for CAD construction, including
output form, strategy, resource limits, diagnostics, and formula reconstruction.
`SemialgebraicContext` stores reusable exact problem state, while
`computation_context` scopes shared caches and exact-computation settings for a
sequence of related operations.

## Gröbner-structured finite-variety CAD

`strategy="auto"` may select `effective_backend="groebner-variety"` when common polynomial equalities define an exact zero-dimensional ideal. The returned decomposition is a variety sub-CAD, so `complete` is false even though formula reconstruction is exact for the requested set. `strategy="collins"` forces a full-space decomposition.

Structured diagnostics include `variety_only`, `equality_dimension`, and `quotient_dimension`. Coordinate eliminants are lazy; exact uncertainty during specialized lifting causes fallback to complete Collins CAD. Pure existential QE may use this backend, while universal or mixed prefixes require complete-space semantics.

For the algorithm, compatibility of coordinate roots, quotient dimension versus point count, fallback contract, and diagnostic examples, see [Gröbner variety CAD](../concepts/groebner_variety_cad.md).
