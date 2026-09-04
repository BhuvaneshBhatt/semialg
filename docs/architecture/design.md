# Architecture and design

## Core interfaces

`semialg` exposes high-level APIs for real semialgebraic reasoning, including
`cad`, quantifier elimination, solution sampling, exact optimization, function
ranges, region operations, and integration.

SymPy Boolean formulas and relational expressions are the primary public
formula representation. Most APIs accept either one formula or a collection of
constraints.

## Exactness policy

Algorithms prefer exact results or explicit failure over uncertified symbolic
or numerical guesses. Specialized solvers may answer inexpensive cases first;
CAD/QE remains the semantic fallback when a complete exact decision is needed
and the computation is tractable.

## Algebraic samples

CAD cells carry explicit rational or algebraic sample objects. Root isolation,
sample comparison, sign evaluation, and approximate display remain separate so
numerical presentation cannot leak into certification.

## Structural presolve and variable ordering

Before complete CAD/QE, `presolve.py` performs only equivalence-preserving structural reductions whose exceptional cases are explicit. It can substitute an innermost existential variable from an affine equality when the divisor is a nonzero symbol-free constant, eliminate constant-coefficient linear existential inequalities by Fourier-Motzkin, and report independent variable-incidence blocks. It never divides by a parameter-dependent expression.

CAD ordering is also quantifier-aware: free variables and variables within one homogeneous quantifier block may be reordered, but an ordering heuristic never crosses an `exists`/`forall` boundary. Brown-style scoring is the automatic low-overhead choice; exhaustive projection-set scoring is an explicit diagnostic option for small systems.

## CAD backends

A conservative Collins-style CAD is the complete baseline. Reduced McCallum,
Lazard, partial-CAD, and formula-aware paths are used when their side conditions
and certification checks succeed. Equational constraints can reduce projection
and lifting work, while the complete path remains available when reduced
reasoning cannot be certified.

## Computation contexts and caching

Each top-level operation creates or reuses an `ExactComputationContext`.
Transient caches share projection, root, sign, comparison, specialization, and
RUR work among nested algorithms. Bounded process-wide caches provide reuse
between independent calls.

## Semantic queries

Simplification and validation routines often ask exact semantic questions such
as whether a constraint is redundant, a region is empty, or two formulas are
equivalent. This is more robust than purely syntactic rewriting, although it can
be more expensive.

## Function ranges

`function_range` treats a range as a semialgebraic image. For an expression `f`
and domain `C`, the solver introduces a value variable and eliminates the source
variables from a graph relation. Common expressions such as `Abs`, square roots,
`Min`, `Max`, and `Piecewise` have guarded exact graph encodings.

## Region integration

Region integration is layered around exact geometry:

```text
recognized region structure
  -> certified cylindrical or explicit bounds
  -> exact integral pieces
  -> symbolic or numeric evaluation when requested
  -> intrinsic Hausdorff measure for certified regular strata
```

Geometric reduction and antiderivative evaluation are reported separately when
symbolic integration cannot close the final integral.

## Internal module organization

Shared exact operations live in small, direct helper modules rather than being
reimplemented by each high-level subsystem:

- `normalization.py` resolves formulas, symbols, parameters, and bounds while
  preserving SymPy symbol identity and assumptions.
- `relations.py` normalizes polynomial relations and constructs relations to
  zero without duplicating operator handling.
- `interval_decomposition.py` owns exact one-dimensional breakpoint extraction,
  interval sampling, and truth decomposition used by measure, integration, and
  optimization.
- `cad/polynomial_utils.py` provides the canonical polynomial key shared by CAD
  projection and lifting code.
- `solve/integer/formula_utils.py` contains common conjunction splitting and
  exact integer-root helpers used by the integer-solving strategies.
- `planner/features.py` provides the stable feature signature used by planner
  strategy memory.

Optimization is organized by responsibility. `optimization_results.py` contains
result and policy models, `optimization_geometry.py` contains polynomial-locus
and geometric helpers, and `optimization_active_sets.py` contains active-set and
KKT construction. Shared projected positive-dimensional critical loci live in
`_critical_loci.py` so optimization and geometry queries use one KKT/singular-locus
implementation. `optimization.py` remains the stable public facade and core
optimizer, while `_optimization_parametric.py` owns parameter-stratified optimum
relations. Function-range orchestration remains in `_optimization_range.py`, with
fast special cases in `_range_special_cases.py` and graph/image QE in
`_function_graph_image.py`. The split is intentionally procedural: algebraic and
CAD inner loops do not gain strategy-object or method-dispatch layers.

Region integration separates the public integration surface from focused
implementation modules. `region_integrate.py` owns ambient reduction and public
entry points; `region_integral_results.py` owns result records,
`_region_integrate_intrinsic.py` owns intrinsic/Hausdorff-dimension handling,
and `_region_integrate_parametric.py` owns parameter-stratified integration.

`solve_semialgebraic` follows the same principle. Input normalization, parameter
analysis, and trivial-result construction are pure helpers, while the public
function remains the orchestration boundary for strategy selection and result
assembly.

## Public import surface

The package root uses a single declarative lazy-export registry. Accessing a
public name imports its owning module on first use and stores the resolved object
in the package module. The registry is the sole routing table for the flat API,
so `__all__` and lazy imports share one source of truth.

## Structural facades

The large public workflow modules are facades rather than monolithic implementation files. Function-property partition/property logic, optimization certification/specializations, decision predicates/solving, and region-integration geometry/CAD reduction live in focused private modules behind the documented public facade paths.

## Performance-sensitive organization

Code-sharing abstractions are kept outside root-isolation, CAD lifting,
algebraic comparison, projection, and RUR inner loops unless profiling shows
that an abstraction is neutral. Structural changes should be benchmarked with
warmed exact operations as well as fresh-process import/API access. Shared helpers are
preferred when they remove duplicate exact work or centralize correctness rules;
object-oriented indirection is not introduced solely for organizational
uniformity.

## Testing principles

Tests prefer semantic equivalence and exact mathematical invariants over string
comparisons. Public APIs are tested directly, while expensive CAD examples are
marked separately so routine feedback remains fast.

## Region module layering

The symbolic and CAD region APIs separate exact symbolic state from numerical
geometry while exposing a compact public namespace. See
[Code quality](code_quality.md#focused-implementation-modules) for the module map
and exception-boundary rules.

## Shared function graph and exact algebraization

Function-domain and function-range analysis share `semialg.function_graph`,
which constructs exact real semialgebraic graphs for supported algebraic heads.
This keeps domain and image semantics aligned, especially for branch-sensitive
rational powers.  Ordinary SymPy `Pow` keeps principal-branch semantics, while
canonical forms produced by `sympy.real_root` are recognized as explicit real
roots. The shared graph layer also supports `Abs`, `sign`, `Heaviside`, finite
`Min`/`Max`, finite `Piecewise`, and semialgebraic Boolean composition.

A preprocessing layer may uniquely back-substitute algebraic equalities before
range analysis.  A separate exact algebraization layer can then convert supported
commensurate trigonometric or exponential/hyperbolic one-variable problems into
finite semialgebraic problems.  Each accepted transformation carries its side
conditions (unit-circle or positive-exponential constraints); transformations
that would lose dependence information are rejected instead of used as
relaxations.

### Structural parametric geometry

Bounded structured regions expose `ParametricCover` objects made of
`ParametricChart` values. These charts are reusable geometry metadata: they
carry the parameter domain and mapping independently of CAD. Region dimension
uses a chart only when domain dimension and map rank are certified; otherwise
it delegates to the complete CAD implementation. Rich boundary queries follow
the same reuse principle by returning the boundary CAD together with exact
cell-level inclusion and active-constraint metadata.

Affine analysis and generic map degree are separate structural layers. Affine
properties are decided by exact linear algebra, while map degree is the generic
algebraic fiber degree of a rational parametrization. Neither abstraction folds
semialgebraic domain restrictions into a generic algebraic statement.
