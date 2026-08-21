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
and geometric helpers, `optimization_active_sets.py` contains active-set and KKT
construction, and `optimization.py` coordinates the public solver and exact
certification paths. The split is intentionally procedural: algebraic and CAD
inner loops do not gain strategy-object or method-dispatch layers.

`solve_semialgebraic` follows the same principle. Input normalization, parameter
analysis, and trivial-result construction are pure helpers, while the public
function remains the orchestration boundary for strategy selection and result
assembly.

## Public import surface

The package root uses cached lazy exports. Accessing a public name imports its
owning module on first use and stores the resolved object in the package module.
This preserves a convenient flat API without maintaining forwarding wrappers or
adding an extra wrapper call to every invocation.

## Performance-sensitive refactoring

Code-sharing abstractions are kept outside root-isolation, CAD lifting,
algebraic comparison, projection, and RUR inner loops unless profiling shows
that an abstraction is neutral. Refactors should be benchmarked with warmed
exact operations as well as fresh-process import/API access. Shared helpers are
preferred when they remove duplicate exact work or centralize correctness rules;
object-oriented indirection is not introduced solely for organizational
uniformity.

## Testing principles

Tests prefer semantic equivalence and exact mathematical invariants over string
comparisons. Public APIs are tested directly, while expensive CAD examples are
marked separately so routine feedback remains fast.
