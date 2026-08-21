# Changelog

All notable changes to `semialg` are documented here.

## 0.2.0b1

**Pre-release status:** Beta. This release is intended for testing and evaluation before final 0.2.0 release. Public APIs and behavior could get compatibility-breaking corrections before the final release.

### Added

- Added first-class semialg `Exists`/`ForAll` integration across expression-level complete QE and transcendental state construction, plus prenex conversion helpers and `ParsedPrenexFormula.quantified_expr`.
- Formalized transcendental result coverage with `ResultSemantics`, including exact, subset, superset, window-scoped, periodic-window, and no-witness semantics.
- Polynomial model comparison with exact worst-case discrepancy, dominance, and equivalence certificates.
- Exact parameter-regime analysis for semialgebraic solvability and real-root-count changes.
- Exact polynomial/geometric probability with density nonnegativity and positive-mass certification.
- Lyapunov-function verification for polynomial continuous-time dynamics.
- Polynomial barrier-certificate verification with exact initial, unsafe-set, and boundary derivative obligations.
- Certified polynomial sensitivity/monotonicity analysis with exact derivative ranges.
- Exact constraint-redundancy analysis with nonredundancy witnesses.
- Feasible-set diagnostics with exact witnesses or irreducible conflicting constraint subsets.

- Added `semialg.applications` workflows for robust parameter/tolerance analysis, certified symbolic-math validation, exact optimization benchmark oracles, polynomial Hurwitz stability regions, polynomial safety/invariant verification, and polynomial response-surface analysis.

- Certified algebraic root functions with exact specialization, comparison, differentiation, regularity checks, and guarded symbolic presentation.
- Typed cylindrical CAD bounds and decomposition-level certification.
- Regular-stratum intrinsic integration with explicit singular-stratum classification.
- Exact multivariate polynomial optimization using stationary points, active-set/KKT systems, singular loci, exact candidate comparison, and CAD certification.
- Cost-controlled range certification and positive-dimensional critical-locus handling for higher-dimensional optimization.
- First-class parameter-stratified and conditional results, including parametric optimization and function-range relations.
- Solve-scoped exact-computation contexts and bounded caches for projection, root isolation, sign determination, algebraic comparisons, specialization, and RUR computations.
- Projection- and lifting-aware CAD variable-order planning, including algebraic-complexity estimates and bounded pilot lifting.
- Multi-level equational-constraint propagation and conservative partial-CAD pruning.
- Repository source-quality verification with `scripts/verify_source_quality.py`.
- First-class symbolic `Exists` and `ForAll` Boolean nodes with bound-variable tracking and capture-avoiding substitution; periodic transcendental reconstruction now uses explicit existential integer indices.

### Changed

- Aligned package metadata and the public `semialg.__version__` with the 0.2.0 release.
- Hardened the PyPI tag workflow so Ruff, source-quality checks, the default test matrix, the slow/performance suite, distribution building, and `twine check` must succeed before publishing.

- Reorganized documentation into a concise landing page, progressive getting-started material, conceptual exactness guidance, task-focused guides, and API-family reference pages; rewrote the demo notebook as an executable progressive tutorial.

- Consolidated exact one-dimensional interval decomposition, root breakpoint,
  interval sampling, and Boolean truth helpers into a shared implementation used
  by measure, integration, and optimization.
- Centralized formula, variable, parameter, and bound normalization around the
  package symbol-resolution layer.
- Centralized polynomial-relation parsing and zero-relation construction.
- Centralized CAD polynomial canonical keys so projection, lifting, and
  sign-invariance code use the same normalization rule.
- Consolidated repeated integer-solver formula helpers and planner feature
  signatures.
- Replaced package-level forwarding wrappers with cached lazy exports, so a
  public object resolves once and subsequent calls invoke the real function or
  class directly.
- Split optimization internals by responsibility into result models, geometry
  and locus utilities, active-set/KKT construction, and high-level
  orchestration without changing the optimization algorithms.
- Decomposed `solve_semialgebraic` request handling into smaller pure helpers for
  normalization, parameter analysis, and trivial-result construction.
- Kept CAD and algebraic hot loops direct; the shared helpers are concentrated
  at normalization, orchestration, and reusable exact-geometry boundaries.
- Refactored source, tests, and documentation to use descriptive names and current-state architecture terminology.
- Strengthened exactness throughout CAD, topology, interval handling, optimization, and integration: certified paths decline unsupported cases rather than treating finite-precision numerical approximations as proofs.
- Improved symbol-name resolution so string arguments reuse the actual SymPy symbols present in an expression and reject ambiguous same-name symbols.
- Reworked plotting adapters and package exports to remove runtime monkeypatch patterns and overridden definitions.
- Expanded optimization preprocessing with equality elimination and more aggressive pruning of strict, duplicate, redundant, and inconsistent active sets.
- Updated documentation for CAD/QE, optimization, parameter stratification, function ranges, intrinsic integration, exactness guarantees, limitations, architecture, and public APIs.

### Fixed

- Avoided unnecessary full three-dimensional CAD for certified monotone odd-degree implicit fiber bounds; exact algebraic root-function reconstruction now handles cases such as `z**3 + x*z + y <= 0` over a base cell with `x >= 0`.
- Full specialization of low-degree algebraic root functions now recovers a certified native radical when possible; for example, the positive root of `y**3 = sqrt(2)` presents exactly as `2**(1/6)`.
- Consolidated duplicate low-degree radical solving in the exact root isolator without changing fallback order.
- Removed the unused planner-history implementation and routed planner feature signatures directly through the live strategy-memory path.
- Reused shared relation construction, conjunction flattening, and parameter symbol resolution instead of retaining parallel private implementations.
- Finished consolidating integer-solver equality splitting and integer univariate-root extraction.
- Incorrect quadratic radical branch reconstruction when the leading coefficient changes sign by CAD base cell.
- Certification paths that could validate a symbolic algebraic-root presentation inconsistent with its ordered-root identity.
- Overlapping explicit `Or` decompositions that could double-count integrals.
- Incomparable symbolic bounds being silently dropped by cylindrical fast paths.
- Invalid or empty explicit cells constructed from contradictory bounds.
- Incorrect regularity classification for uncertified algebraic sections.
- String-variable assumption mismatches across decision, sampling, optimization, integration, reasoning, and root-classification APIs.
- Disjunctive optimization overstating global certification when a feasible branch was uncertified.
- Fixed-precision fallbacks in exact root ordering, sign determination, isolating intervals, sector sampling, interval simplification, and optimizer candidate comparison.
- CAD boundary reconstruction for closed regions such as the unit disk.
- Substitution semantics of `root_of(...)` by making its fiber variable binder-aware.

- Package-level `equivalent()` shortcuts conflating distinct same-name SymPy symbols or consuming generator-valued variable lists.
- Shared bound normalization silently accepting undeclared variables, malformed or duplicate bounds, and reversed intervals that could yield negative measure.
- Standard-region integration creating assumption-distinct symbols from string variable names.
- Floating-point ordering of exact endpoints in standard-region intersections, lazy CAD root stacks, and algebraic-root RUR fallback evaluation.
- One-dimensional component merging shrinking a union when a later interval was contained in an earlier interval.
- Remaining approximate interval-order validation in CAD invariants and implicit-geometry cell ordering.
- Parametric-region integration creating assumption-distinct ambient symbols from string names.
- `ParametricRegion` accepting undeclared, missing, or duplicate parameter limits and non-positive multiplicity.
- Standard-region constructors accepting provably reversed intervals, negative radii, reversed spherical-shell radii, or inconsistent ambient dimensions.
- Added direct behavioral coverage for every exported public function, with a regression guard that flags newly exported functions lacking a test call.
- Added interval/Boolean-region metamorphic tests for union, intersection, and difference measure identities.

### Compatibility

Version 0.2.0 contains intentional API and internal naming cleanups. Code that depended on obsolete internal names should be updated to the current public API documented in `docs/api_overview.md`.

### Internal refactoring

- Consolidated exact-algebraic and CAD process-local caches on one shared `BoundedLRU` implementation.
- Consolidated integer-solver recoverable exceptions and expression-complexity scoring in `solve/integer/_common.py`.
- Split decision input preparation, metadata collection, output selection/diagnostics, and witness validation into private modules while preserving the public decision API.
- Cached solution metadata and made structural CAD metadata demand-driven for explicit lightweight outputs; default structured solutions retain supported cells/cylindrical metadata, while connectivity is computed only when requested or needed for component sampling.
- Distinguished cheap relation-to-zero-RHS conversion from canonical polynomial relation normalization with explicit names and compatibility aliases.
- Added shared variable-normalization policies for problem variables, sampling variables, and explicit symbol sequences.
- Evaluated splitting `optimization.py`; deferred it because the current range/image and optimization internals are mutually coupled, and a safe split requires a larger dependency-boundary redesign rather than a mechanical file move.
