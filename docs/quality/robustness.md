# Robustness and exactness

`semialg` treats exactness and certification as explicit API contracts. Exact
algorithms either produce a result justified by symbolic or algebraic reasoning
or decline the computation; finite-precision numerics are not used as proof.

## Symbol identity

Public APIs that accept variable or parameter names resolve those names against
the symbols already present in the input expressions. This preserves SymPy
assumptions and object identity. If distinct symbols share the same printed
name, callers must pass the desired `Symbol` object explicitly.

## Exact algebraic decisions

The certified algebraic layer uses exact operations for:

- real-root isolation and deduplication;
- rational isolating intervals and interval refinement;
- algebraic sign and order comparisons;
- CAD sector-sample construction;
- root-order and sample-containment certificates;
- one-dimensional endpoint ordering and component merging;
- virtual-substitution witness ordering;
- lazy-CAD root ordering and algebraic-root RUR fallback ordering.

Certified paths do not fall back to `nroots()` or `evalf()` when an exact
operation is unavailable. Low-degree fibers with exact algebraic coefficients
are solved symbolically when possible; otherwise the operation is declined.


## Bounds and explicit-region validation

Shared bound normalization rejects undeclared variables, duplicate entries, malformed
bound tuples, and exactly reversed endpoints. Public measure/integration APIs therefore
do not silently ignore misspelled bounds or turn a reversed interval into a negative
measure.

One-dimensional component reconstruction compares algebraic endpoints exactly and
merges overlapping intervals by their true union, including containment and endpoint
closure.

## Root-function semantics

`root_of(p, y, k)` treats `y` as a bound fiber variable. Substitution may
specialize base parameters without replacing the bound variable inside the
defining polynomial. Once the base parameters are specialized, the selector can
resolve to the requested exact real root.

`AlgebraicRootFunction` keeps CAD root identity, base-cell validity, ordering,
and regularity certification separate from optional radical presentation.

## CAD and topology certificates

Reduced CAD paths are accepted only when their side conditions and invariants
are certified. Cell-bound verification uses exact comparisons, and topology
incidence specializes algebraic root functions before evaluating relations.

Decomposition certification distinguishes three claims:

1. individual cell bounds are valid;
2. cells are pairwise disjoint where required;
3. the cells cover the represented solution set.

These claims are not inferred from one another.

## Optimization certification

Exact optimization compares algebraic candidate values exactly. A disjunctive
optimum is globally certified only when every feasible branch has been
certified. Cost-controlled CAD range certification may be declined in automatic
mode; `certification="complete"` explicitly requests the full exact fallback.

## Computation context

A top-level solve uses an `ExactComputationContext` to share transient caches
across nested CAD, QE, optimization, range, root, sign, specialization, and RUR
operations. Bounded process-level caches provide reuse across independent calls
without allowing a single solve to retain unbounded state.

## Conservative failure

Expected strategy failure, unavailable exact arithmetic, and failed
certification should be represented explicitly. Internal programming errors
must not be converted into successful fallback results. Some heuristic
subsystems still contain broad exception handlers; these should be narrowed as
their failure contracts are made explicit.

## Validation

The test suite includes exactness, symbol-identity, root-ordering, topology,
optimization-certificate, parameter-stratification, and CAD consistency tests.
Computationally expensive CAD cases are marked `slow` so CI can run them with an
appropriate time budget.

## Region and API validation

The standard-region layer enforces exact, local construction invariants so malformed
regions cannot silently reach integration. The parametric-region layer preserves
ambient symbol identity and validates complete one-to-one parameter limits and
positive multiplicity. Deterministic metamorphic tests verify interval measure
identities for union, intersection, and difference.

A public-function coverage test scans the test suite and compares direct function
calls with the exported function set. This is a release guard against adding public
functions without at least one behavioral contract test.


## Algebraic root presentation

Low-degree root presentation is interval-certified. When SymPy cannot build a
native `RootOf` over an algebraic coefficient domain, semialg keeps the ordered
root as an exact isolating-interval object. A readable radical is returned only
when exactly one certified real radical candidate lies in that interval. This
prevents presentation logic from changing CAD root identity.
