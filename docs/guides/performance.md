# Performance guide

Exact real algebraic computation can be expensive. CAD in particular has severe worst-case complexity, so performance depends strongly on formulation, variable order, degree growth, and how much of a decomposition must actually be built.

## Start with the default policy

The default planners use structural scoring, estimated lifting/root counts, algebraic-degree and coefficient-height information, bounded pilot lifting, equational constraints, partial lifting, and solve-scoped caches. Start with defaults unless you have evidence that a particular order or certification policy is poor.

## Reduce the problem before CAD

Cheap symbolic simplification can have a large effect:

- remove redundant constraints;
- expose equality constraints explicitly;
- eliminate variables with simple equalities when safe;
- factor polynomials when that clarifies Boolean structure;
- avoid introducing unnecessary auxiliary variables.

Optimization does some equality-based dimension reduction automatically.

## Variable order matters

Different CAD variable orders can change projection degree, coefficient height, root-isolation burden, and cell count dramatically; on difficult inputs the runtime and memory difference can be orders of magnitude. Variable order should therefore be one of the first things inspected when a CAD-backed computation is unexpectedly expensive.

`semialg`'s automatic planner combines structural incidence/degree scores with projection statistics, estimated lifting/root counts, and bounded pilot lifting for leading candidates. CAD-driven region integration also considers alternative coordinate orders because a good order can turn algebraic/root-function bounds into simple polynomial bounds. Expert code can inspect order suggestions through `semialg.heuristics.suggest_variable_order` and `suggest_cad_variable_order`.

If you supply an order manually, benchmark it on the actual problem family rather than assuming a syntactic heuristic will generalize. Preserve logical structure: variables in one quantifier block may be reorderable, but moving variables across alternating `Exists`/`ForAll` blocks is not semantics-preserving. For reconstructed free-variable cells, also consider whether a particular cylindrical orientation is useful to the downstream task.

## Equational constraints and partial CAD

A formula containing an equation such as `f == 0` can often be decomposed more cheaply than a fully sign-invariant CAD for every polynomial. semialg propagates equational constraints across levels and can avoid lifting cells that cannot affect the formula's truth.

Existential satisfiability now uses exhaustive lazy CAD before materializing a full
cell decomposition. It stops immediately on a certified true branch, while UNSAT is
returned only after all relevant branches have been exhausted. Inside a single
quantifier block the lazy driver may reorder variables with the measured CAD order
planner because such reordering preserves semantics.

Collins projection also works with exact squarefree factors rather than projecting a
reducible product and all of its factors simultaneously. Factor signs plus pairwise
resultants retain the required sign/root events and reduce redundant coefficients and
discriminants. Reduced McCallum/Lazard paths retain their conservative certification
and Collins fallback contract.

Write logical structure explicitly rather than hiding useful equalities inside opaque transformations.

## Reuse within a solve

A top-level exact operation creates a computation context that reuses projection data, sign determinations, root comparisons, specializations, and RUR computations across nested work. Keep logically related work inside the high-level operation rather than manually recreating equivalent low-level calls when possible.

## Optimization certification policies

`semialgebraic_minimize` and `semialgebraic_maximize` support:

- `certification="auto"`: use cost estimates to decide whether expensive range CAD is justified;
- `certification="candidate"`: avoid the full range fallback while retaining cheaper candidate certification paths;
- `certification="complete"`: permit the complete exact range computation even when estimated expensive.

`range_cost_limit` controls the automatic policy. Raising it can make a computation much more expensive; it is not simply a precision setting.

## Separate exploration from proof

For plotting or exploratory sampling, explicitly numerical workflows may be much cheaper. Do not force a complete CAD merely to obtain points for visualization. Conversely, do not treat exploratory numerical samples as proof.

## Diagnose before optimizing code

When a problem is slow, record:

1. variable count and order;
2. polynomial degrees;
3. equality constraints;
4. projection/lifting estimates or diagnostics;
5. whether the cost comes from candidate generation, root isolation, lifting, or global certification.

A different mathematical formulation often matters more than micro-optimizing Python code.

## Performance expectations

Low-dimensional polynomial problems are the intended sweet spot. Higher-dimensional or high-degree problems can become expensive even when the final answer is simple. See [Limitations](../limitations.md) and [Errors and failure modes](errors_and_failure_modes.md).

## Structural cache keys

Projection and exact-algebraic caches use immutable SymPy/`Poly` structure for identity rather than serializing expressions with `sstr` or `srepr`. This avoids repeated expression-to-string conversion in hot cache and deduplication paths and preserves exact symbol identity, including assumption-distinct symbols that share a printed name. Human-readable strings are still used for diagnostics and provenance, but not as mathematical identity keys.

Collins and reduced CAD projection also share the same low-level polynomial normalization and projection operations. This reduces repeated `Poly -> Expr -> Poly` conversion while keeping the projection algorithms themselves separate and auditable.

## Tunable resource and heuristic limits

Several resource limits are explicit
options so applications can trade runtime and memory for aggressiveness without
patching semialg internals:

- `cad(..., max_preprocess_aux_vars=3)` bounds existential auxiliaries introduced
  by semialgebraic preprocessing; `None` disables this guard.
- `semialgebraic_minimize(..., max_boolean_branches=32)` and
  `semialgebraic_maximize(...)` bound DNF branch expansion before per-branch
  optimization.
- `semialg.simplify.simplify_semialgebraic_formula` exposes
  `max_dnf_branches`, `max_implication_atoms`, `max_implication_vars`, and
  `max_implication_degree`. These affect how aggressively semantic redundancy is
  searched for, not the meaning of the returned formula.
- `semialg.planner.candidate_variable_orders` exposes the returned candidate
  `limit`, `exhaustive_var_limit`, `projection_var_limit`, and
  `projection_shortlist`.
- `sample_point`/`sample_points` expose `default_sampling_radius`,
  `random_attempts_min`, `random_attempts_each`, `max_random_denominator`,
  and `numeric_precision`; an explicit `random_attempts` still overrides the
  automatic attempt policy.
- `sign_at` and `sign_vector` expose `numeric_precision` for their explicitly
  inexact (`exact=False`) fallback.

Process-local performance-cache capacities are expert controls rather than root
APIs. Use `semialg.algebraic.configure_algebraic_cache_limits(...)` for root,
sign, comparison, specialization, and RUR caches, and
`semialg.cad_algorithms.configure_cad_cache_limits(...)` for projection-tower,
squarefree-base, projection-step, and complete-CAD caches. Resizing a cache keeps
its newest entries up to the new capacity and does not change mathematical
semantics.

## Cache lifecycle and teardown

Long-running workers and large test processes can release every process-local
`semialg` performance cache through `clear_caches()`. The call clears the
algebraic root/sign/RUR caches, CAD projection caches, GTZ canonical-basis
caches, optimization/planner/convexity/integration LRUs, integer Groebner
recursion cache, simplification implication cache, and solution-metadata cache.
It runs Python garbage collection by default. Set `include_sympy=True` only at
an explicit process lifecycle boundary when it is also appropriate to clear
SymPy's global expression cache.

`cache_report()` returns cache sizes, configured bounds, and available hit/miss
statistics without exposing cached mathematical objects. Cache state is a
performance detail only: clearing it must not change exact results or
certificate validity. The cache lifecycle tests exercise repeated workloads,
checks every bounded cache stays within its declared capacity, clears all
caches through the unified hook, and then replays the same certificates.

## Performance regression testing

Prefer structural assertions for routine regression tests. A test that proves a Sturm fallback was not entered, that a projection cache gained a hit, or that a specialized backend avoided generic QE is more portable than a narrow elapsed-time bound. Use wall-clock tests only for coarse end-to-end regressions and keep them in the opt-in performance suite.

If an accelerated backend bypasses a fallback cache, test the accelerated backend and the fallback cache as separate contracts. Do not require cache counters to move on a path that no longer uses that cache.
