# Errors and failure modes

Several very different outcomes can look like “the solver did not return the answer I expected.” Distinguishing them is important.




## Primary API overview

This table is the substantive coverage target for the primary APIs assigned to this reference page. Each entry states the API's primary role; the family contract and detailed sections below explain shared algorithms, exactness guarantees, and limitations. It is maintained together with `docs/reference/primary_api_manifest.toml`, and documentation tests require every root-level primary API to map here rather than merely appearing in the generated public index.

| API | Kind | Role / return |
|---|---|---|
| `ResourceLimitError` | class | A computation stopped because a configured resource limit was reached. |
| `SemialgError` | class | Base class for semialg-specific failures. |
| `SemialgStrategyFailure` | class | A speculative exact/symbolic strategy could not handle the input. |
| `UnsupportedFragmentError` | class | The input is valid, but outside the symbolic fragment a strategy supports. |

## Invalid input

Examples include malformed bounds, bounds on undeclared variables, provably reversed intervals, negative geometric radii, inconsistent ambient dimensions, duplicate/missing parametric limits, and ambiguous same-name symbols.

These are input-contract errors and should fail early, usually with `ValueError`, rather than being interpreted as empty geometry.

## Infeasible or empty problem

A valid formula may simply have no real solution. That is a mathematical result, not an algorithm failure. Use structured decision/solution results when you need to distinguish this state programmatically.

## Unbounded objective or range

An optimization problem can be feasible but unbounded. This is different from an unsupported computation. Likewise, an infimum can be finite but unattained on an open set.

Check structured optimization fields such as `value`, `attained`, and `certified` rather than inferring status from a single expression.

## Unsupported exact case

Some valid semialgebraic problems lie outside a specialized fast path or exact representation implemented by semialg. Certified code should decline such a step rather than silently substitute a floating-point decision.

A different backend, a simpler formulation, or complete CAD may still solve the problem.

## Resource or cost limit

`certification="auto"` may decline an expensive range-CAD fallback based on its cost model. That does not mean the mathematical claim is false. It means the requested global certificate was not attempted under the current cost policy.

Use `certification="complete"` only when you accept the potentially much larger computation.

For `cad`, `max_preprocess_aux_vars` is an enforced preprocessing limit. Exceeding it
raises `ResourceLimitError` for direct or strict calls. A non-strict call with
`return_result=True` instead returns `status="unknown"`, preserves the input formula,
and records the limit in diagnostics. This formula is not a completed decomposition;
converting the result with `as_function()` raises `ResourceLimitError`.

The CAD backend does not implement the `max_cells` and `timeout` controls. Supplying
either control raises `NotImplementedError` for direct or strict calls, or returns an
unknown structured result with an explanatory diagnostic. These controls must not be
interpreted as enforced computation budgets. Cached results do not bypass this validation.

## Candidate found, certification incomplete

Optimization may find exact KKT/active-set candidates before proving global optimality. A candidate value and a global certificate are separate pieces of information. Inspect `OptimizationResult.certified`.

## Symbolic comparison undecidable by the current method

Exact root and bound logic can encounter expressions whose order cannot be established by the available exact comparator. The correct behavior is conservative failure or fallback to another exact representation—not a fixed-precision guess.

## Numerical mode

Some sampling/plotting paths allow explicitly numerical operation. These results are appropriate for exploration but should not be interpreted as exact certificates.


## Exception hierarchy

Package-specific failures derive from `semialg.errors.SemialgError`. Important subclasses distinguish unsupported fragments/strategy failure, backend failure, formula normalization, algebraic solving, quantifier elimination, reconstruction, certification, exact evaluation, dimension mismatches, and configured resource limits. `semialg.exceptions` provides the same exception classes as a dedicated exception import surface.

Fallback code should catch the narrowest expected strategy exceptions it can justify. Programming defects such as `AssertionError` should propagate rather than being converted into an apparently harmless fallback.

## Debugging checklist

When reporting a failure, include:

- the exact SymPy formula and variable objects;
- the installed semialg and SymPy versions;
- the API and options used;
- whether strings or `Symbol` objects were supplied;
- the structured result/diagnostics if available;
- whether the problem succeeds under a different certification policy or variable order.

See also [Exactness and certification](../concepts/exactness_and_certification.md) and [Performance](performance.md).

## Exception reference

See the complete [exception hierarchy](../reference/exceptions.md) for the canonical classes and guidance on which level to catch.
