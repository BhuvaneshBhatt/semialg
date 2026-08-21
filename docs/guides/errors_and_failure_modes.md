# Errors and failure modes

Several very different outcomes can look like “the solver did not return the answer I expected.” Distinguishing them is important.

## Invalid input

Examples include malformed bounds, bounds on undeclared variables, provably reversed intervals, negative geometric radii, inconsistent ambient dimensions, duplicate/missing parametric limits, and ambiguous same-name symbols.

These are input-contract errors and should fail early, usually with `ValueError`, rather than being interpreted as empty geometry.

## Infeasible or empty problem

A valid formula may simply have no real solution. That is a mathematical result, not an algorithm failure. Use structured decision/solution results when you need to distinguish this state programmatically.

## Unbounded objective or range

An optimization problem can be feasible but unbounded. This is different from an unsupported computation. Likewise, an infimum can be finite but unattained on an open set.

Check structured optimization fields such as `value`, `attained`, and `certified` rather than inferring status from a single expression.

## Unsupported exact case

Some valid semialgebraic problems lie outside a specialized fast path or exact representation currently implemented by semialg. Certified code should decline such a step rather than silently substitute a floating-point decision.

A different backend, a simpler formulation, or complete CAD may still solve the problem.

## Resource or cost limit

`certification="auto"` may decline an expensive range-CAD fallback based on its cost model. That does not mean the mathematical claim is false. It means the requested global certificate was not attempted under the current cost policy.

Use `certification="complete"` only when you intentionally accept the potentially much larger computation.

## Candidate found, certification incomplete

Optimization may find exact KKT/active-set candidates before proving global optimality. A candidate value and a global certificate are separate pieces of information. Inspect `OptimizationResult.certified`.

## Symbolic comparison undecidable by the current method

Exact root and bound logic can encounter expressions whose order cannot be established by the available exact comparator. The correct behavior is conservative failure or fallback to another exact representation—not a fixed-precision guess.

## Numerical mode

Some sampling/plotting paths allow explicitly numerical operation. These results are appropriate for exploration but should not be interpreted as exact certificates.

## Debugging checklist

When reporting a failure, include:

- the exact SymPy formula and variable objects;
- semialg and SymPy versions;
- the API and options used;
- whether strings or `Symbol` objects were supplied;
- the structured result/diagnostics if available;
- whether the problem succeeds under a different certification policy or variable order.

See also [Exactness and certification](../concepts/exactness_and_certification.md) and [Performance](performance.md).
