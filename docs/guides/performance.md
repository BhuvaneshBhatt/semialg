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

Different CAD variable orders can change projection degree and cell count dramatically. semialg's automatic planner combines projection statistics with estimated lifting/root counts and bounded pilot lifting for leading candidates.

If you supply an order manually, benchmark it on the actual problem family rather than assuming a syntactic heuristic will generalize.

## Equational constraints and partial CAD

A formula containing an equation such as `f == 0` can often be decomposed more cheaply than a fully sign-invariant CAD for every polynomial. semialg propagates equational constraints across levels and can avoid lifting cells that cannot affect the formula's truth.

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
