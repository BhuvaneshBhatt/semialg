# CAD and quantifier elimination

Many `semialg` operations reduce questions to first-order formulas over the real numbers. Examples include satisfiability, implication, range/image computation, optimization certificates, subset checks, and parameter conditions.

## Exactness

CAD/QE is used as a correctness engine for polynomial equations and inequalities over the reals. Fast paths may avoid CAD, but exact public paths should not replace an unresolved algebraic sign/order with a fixed-precision numerical guess.


## First-class quantified expressions

For programmatic formulas, prefer semialg's `Exists` and `ForAll` nodes rather
than manually building quantifier-prefix tuples.

```python
import sympy as sp
from semialg import Exists, ForAll
from semialg.solve import reduce_complete_expr

x, y = sp.symbols("x y", real=True)

formula = ForAll(x, Exists(y, sp.Eq(x + y, 0)))
reduce_complete_expr(formula)
# True
```

Text interfaces such as `reduce_complete_text("forall x. ...")` remain
supported.  The solver lowers both forms to the same internal prenex blocks.

## Variable ordering

CAD is highly sensitive to variable ordering. The current planner combines inexpensive Brown/SOTD/NDRR-style scores with more expensive scoring on a small candidate shortlist. For small problems it considers the actual projection tower and estimates lifting cost using real fiber-root counts and an estimated cell count, not projection size alone. The shortlist score also records projected coefficient-height growth and an algebraic-degree proxy. The two best projection-scored orders receive a bounded pilot lift over a few representative stacks; the measured result refines rather than replaces the cheaper estimate.

## Caching

Repeated exact work is memoized at two levels. `ExactComputationContext` is a transient per-operation cache automatically shared by nested CAD/QE/optimization calls, so repeated projection, root, sign, comparison, specialization, and RUR queries within one solve are reused without polluting global state. Existing bounded process-local LRUs remain a second-level cache. CAD caches cover projection towers, projection steps, and squarefree bases; algebraic caches cover root isolation, exact signs, sample/root comparisons, algebraic-root specialization/evaluation, and RUR construction. Use `computation_context()` explicitly when several related public calls should share one solve-local context.

## Equational constraints and partial CAD

Lazy CAD performs prefix truth evaluation and short-circuits subtrees whose truth value is already determined. Conjunctively necessary equations can produce lower-level necessary resultants; these derived ECs may prune prefixes before later variables are lifted. For existential variables, a necessary EC can justify section-only lifting. Universal levels remain conservative and are not restricted to equality sections.

The implementation intentionally uses only logically necessary derived constraints. A projection polynomial or resultant is not treated as an EC merely because it appears algebraically.

## Reduced projection and fallback

Reduced/EC-aware projection paths carry side-condition/certification information. If the reduced path cannot establish the required invariance, the solver falls back to the complete Collins-style path rather than assuming the reduced decomposition is complete.

## Practical implication

A common execution pattern is:

```text
normalize formula
  -> specialized solver/virtual substitution/RUR if applicable
  -> choose CAD order using structural + arithmetic + estimated lifting cost
  -> bounded pilot lift of the top orders
  -> context-scoped + bounded cached projection
  -> lazy/EC-aware lifting and prefix pruning
  -> certification
  -> conservative fallback to complete CAD when needed
```
