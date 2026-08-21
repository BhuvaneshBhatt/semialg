# Decision and QE reference

## `is_satisfiable(formula, variables=None, *, domain="reals", strategy=None, return_result=False)`

Decides whether a real assignment satisfies `formula`. With `return_result=False`, returns a Boolean. With `return_result=True`, returns a `SatisfiabilityResult` containing structured status and a validated witness when available.

## `is_tautology(...)`

Decides whether a formula holds for all assignments of the declared variables.

## `implies(assumptions, conclusion, variables=None, *, ..., return_result=False)`

Checks whether the assumptions imply the conclusion. A structured `ImplicationResult` can contain a validated counterexample when implication fails.

## `equivalent(lhs, rhs, variables=None, *, ..., return_result=False)`

Checks logical equivalence over the declared real variables. Printed-expression equality is not used as a substitute for symbolic identity or proof.

## `qe_by_complete_cad(...)`

Runs the complete-CAD quantifier-elimination backend for supported real polynomial formulas. Use this lower-level API when you specifically need the QE result rather than a decision wrapper.

## Variables and domains

Variables may be SymPy symbols or, in supported APIs, unambiguous string names. See [Symbol handling](../guides/symbol_handling.md).

## Result semantics

A Boolean decision is a mathematical conclusion. A witness/counterexample returned by structured decision APIs is validated against the original formula before exposure.

For implementation and strategy details, see [Decision procedures](../decision_procedures.md) and [Exactness and certification](../concepts/exactness_and_certification.md).
