# Decision and QE reference
## Family contract

**computer algebra systeml return.** Decision functions return exact truth values or structured exact decision results for first-order formulas over the reals. QE operations return a quantifier-free formula equivalent over the declared real variables.

**Exactness and certification.** Successful certified paths do not use fixed-precision numerical sign guesses. Specialized methods may decline and allow another exact backend to run.

**Algorithm selection.** The planner may use normalization/presolve, quadratic virtual substitution for supported low-degree fragments, zero-dimensional algebraic solving where applicable, and complete or reduced CAD/QE.

**Complexity and limitations.** Complete real QE has severe worst-case complexity. A specialized backend accepting a formula is a performance choice, not a weaker semantics. See [How semialg chooses an algorithm](../concepts/algorithm_selection.md).




## Primary API overview

This table is the substantive coverage target for the primary APIs assigned to this reference page. Each entry states the API's primary role; the family contract and detailed sections below explain shared algorithms, exactness guarantees, and limitations. It is maintained together with `docs/reference/primary_api_manifest.toml`, and documentation tests require every root-level primary API to map here rather than merely appearing in the generated public index.

| API | Kind | Role / return |
|---|---|---|
| `apply_quantifiers` | function | Wrap ``formula`` in a prenex quantifier prefix. |
| `computation_context` | function | Create/reuse a context for one complete exact solve operation. |
| `equivalent` | function | Return whether two semialgebraic formulas define the same real set. |
| `Exists` | class | Existentially quantify one or more variables in a Boolean formula. |
| `ForAll` | class | Universally quantify one or more variables in a Boolean formula. |
| `function_domain` | function | Return exact recognized real-domain constraints for supported ``expr`` forms. |
| `implies` | function | Return whether ``assumptions`` imply ``conclusion`` over the reals. |
| `is_real_valued` | function | Return whether supported domain conditions follow from assumptions. |
| `is_satisfiable` | function | Return whether a real semialgebraic formula has a satisfying point. |
| `is_tautology` | function | Return whether a real semialgebraic formula is true for all variables. |
| `prove_negative` | function | Return whether the expression is certified negative on the stated domain. |
| `prove_nonnegative` | function | Return whether the expression is certified nonnegative on the stated domain. |
| `prove_nonpositive` | function | Return whether the expression is certified nonpositive on the stated domain. |
| `prove_positive` | function | Return whether the expression is certified positive on the stated domain. |
| `reduce_formula` | function | Reduce a parsed real formula using the selected exact decision strategy. |
| `resolve_formula` | function | Resolve a parsed formula and return its exact solution representation. |
| `SatisfiabilityResult` | class | Structured result for a real satisfiability query. |
| `SemialgebraicSolution` | class | Structured solution summary for a semialgebraic constraint system. |
| `SemialgOptions` | class | Shared options accepted by high-level semialgebraic APIs. |
| `simplify_boole` | function | Simplify a semialgebraic Boolean formula over the real numbers. |
| `simplify_piecewise` | function | Simplify a Piecewise expression using semialgebraic branch conditions. |
| `simplify_system` | function | Simplify a real semialgebraic system with CAD/QE-backed checks. |
| `simplify_under_assumptions` | function | Simplify real expressions using provable assumptions. |
| `solve_semialgebraic` | function | Reduce, sample, and summarize a semialgebraic system over the reals. |
| `SolveDomain` | class | Enumeration of supported high-level solving domains. |

## `is_satisfiable(formula, variables=None, *, domain="reals", strategy=None, return_result=False)`

Decides whether a real assignment satisfies `formula`. With `return_result=False`, returns a Boolean. With `return_result=True`, returns a `SatisfiabilityResult` containing structured status and a validated witness when available.

## `is_tautology(...)`

Decides whether a formula holds for all assignments of the declared variables.

## `implies(assumptions, conclusion, variables=None, *, ..., return_result=False)`

Checks whether the assumptions imply the conclusion. A structured `ImplicationResult` can contain a validated counterexample when implication fails.

## `equivalent(lhs, rhs, variables=None, *, ..., return_result=False)`

Checks logical equivalence over the declared real variables. Printed-expression equality is not used as a substitute for symbolic identity or proof.

## `qe_by_complete_cad(...)`

Runs the complete-CAD quantifier-elimination backend for supported real polynomial formulas. It returns the eliminated Boolean formula by default. Set `return_result=True` when CAD cells, witnesses, truth-propagation state, or diagnostics are needed. The companion expert APIs `qe_from_cad`, `qe_prenex`, `qe_prenex_suffix`, `qe_blocks`, `qe_parsed`, and `qe_text` use the same convention.

## Variables and domains

Variables may be SymPy symbols or, in supported APIs, unambiguous string names. See [Symbol handling](../guides/symbol_handling.md).

## Result semantics

A Boolean decision is a mathematical conclusion. A witness/counterexample returned by structured decision APIs is validated against the original formula before exposure.

For implementation and strategy details, see [Decision procedures](../decision_procedures.md) and [Exactness and certification](../concepts/exactness_and_certification.md).
