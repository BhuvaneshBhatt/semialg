# Parameters and conditional results reference
## Family contract

**Mathematical return.** Parameter APIs return guarded exact branches, parameter stratifications, and certificates describing where each symbolic result is valid.

**Exactness and certification.** A representative sample from a parameter cell is not promoted to an unconditional symbolic answer. Guards, coverage, and disjointness are explicit parts of the result model.

**Algorithm.** Parameter space may be decomposed by CAD/projection conditions; operation-specific code then attaches exact optimization, range, integration, root-count, or solvability results to those guards.

**Complexity and limitations.** Forcing a quantifier-free presentation may require an additional complete QE and can be substantially more expensive than retaining an exact quantified relation.

## Why stratified results exist

A parameterized problem can have qualitatively different answers on different parameter regions. semialg represents this explicitly rather than evaluating one representative parameter point and pretending its answer is globally valid.

## `ParameterStratifiedResult`

Represents guarded branches together with partition/certification metadata. Branch conditions can be checked for disjointness and coverage.

## `conditional_result`

Constructs a guarded conditional result from explicit cases.

## `verify_parameter_stratification`

Checks the structural/logical validity of a proposed parameter partition in supported settings.

## Root conditions

- `solvability_conditions`
- `root_count_conditions`
- `classify_real_roots`

These expose parameter-dependent existence/count/classification information.

`SolvabilityConditionsResult` deliberately has Boolean truth only for the two
unconditional cases: an identically true solvability condition is truthy and an
identically false condition is falsey. A genuinely parameter-dependent condition
raises `TypeError` in `bool(result)`; inspect `result.formula`,
`result.is_conditional`, `result.is_unconditionally_solvable`, or
`result.is_never_solvable` instead. This prevents a condition such as `a >= 0`
from being mistaken for unconditional solvability.

## Parametric optimization and ranges

`semialgebraic_minimize`, `semialgebraic_maximize`, and `function_range` accept `parameters=[...]` with `return_stratified=True`.

Their branch values can be exact first-order relations with explicit quantifier prefixes. `quantifier_free=False` means the exact relation has not been subjected to an additional QE merely to simplify its presentation; it does **not** mean the relation is numerical or approximate.

## Symbol identity

Parameter names supplied as strings are resolved against the original symbols in the problem. Ambiguous same-name symbols are rejected. See [Symbol handling](../guides/symbol_handling.md).

See [Parameter stratification](../parameter_stratification.md) for detailed examples.


`integrate_over_region` and `semialgebraic_measure` also accept `parameters=[...]` with `return_stratified=True` for parameter-dependent formula regions. Empty fibers are represented by an exact zero-valued branch, and feasible strata carry the symbolic integral/measure expression.

For parametric optimization and function ranges, `eliminate_quantifiers=True` explicitly requests the additional complete-CAD QE pass; quantified exact relations remain the default.
