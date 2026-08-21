# Parameters and conditional results reference

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

## Parametric optimization and ranges

`semialgebraic_minimize`, `semialgebraic_maximize`, and `function_range` accept `parameters=[...]` with `return_stratified=True`.

Their branch values can be exact first-order relations with explicit quantifier prefixes. `quantifier_free=False` means the exact relation has not been subjected to an additional QE merely to simplify its presentation; it does **not** mean the relation is numerical or approximate.

## Symbol identity

Parameter names supplied as strings are resolved against the original symbols in the problem. Ambiguous same-name symbols are rejected. See [Symbol handling](../guides/symbol_handling.md).

See [Parameter stratification](../parameter_stratification.md) for detailed examples.
