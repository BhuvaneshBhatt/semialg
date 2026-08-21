# Optimization and range reference

## `semialgebraic_minimize`

```text
semialgebraic_minimize(
    objective, constraints=None, variables=None, *, domain="reals",
    return_result=True, certification="auto", range_cost_limit=2500,
    recursion_limit=4, parameters=None, return_stratified=False
)
```

Computes an exact minimum/infimum for supported semialgebraic problems. The pipeline can use equality reduction, stationary/KKT systems, active-set pruning, finite RUR solving, positive-dimensional critical-locus recursion, exact candidate comparison, and CAD certification.

`OptimizationResult` distinguishes:

- `value`: exact optimum/infimum value when obtained;
- `points`: known attaining optimizer points;
- `attained`: whether the value is attained;
- `certified`: whether global optimality was established exactly;
- `method`: diagnostic method information.

## `semialgebraic_maximize`

The corresponding maximum/supremum API with the same certification model.

## Certification policy

- `"candidate"`: avoid the expensive full range fallback;
- `"auto"`: use a symbolic cost model to decide whether that fallback is reasonable;
- `"complete"`: permit complete exact range certification regardless of the automatic cost estimate.

## `polynomial_locus_dimension(equations, variables)`

Returns exact dimension information used to distinguish empty, zero-dimensional, and positive-dimensional polynomial critical loci.

## `function_range`

```text
function_range(
    expression, constraints=None, variables=None, *, value_symbol=None,
    domain="reals", method="qe", return_result=False,
    parameters=None, return_stratified=False
)
```

Computes a semialgebraic image/range condition by introducing a graph relation and eliminating source variables.

## Parametric results

With `parameters=[...]` and `return_stratified=True`, optimization/range APIs can return guarded exact parameter-dependent relations rather than forcing a second QE solely to create a compact presentation.

See [Optimization](../optimization.md), [Function range](../function_range.md), and [Performance](../guides/performance.md).
