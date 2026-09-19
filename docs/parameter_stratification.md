# Parameter-stratified and conditional results

Many exact symbolic answers depend on parameter regions. `semialg` represents these answers explicitly rather than attaching an unconditional expression to a representative sample.

## Public API

- `ConditionalBranch`
- `ParameterStratifiedResult`
- `ParameterStratificationCertificate`
- `conditional_result`
- `verify_parameter_stratification`
- `solvability_conditions(..., return_stratified=True)`
- `root_count_conditions(..., return_stratified=True)`
- `ParameterizedCylindricalDecomposition.as_stratified_result()`
- `semialgebraic_minimize(..., parameters=[...], return_stratified=True)`
- `semialgebraic_maximize(..., parameters=[...], return_stratified=True)`
- `function_range(..., parameters=[...], return_stratified=True)`
- `integrate_over_region(..., parameters=[...], return_stratified=True)`
- `semialgebraic_measure(..., parameters=[...], return_stratified=True)`
- `ParametricOptimizationResult`
- `ParametricFunctionRangeResult`

A branch consists of an exact semialgebraic guard and a value valid under that guard.

```python
import sympy as sp
from semialg.conditional import conditional_result, verify_parameter_stratification

a = sp.Symbol("a")
result = conditional_result(
    [a],
    [(a < 0, -a), (a >= 0, a)],
)

result.select({"a": -3})
# 3

result.as_piecewise()
# Piecewise((-a, a < 0), (a, True))

verify_parameter_stratification(result).verify()
# True
```

String assignment keys are resolved to the actual parameter symbols carried by the result, so symbol assumptions are preserved.

## Certification model

Three claims are checked independently:

1. every branch value is certified on its guard;
2. guards are pairwise disjoint;
3. guards cover the requested parameter domain.

`verify_parameter_stratification` returns a `ParameterStratificationCertificate` containing these checks, overlap conditions, and any uncovered condition.

A sampled CAD fiber is evidence for a parameter cell, not automatically a symbolic formula valid on the whole cell. `ParameterizedCylindricalDecomposition.as_stratified_result()` therefore exposes guarded `ParameterStratum` objects rather than promoting representative fibers to unconditional answers.

## Optimization and range relations

Parameterized optimization and range calls return guarded exact first-order result objects. They do not assume that a representative fiber is constant over its parameter cell. `ParametricOptimizationResult.formula` together with `.quantifiers` characterizes the infimum/supremum relation, including nonattainment; `ParametricFunctionRangeResult.formula` plus `.quantifiers` characterizes image membership. Both expose `quantifier_free=False` to make the representation explicit. This design keeps stratified calls predictable: requesting a first-class parameter result does not implicitly force a second full CAD elimination.


## Optional relation elimination

Parametric optimization and function-range results retain their exact quantified first-order definitions by default. Set `eliminate_quantifiers=True` to perform the additional complete-CAD elimination explicitly. This can be substantially more expensive, so it is not implied by `return_stratified=True`.

## Parametric integration

For semialgebraic regions whose inequalities contain symbolic parameters, region integration stratifies parameter space, reduces the region with symbolic parameter-dependent bounds, and attaches exact integral expressions to the certified guards. Empty fibers are represented by a zero-valued branch, giving complete parameter-space coverage. This is different from `ParametricRegion`, which represents an explicit geometric parametrization.

## Representative parameter samples

A representative sample is auxiliary data, not a proof that an arbitrary convenient parameter value belongs to a stratum. When CAD supplies a cell sample, semialg validates it against that cell. When a finer parameter decomposition is unavailable, the fallback searches for and validates an exact representative of the full parameter condition. It never assumes that zero is feasible. If no representative can be certified for a nonempty stratum, the operation raises `NotImplementedError` rather than specializing the problem at an invalid point.

Sampled fibers are constructed through the same internal path for CAD cells and fallback strata, so specialization and validation semantics remain consistent.
