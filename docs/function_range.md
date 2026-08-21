# Function range computation

`function_range` computes the real range of an expression over a semialgebraic domain.

## Range as a semialgebraic image

The range of `f(x)` over a domain `C(x)` is represented by an existential
formula.  Programmatically, semialg represents it with `Exists`:

```python
from semialg import Exists

image_relation = Exists(x, sp.And(C, graph_f))
```

Mathematically this is $\exists x\,[C(x) \land \operatorname{graph}_f(x,t)]$,
where `t` is a value symbol. Eliminating the original variables gives a condition on `t`.

```python
import sympy as sp
from semialg import function_range

x, t = sp.symbols("x t", real=True)

function_range(x, sp.Or(x <= -1, x >= 1), [x], value_symbol=t)
# Abs(t) >= 1
```

The result may be disconnected, so the primary answer is a formula, not just a pair of bounds.

## Polynomial and rational expressions

For rational expressions, `semialg` clears denominators:

```text
t = p(x)/q(x)
```

becomes:

```text
t*q(x) - p(x) == 0 and q(x) != 0
```

Example:

```python
function_range(1/x, x > 0, [x], value_symbol=t)
# t > 0
```

## Semialgebraic expression graphs

The current implementation supports common semialgebraic expressions by introducing graph constraints.

```python
function_range(sp.Abs(x), True, [x], value_symbol=t)
# t >= 0

function_range(sp.sqrt(1 - x**2), True, [x], value_symbol=t)
# (t >= 0) & (t <= 1)

function_range(sp.Max(x, 0), sp.And(x >= -1, x <= 2), [x], value_symbol=t)
# (t >= 0) & (t <= 2)

function_range(sp.Min(x, 1), sp.And(x >= 0, x <= 3), [x], value_symbol=t)
# (t >= 0) & (t <= 1)
```

Simple `Piecewise` expressions are also supported when their branches can be represented semialgebraically.

## Metadata

Use `return_result=True` for metadata.

```python
r = function_range(
    2*x + 1,
    sp.And(x > 0, x < 1),
    [x],
    value_symbol=t,
    return_result=True,
)

r.range_condition
r.lower_bound
r.upper_bound
r.lower_bound_attained
r.upper_bound_attained
r.is_interval
r.interval_count
```

The primary answer is `range_condition`; `lower_bound` and `upper_bound` are summaries and may lose information for disconnected ranges.

## Parameter-stratified ranges

When an expression or its domain depends on symbolic parameters, pass `parameters=[...]` and `return_stratified=True`. The result is a `ParameterStratifiedResult`; each branch contains an exact `ParametricFunctionRangeResult` guarded by a semialgebraic parameter condition.

```python
a = sp.Symbol("a", real=True)

r = function_range(
    x + a,
    sp.And(x >= 0, x <= 1),
    [x],
    parameters=[a],
    return_stratified=True,
)
```

`ParametricFunctionRangeResult.formula` is an exact first-order image relation.
For new expression-facing code, wrap that relation with `Exists` (or use
`apply_quantifiers`) when a single first-class quantified expression is desired;
`quantifiers` remains available as the normalized internal elimination prefix. The relation is intentionally not forced through a second potentially expensive CAD elimination: `quantifier_free` is therefore `False`. This preserves an exact first-class parametric answer without making a range query unexpectedly perform a much larger QE problem.

## Limitations

`function_range` currently targets polynomial, rational, and common semialgebraic expression graphs. Arbitrary transcendental expressions are outside the exact real-closed-field setting and are not generally supported.
