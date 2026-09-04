# Tutorial: exact parameterized computation

A parameter can change not only the numerical value of an answer but also the topology of the feasible set, the number/order of algebraic roots, and whether an optimum is attained. `semialg` represents those changes with exact semialgebraic guards.

```python
import sympy as sp

x = sp.Symbol("x", real=True)
a = sp.Symbol("a", real=True)
```

## 1. A parameterized interval

Start with

\[
0\le x\le a.
\]

For \(a<0\), the region is empty. For \(a\ge0\), it is the interval \([0,a]\).

### Function range

```python
from semialg import function_range

r = function_range(
    x,
    (x >= 0) & (x <= a),
    [x],
    parameters=[a],
    return_stratified=True,
)
```

The exact answer must distinguish the empty and nonempty parameter regimes. If a quantifier-free relation is required, request `eliminate_quantifiers=True`; supported affine one-dimensional families first use specialized direct reconstruction before falling back to complete QE.

### Optimization

```python
from semialg import semialgebraic_minimize, semialgebraic_maximize

mn = semialgebraic_minimize(
    x,
    (x >= 0) & (x <= a),
    [x],
    parameters=[a],
    return_stratified=True,
)

mx = semialgebraic_maximize(
    x,
    (x >= 0) & (x <= a),
    [x],
    parameters=[a],
    return_stratified=True,
)
```

For \(a\ge0\), the minimum is \(0\) and maximum is \(a\). The parameter guard is part of the exact answer.

### Measure

```python
from semialg import semialgebraic_measure

m = semialgebraic_measure(
    (x >= 0) & (x <= a),
    [x],
    parameters=[a],
    return_stratified=True,
)
```

The empty regime has measure zero; the nonempty regime has length \(a\).

### Integration

```python
from semialg import integrate_over_region

I = integrate_over_region(
    x**2,
    (x >= 0) & (x <= a),
    [x],
    parameters=[a],
    return_stratified=True,
)
```

On \(a\ge0\),

\[
\int_0^a x^2\,dx=\frac{a^3}{3}.
\]

## 2. Algebraic endpoints: \(x^2\le a\)

Consider

\[
S_a=\{x:x^2\le a\}.
\]

This is qualitatively different:

- \(a<0\): empty;
- \(a=0\): one point;
- \(a>0\): interval \([-\sqrt a,\sqrt a]\).

The critical parameter value \(a=0\) is where two algebraic boundary roots coalesce.

```python
length = integrate_over_region(
    1,
    x**2 <= a,
    [x],
    parameters=[a],
    return_stratified=True,
)

length.select({a: -1})
# 0

length.select({a: 4})
# 4

length.select({a: 9})
# 6
```

For the positive parameter cell, the integration bounds are algebraic root functions. The parameter CAD guarantees that their number and ordering are stable on the cell before the antiderivative is evaluated at those exact endpoints.

## 3. Why parameter strata matter

Sampling a single parameter value is not a proof that the same formula remains valid everywhere nearby. Changes occur at discriminants, root collisions, feasibility boundaries, and other projection conditions.

A stratified result therefore records a collection of exact guards and branch values rather than extrapolating from representative fibers.

## 4. Quantified relation versus `Piecewise`

For range and optimization problems, the cheapest exact answer may be a guarded first-order relation with explicit quantifiers. Setting `return_stratified=True` does not automatically request another expensive QE solely for presentation.

Use `eliminate_quantifiers=True` when a quantifier-free relation is required. The implementation tries specialized exact reconstruction first for supported families and falls back to complete QE when necessary.

## 5. Selecting and verifying branches

Where supported:

```python
result.select({a: 4})
```

selects the exact branch containing that parameter assignment.

`verify_parameter_stratification(...)` can check coverage and pairwise disjointness of guards for generic stratified results.

## 6. Current scope

Parameterized one-variable integration with algebraic CAD/root-function endpoints is supported for the implemented reducible cases. Arbitrary multidimensional parameter-dependent algebraic integration is outside the general supported capability.

For the representation model, see [Understanding result objects](../concepts/result_objects.md). For the lower-level parameter API, see [Parameter stratification](../parameter_stratification.md).
