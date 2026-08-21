# Getting started

This guide introduces the main public workflows without requiring knowledge of CAD internals.

## 1. Create exact variables and formulas

semialg operates on SymPy expressions and relations.

```python
import sympy as sp

x, y = sp.symbols("x y", real=True)
region = (x**2 + y**2 <= 1) & (x >= 0)
```

Passing the actual `Symbol` objects is the clearest style. String names are also accepted by many APIs and are resolved against symbols already present in the problem. See [Symbol handling](guides/symbol_handling.md).

## 2. Ask decision questions

```python
from semialg import equivalent, implies, is_satisfiable

is_satisfiable(region & (y > 0), [x, y])
# True

implies(x >= 2, x**2 >= 4, [x])
# True

equivalent(x**2 <= 1, (x >= -1) & (x <= 1), [x])
# True
```

Use `return_result=True` when you need structured evidence such as a witness or counterexample.

```python
res = implies(x >= 0, x > 0, [x], return_result=True)
res.valid
# False
res.counterexample
# typically {x: 0}
```

## 3. Solve a semialgebraic system

```python
from semialg import solve_semialgebraic

sol = solve_semialgebraic(
    sp.Eq(x**2 + y**2, 1) & sp.Eq(x, y) & (x > 0),
    [x, y],
    count=4,
)
sol.samples
# contains the exact point {x: sqrt(2)/2, y: sqrt(2)/2}
```

`solve_semialgebraic` can return samples, formulas, cells, and structured solution metadata depending on the requested output. See [Solving and sampling](reference/solving_and_sampling.md).

## 4. Optimize exactly

```python
from semialg import semialgebraic_minimize

opt = semialgebraic_minimize(
    x**2 + y**2,
    x + y >= 1,
    [x, y],
)
opt.value
# 1/2
opt.attained
# True
opt.certified
# True for an exactly certified global result
```

Open sets can have an exact infimum that is not attained. semialg records this separately rather than conflating “minimum” with “infimum”.

## 5. Compute a function range

```python
from semialg import function_range

t = sp.Symbol("t", real=True)
function_range(x, sp.Or(x <= -1, x >= 1), [x], value_symbol=t)
# Abs(t) >= 1   (or an equivalent formula)
```

Range computation is an image/QE problem and may be more expensive than finding candidate extrema.

## 6. Work with regions

```python
from semialg import region_intersection, region_union

left = (x >= -1) & (x <= 0)
right = (x >= 0) & (x <= 1)

region_union(left, right)
region_intersection(left, right)
```

For common geometric objects, use standard-region classes such as `IntervalRegion`, `BoxRegion`, `BallRegion`, and `ParametricRegion`. Their constructor invariants are summarized in [Region invariants](guides/region_invariants.md).

## 7. Integrate and measure

```python
from semialg import integrate_over_region, semialgebraic_measure

semialgebraic_measure(x**2 + y**2 <= 1, [x, y])
# pi

integrate_over_region(x**2 + y**2, x**2 + y**2 <= 1, [x, y])
# pi/2
```

Intrinsic measure is available for supported lower-dimensional regular strata:

```python
semialgebraic_measure(
    sp.Eq(x**2 + y**2, 1),
    [x, y],
    measure_dimension="intrinsic",
)
# 2*pi
```

## 8. Parameters and guarded answers

Some answers change qualitatively with parameters. semialg can represent these as guarded branches rather than extrapolating from a representative sample.

```python
from semialg import classify_real_roots

a = sp.Symbol("a", real=True)
classification = classify_real_roots(x**2 - a, x, parameters=[a])
```

See [Parameters and conditional results](reference/parameters.md).

## 9. Understand the status of an answer

Before building downstream logic around a result, ask:

- Is the representation exact?
- Was the mathematical conclusion certified?
- Is this only a candidate or heuristic?
- Was numerical approximation explicitly requested?

These distinctions are explained in [Exactness and certification](concepts/exactness_and_certification.md).

## 10. When a problem is slow

CAD complexity can grow rapidly with variable count and degree. Before forcing a complete computation, see the [Performance guide](guides/performance.md) and [Errors and failure modes](guides/errors_and_failure_modes.md).
