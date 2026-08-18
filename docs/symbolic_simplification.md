# Symbolic simplification

The simplification layer uses CAD/QE-backed semantic checks to remove contradictions, unreachable branches, and provably redundant conditions.

## Boolean simplification

```python
import sympy as sp
from semialg import simplify_boole

x, y = sp.symbols("x y", real=True)

simplify_boole(sp.And(x > 0, x >= 0), [x])
# x > 0

simplify_boole(sp.Or(x < 0, x >= 0), [x])
# True

simplify_boole(sp.And(x**2 < 0, y > 0), [x, y])
# False

# Relational atom canonicalization
simplify_boole(2*x - 2 >= 0, [x])
# x >= 1

simplify_boole(-x <= -1, [x])
# x >= 1

# Univariate formulas can be converted to interval normal form.
simplify_boole(sp.And(x**2 <= 1, x >= 0), [x])
# (x >= 0) & (x <= 1)

simplify_boole(sp.Or(x**2 > 1, sp.And(x >= -1, x <= 1)), [x])
# True

result = simplify_boole(sp.And(x > 1, x >= 0), [x], return_result=True)
result.formula
# x > 1
result.method
# 'semantic_boolean_simplification'
```

## System simplification

```python
from semialg import simplify_system

simplify_system([x > 0, x >= 0], [x])
# x > 0

simplify_system([x**2 < 0, y > 0], [x, y])
# False
```

## Piecewise simplification

```python
from semialg import simplify_piecewise

expr = sp.Piecewise((1, x**2 < 0), (2, x > 0), (3, True))
simplify_piecewise(expr, [x])
# Piecewise((2, x > 0), (3, True))
```

## Simplification under assumptions

```python
from semialg import simplify_under_assumptions

simplify_under_assumptions(sp.Abs(x), x >= 0, [x])
# x

simplify_under_assumptions(sp.sqrt((x - 1)**2), x >= 1, [x])
# x - 1

simplify_under_assumptions(sp.Max(x, y), x >= y, [x, y])
# x

simplify_under_assumptions(sp.Min(x, y), x >= y, [x, y])
# y

simplify_under_assumptions(sp.sqrt(x**2*y**2), (x >= 0) & (y <= 0), [x, y])
# -x*y

simplify_under_assumptions(sp.log(sp.exp(x)), True, [x])
# x

simplify_under_assumptions(sp.log(x**2), x > 0, [x])
# 2*log(x)

simplify_under_assumptions((x**2 - 1)/(x - 1), x > 1, [x])
# x + 1
```

Rational cancellation is domain-sensitive. By default, cancellation is applied only when the original denominator is provably nonzero under the active assumptions. If you want a conditional rewrite, request side conditions:

```python
result = simplify_under_assumptions((x**2 - 1)/(x - 1), True, [x], return_conditions=True)
result.expression
# x + 1
result.conditions
# (Ne(x - 1, 0),)
```

## Assumption-based expression simplification

The assumption simplifier handles `simplify_under_assumptions` with shifted and product square-root rewrites, real `log(exp(x))`, domain-aware `log(x**2)`, safe rational cancellation, and `AssumptionSimplificationResult` side-condition reporting.

## Scope

These simplifiers are semantic but not fully canonical. `simplify_boole` performs conservative relation canonicalization, exact univariate interval simplification, and implication-based redundancy removal, but it still does not attempt to compute a unique shortest formula for every semialgebraic set. `simplify_under_assumptions` intentionally avoids silent domain-changing rewrites unless the needed condition is already implied by the assumptions.

## Piecewise and system simplification

The semantic simplification layer includes the semantic simplification layer beyond Boolean formulas.

### `simplify_piecewise`

`simplify_piecewise(expr, variables, assumptions=..., return_result=True)` returns a `PiecewiseSimplificationResult` when requested. The simplifier:

- removes branches whose effective condition is unsatisfiable after earlier branches are taken into account;
- records skipped branches in `removed_unreachable`;
- simplifies each branch value under the assumptions plus the branch condition;
- merges adjacent semantically equal branches when possible;
- collapses to a single expression when assumptions force one branch to cover all feasible inputs.

Example:

```python
expr = sp.Piecewise((sp.sqrt(x**2), x >= 0), (-x, True), evaluate=False)
simplify_piecewise(expr, [x])
# Piecewise((x, x >= 0), (-x, True))

simplify_piecewise(expr, [x], assumptions=x >= 0)
# x
```

### `simplify_system`

`simplify_system` has a stronger system-cleanup pass before semantic redundancy removal. It detects simple acyclic equalities such as `y == x + 1`, substitutes them into the remaining constraints, and records the substitutions in `SimplifiedSystem.substitutions`. Passing `eliminate_equalities=True` removes those defining equalities from the returned constraint system when doing so is safe for the requested output.

The function also accepts `output="constraints"` and `output="result"` in addition to the default formula output.

```python
result = simplify_system([sp.Eq(y, x + 1), y >= 2], [x, y], return_result=True)
result.substitutions
# {y: x + 1}

simplify_system([sp.Eq(y, x + 1), y >= 2], [x, y], eliminate_equalities=True)
# x >= 1, or an equivalent one-dimensional condition
```

## Structured sign proofs

The sign-proving helpers preserve their historical Boolean return values by
 default while also supporting `return_result=True` for diagnostics and
 counterexamples:

```python
result = prove_positive(x**2, [x], return_result=True)
assert not result.proven
assert result.counterexample == {x: 0}
```

The structured `SignProofResult` records the requested relation, normalized
assumptions, variables, proof/refutation method, an optional certificate, and an
optional validated counterexample. Cheap syntactic certificates are tried before
falling back to the general implication/decision layer: constant signs, obvious
sums of squares/even-power products, negative sums of squares, and positive
constant-plus-squares forms.
