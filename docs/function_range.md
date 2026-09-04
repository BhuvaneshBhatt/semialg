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

When an expression or its domain depends on symbolic parameters, pass `parameters=[...]` and `return_stratified=True`. The result is a `ParameterStratifiedResult`; each branch contains an exact `ParametricFunctionRangeResult` guarded by a semialgebraic parameter condition. The branch relation is intentionally left quantified by default. Add `eliminate_quantifiers=True` when a quantifier-free formula is needed and the extra complete-CAD/QE cost is acceptable.

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
Expression-facing code, wrap that relation with `Exists` (or use
`apply_quantifiers`) when a single first-class quantified expression is desired;
`quantifiers` remains available as the normalized internal elimination prefix. The relation is intentionally not forced through a second potentially expensive CAD elimination: `quantifier_free` is therefore `False`. This preserves an exact first-class parametric answer without making a range query unexpectedly perform a much larger QE problem.

## Limitations

`function_range` targets polynomial, rational, and common semialgebraic expression graphs. Arbitrary transcendental expressions are outside the exact real-closed-field setting and are not generally supported.

## Exact algebraization of selected transcendental forms

`function_range` remains grounded in semialgebraic computation, but an input
expression need not be syntactically polynomial when it admits a certified
finite reduction to a semialgebraic image problem.  Current exact reductions
include commensurate trigonometric expressions and commensurate
exponential/hyperbolic expressions in one real source variable.

```python
import sympy as sp
from semialg import function_range

x, t = sp.symbols("x t", real=True)

function_range(sp.sin(x) + sp.cos(2*x), variables=[x], value_symbol=t)
# equivalent to -2 <= t <= 9/8

function_range(sp.exp(x) + sp.exp(-x), variables=[x], value_symbol=t)
# t >= 2

function_range(sp.cosh(2*x), variables=[x], value_symbol=t)
# t >= 1
```

For trigonometric expressions the reducer finds a common frequency `g`, writes
all `sin(n*g*x)` and `cos(n*g*x)` terms as polynomials in `s = sin(g*x)` and
`c = cos(g*x)`, and adds `s**2 + c**2 = 1`.  Because real `x` covers the whole
unit circle, this is exact when no untransformed dependence on `x` remains.

For exponential/hyperbolic expressions it writes all rates as integer multiples
of a common `g`, substitutes `u = exp(g*x)`, and adds `u > 0`.  The map from real
`x` to positive `u` is bijective, so the transformed range problem is again
exact.

The implementation does **not** replace unrelated transcendental values by
independent bounded variables.  For example, `sin(x) + sin(sqrt(2)*x)` is not
turned into an unconstrained two-circle problem, and `sin(x) + x` is not treated
as though `sin(x)` were independent of `x`.

## Supported semialgebraic graph algebra

`function_range` can eliminate an exact graph representation for rational expressions and recursively nested semialgebraic heads. The graph layer covers:

| Expression/formula form | Exact graph encoding | Important domain semantics |
|---|---|---|
| Polynomial or rational expression | Yes | Denominators are constrained nonzero. |
| `Abs(f)`, `sign(f)`, and `Heaviside(f, h0)` | Yes | Arguments, including the zero value `h0`, may themselves use supported graph forms. |
| Rational principal powers `f**(p/q)` | Yes | SymPy principal-branch real locus is preserved; negative bases are not silently reinterpreted as real odd roots. |
| `real_root(f, q)` canonical forms and integer powers | Yes | Odd-denominator real-root powers use a direct exact polynomial graph, including negative integer powers away from zero. |
| Finite `Min` / `Max` | Yes | Arguments may be nested supported graph expressions. |
| Finite `Piecewise` | Yes | Branch order is respected; supported semialgebraic function heads may occur in branch conditions. |
| Semialgebraic function heads inside constraints | Yes | Relational residuals are algebraized with auxiliary graph variables; `And`, `Or`, `Not`, `Xor`, `Implies`, and `Equivalent` compositions are supported. |
| General transcendental functions | No | Separate recognized exact reductions exist for selected commensurate trig/exponential families. |

The graph conversion introduces auxiliary variables only when needed. Cheap direct univariate image formulas are attempted before CAD/QE, so adding graph support does not force every supported range problem through the most expensive backend.
