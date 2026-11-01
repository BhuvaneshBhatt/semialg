# Function range computation

Use `function_range` to find which real outputs are reachable under polynomial or supported semialgebraic constraints. This is useful for inverse feasibility, reachable-output analysis and exact validation of an optimizer's proposed value. A range can be disconnected, so the primary answer is an output formula.

<!-- examples: executable -->

## Exact image and interval bounds

```python
import sympy as sp
from semialg import function_range
from semialg.decision import equivalent

x, t, z = sp.symbols("x t z", real=True)
image = function_range(1 / x, sp.Ne(x, 0), [x], value_symbol=t, return_result=True)
assert equivalent(image.formula, sp.Ne(t, 0), (t,))
assert image.is_interval is False and image.interval_count == 2
assert image.diagnostics["range_kind"] == "exact_image"
```

The exact modes `qe` (the default), `cad` and `auto` describe attainable values. They use exact image shortcuts and graph projection, and raise `NotImplementedError` when the required exact projection is unsupported. They do not fill image holes with an enclosing interval.

```python
source = sp.Eq((x - 1) * (x - 3), 0)
image = function_range(x, source, [x, z], value_symbol=t, return_result=True)
assert image.formula.as_set() == sp.FiniteSet(1, 3)
bounds = function_range(
    x,
    source,
    [x, z],
    value_symbol=t,
    method="bounds",
    return_result=True,
)
assert bounds.formula.as_set() == sp.Interval(1, 3)
assert bounds.diagnostics["range_kind"] == "interval_enclosure"
```

`bounds` and its alias `optimization` request the interval between the extremal bounds, with endpoint attainment information. That interval can contain values absent from the image. The result's interval metadata describe the returned formula. See [exact images and bounds](range-images.md) for the distinction.

## Graph projection and source definedness

The image of `f` over a source condition `C` is the existential relation `exists x: C(x) and graph_f(x,t)`. Exact graph construction preserves the source's definition guards. A rational equation `t = p/q` uses both `t*q - p = 0` and `q != 0`. Inverse-map substitutions must preserve their own poles as well.

```python
positive = function_range(1 / x, x > 0, [x], value_symbol=t)
assert equivalent(positive, t > 0, (t,))
restricted = function_range(x, sp.Or(x <= -1, x >= 1), [x], value_symbol=t)
assert equivalent(restricted, sp.Or(t <= -1, t >= 1), (t,))
```

If SymPy already erased a source operation before the query, the package cannot recover it. Supply its excluded points or other original source conditions explicitly. A principal square root and a real odd root are different graphs; use `sympy.real_root` when a real odd root is intended.

## Supported graphs and exact reductions

| Source form | Representation | Conditions retained |
| --- | --- | --- |
| Polynomial or rational expression | Direct polynomial graph | Nonzero denominators |
| `Abs`, `sign`, `Heaviside`, finite `Min`/`Max` | Auxiliary graph relations | Supported nested arguments and zero-value choices |
| Rational powers and canonical `real_root` forms | Exact algebraic graphs | Real locus and principal or selected real branch |
| Finite `Piecewise` | Ordered branch graph | Branch ownership and supported semialgebraic conditions |
| Selected commensurate trigonometric families | Circle variables with `s**2 + c**2 = 1` | An exact reduction that removes all remaining source dependence |
| Selected commensurate exponential/hyperbolic families | `u = exp(g*x)` with `u > 0` | Exact transformed source constraints |
| General mixed transcendental expressions | No general semialgebraic graph | A supported exact reduction is required |

```python
semicircle = function_range(sp.sqrt(1 - x * x), True, [x], value_symbol=t)
assert equivalent(semicircle, sp.And(t >= 0, t <= 1), (t,))
absolute = function_range(sp.Abs(x), True, [x], value_symbol=t)
assert equivalent(absolute, t >= 0, (t,))
```

Commensurate reductions are exact only when the transformed objective and constraints no longer depend on the original source variable. They do not replace unrelated transcendental values by independent bounded coordinates. Incommensurate phases and expressions such as `x + sin(x)` need other analytic arguments. Restricted source intervals can also require phase/domain information that a full circle or positive-axis substitution would lose. See [transcendental solving scope](guides/transcendental_scope.md).

## Metadata

`return_result=False` returns the mathematical condition. `return_result=True` returns a `FunctionRangeResult` with `formula`/`range_condition`, extremal summaries, endpoint attainment, interval structure, method and diagnostics.

```python
result = function_range(
    2 * x + 1,
    sp.And(x > 0, x < 1),
    [x],
    value_symbol=t,
    return_result=True,
)
assert result.range_condition.as_set() == sp.Interval.open(1, 3)
assert result.lower_bound == 1 and result.upper_bound == 3
assert not result.lower_bound_attained and not result.upper_bound_attained
```

The lower and upper summaries alone cannot recover a disconnected image. Consumers should use the formula and `diagnostics["range_kind"]` to determine whether they received an exact image or an explicitly requested interval enclosure.

## Declaring parameters

Symbol roles must be explicit. Ordinary variable normalization can add free symbols from the objective and constraints to the source-variable list. Passing `variables=[x]` alone does not declare every other symbol a free parameter. For a parameter-dependent image, pass `parameters=[...]` with `return_stratified=True`.

```python
a = sp.Symbol("a", real=True)
family = function_range(
    x + a,
    sp.And(x >= 0, x <= 1),
    [x],
    value_symbol=t,
    parameters=[a],
    return_stratified=True,
    eliminate_quantifiers=True,
)
for av in (-2, 0, 3):
    specialized = sp.Or(
        *(sp.And(branch.condition, branch.value.formula).subs(a, av) for branch in family.branches)
    )
    assert equivalent(specialized, sp.And(t >= av, t <= av + 1), (t,))
```

The default `eliminate_quantifiers=False` retains an exact first-order relation in each `ParametricFunctionRangeResult`, along with its normalized `quantifiers` prefix. Its `formula` is the matrix to which that prefix applies. Use `apply_quantifiers` to build a first-class quantified expression, or opt into elimination as above. The extra QE cost is separate from exactness, and a compact quantifier-free result is not promised for every family. Candidate degree-loss strata and source-domain exceptions must remain in the parameter conditions.

## Analytic ranges and optimization

For general supported transcendental image analysis, use [funcprops](https://github.com/BhuvaneshBhatt/funcprops/blob/main/docs/range-analysis.md), which supplies explicit implicit-critical-value, periodic, tail and partial-image interfaces. For optimizer selection, infima, minimizing points and optimization certificate replay, use [symbopt](https://github.com/BhuvaneshBhatt/symbopt/blob/main/docs/guide/ranges-and-extrema.md). semialg supplies exact real algebraic projection and side conditions; it does not require either consumer package.

The image, its closure, its extremal bounds and its attained extrema are different mathematical objects. More decimal precision does not replace a proof of source definedness, exhaustive root coverage or endpoint membership. See [package ownership](package_ownership.md), [optimization](optimization.md) and the [range reference](reference/optimization_and_range.md).
