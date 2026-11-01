# Exact images and interval bounds

<!-- examples: executable -->

Use `function_range` to determine which output values are reachable. The exact modes `qe`, `cad` and `auto` preserve missing values and disconnected components. They raise `NotImplementedError` if the available exact graph projection cannot finish. They do not replace a failed image computation with its enclosing interval.

```python
import sympy as sp
from semialg import function_range

x, y, z = sp.symbols("x y z", real=True)
image = function_range(1 / x, sp.Ne(x, 0), [x], value_symbol=y, return_result=True)
assert image.formula.as_set() == sp.S.Reals - sp.FiniteSet(0)
assert image.interval_count == 2 and image.is_interval is False
assert image.diagnostics["range_kind"] == "exact_image"

condition = sp.Eq((x - 1) * (x - 3), 0)
image = function_range(x, condition, [x, z], value_symbol=y, method="qe", return_result=True)
assert image.formula.as_set() == sp.FiniteSet(1, 3)
```

`method="bounds"` and its alias `"optimization"` explicitly request the interval between the infimum and supremum. Its endpoints and attainment data come from optimization, but interior values can be missing from the actual image. With `return_result=True`, `diagnostics["range_kind"]` is `"interval_enclosure"`. Interval metadata describe the returned formula.

```python
bounds = function_range(x, condition, [x, z], value_symbol=y, method="bounds", return_result=True)
assert bounds.formula.as_set() == sp.Interval(1, 3)
assert bounds.diagnostics["range_kind"] == "interval_enclosure"
```

The rational inverse shortcut keeps poles of the inverse map before simplifying substituted constraints. For `x=1/y`, the guard `y != 0` is separate from the substituted source condition `1/y != 0`. Unsupported non-rational inverse branches and parameter cases that can change the map's degree are left to other graph or parameter-stratification routes. The shortcut does not prove that an inverse expression is real by omitting an imaginary unit.

## Parameters and consumer packages

Ordinary variable normalization can add undeclared free symbols from the objective and source constraints to the source coordinates. Declare coefficient parameters through `parameters=[...]` and `return_stratified=True`; the [function-range guide](function_range.md#declaring-parameters) gives an executable example. A quantified parameter relation can be exact even when it has no compact quantifier-free presentation.

For supported general transcendental images, use [funcprops](https://github.com/BhuvaneshBhatt/funcprops/blob/main/docs/range-analysis.md). For minimizing points and optimization certificates, use [symbopt](https://github.com/BhuvaneshBhatt/symbopt/blob/main/docs/guide/ranges-and-extrema.md). A partial image, an interval enclosure, an infimum and an attained minimum have separate contracts. semialg does not treat a consumer's sampled transcendental values as an exact algebraic graph.
