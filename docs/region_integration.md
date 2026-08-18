# Region integration and measure

`semialg` provides a layered, exact-first region integration framework.

## Public API

- `reduce_region_integral`
- `integrate_over_region`
- `semialgebraic_measure`

## Architecture

The integration engine is layered:

1. recognize standard shapes where exact moment formulas are available;
2. reduce supported semialgebraic regions to iterated-integral pieces;
3. evaluate symbolically, numerically, or automatically depending on `method`;
4. support intrinsic-dimensional integration for selected lower-dimensional regions.

## Reducing to iterated integrals

```python
import sympy as sp
from semialg import reduce_region_integral

x, y = sp.symbols("x y", real=True)

red = reduce_region_integral(1, x**2 + y**2 <= 1, [x, y])
```

A reduced integral consists of one or more pieces, each with an integrand, limits, and a sign. For the unit disk, the result is equivalent to:

```python
sp.Integral(1, (y, -sp.sqrt(1 - x**2), sp.sqrt(1 - x**2)), (x, -1, 1))
```

## Symbolic, numeric, and auto modes

```python
from semialg import integrate_over_region

integrate_over_region(1, x**2 + y**2 <= 1, [x, y], method="symbolic")
# pi

integrate_over_region(1, x**2 + y**2 <= 1, [x, y], method="numeric")
# approximate numeric value

integrate_over_region(1, x**2 + y**2 <= 1, [x, y], method="auto")
# symbolic if possible, numeric otherwise
```

`method="symbolic"` is the default. It requires exact symbolic evaluation and raises `NotImplementedError` if any reduced piece remains unevaluated.

## Measure

```python
from semialg import semialgebraic_measure

semialgebraic_measure(x**2 + y**2 <= 1, [x, y])
# pi

# alias for semialgebraic_measure
```

## Intrinsic-dimensional measure

By default, measure and integration use ambient Lebesgue measure in the supplied variables.

```python
semialgebraic_measure(sp.Eq(x**2 + y**2, 1), [x, y])
# 0 under ambient 2D measure

semialgebraic_measure(sp.Eq(x**2 + y**2, 1), [x, y], measure_dimension="intrinsic")
# 2*pi in supported cases

integrate_over_region(x**2, sp.Eq(x**2 + y**2, 1), [x, y], measure_dimension=1)
# pi in supported cases
```

## Supported standard shapes and forms

The current exact layer supports many common cases, including:

- intervals;
- axis-aligned boxes;
- the standard 2D unit simplex;
- origin-centered disks and annuli;
- axis-aligned ellipses;
- selected vertical-slice regions;
- selected graph curves and circles for intrinsic one-dimensional measure.

## Limitations

Full arbitrary CAD-cell-to-bounds conversion is not yet complete. Higher-dimensional intrinsic measure is supported only for selected parametrizable cases. For unsupported cases, the package should fail conservatively rather than returning an unreliable expression.
