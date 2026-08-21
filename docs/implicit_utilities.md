# Implicit-region utilities

This page documents utilities inspired by symbolic region-processing workflows.

## `semialgebraic_level_function`

For Boolean combinations of inequalities, this returns a scalar expression `F` such that the condition with strict inequalities relaxed is equivalent to `F <= 0` in supported cases.

```python
import sympy as sp
from semialg import semialgebraic_level_function

x, y = sp.symbols("x y", real=True)

semialgebraic_level_function(sp.And(x**2 + y**2 <= 1, y >= 0), [x, y])
# Max(x**2 + y**2 - 1, -y)
```

## `decompose_implicit_formula`

This decomposes a Boolean formula into disjunctive pieces. Each piece records inequality polynomials and equality polynomials.

```python
from semialg import decompose_implicit_formula

pieces = decompose_implicit_formula(sp.Or(x < 0, sp.And(x >= 1, sp.Eq(y, 0))), [x, y])

for piece in pieces:
    print(piece.inequalities, piece.equalities)
```

Each piece corresponds to:

```text
f1 <= 0 and ... and fm <= 0 and g1 == 0 and ... and gn == 0
```

after inequality expansion, strict-to-weak relaxation, and supported normalization.

## `extract_symbolic_box_bounds`

This detects explicit box constraints.

```python
from semialg import extract_symbolic_box_bounds

extract_symbolic_box_bounds(sp.And(x >= 0, x <= 1, y >= -2, y <= 3), [x, y])
# bounds for x in [0, 1], y in [-2, 3]
```

## `decompose_cylindrical_formula_to_vertical_bounds_2d`

This parses supported CAD-like cylindrical formulas and ordinary vertical-slice Boolean formulas into 2D vertical bounds.

```python
from semialg import decompose_cylindrical_formula_to_vertical_bounds_2d

cells = decompose_cylindrical_formula_to_vertical_bounds_2d(
    sp.And(x >= 0, x <= 1, y >= x, y <= 1),
    [x, y],
)
```

The output records `x` intervals or points and associated lower/upper bounds for `y`. This is used by `reduce_region_integral`, `integrate_over_region`, and `region_boundary` in supported cases.

## Limitation

The 2D syntactic parser is intentionally limited. For complete CAD output, use `extract_cylindrical_solution`, which provides typed arbitrary-dimensional bounds and certified algebraic root functions.
