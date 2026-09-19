# CAD-driven integration variable ordering

Let the integration planner choose a cylindrical order that avoids unnecessary radicals and yields simpler exact iterated bounds.

## Executable example

```python
import sympy as sp

from semialg import integrate_over_region, reduce_region_integral

x, y = sp.symbols("x y", real=True)
condition = (y >= 0) & (y <= 1) & (x >= y**2) & (x <= y)

reduced = reduce_region_integral(1, condition, [x, y])
value = integrate_over_region(1, condition, [x, y])
piece = reduced.pieces[0]

print("method:", reduced.method)
print("selected order:", piece.diagnostics["integration_variable_order"])
print("limits:", piece.limits)
print("integral:", value)

assert reduced.method == "coordinate_permuted_cylindrical_integration"
assert piece.diagnostics["integration_variable_order"] == (y, x)
assert piece.limits == ((x, y**2, y), (y, 0, 1))
assert value == sp.Rational(1, 6)
```

The matching executable file is `examples/gallery/16_cad_integration_variable_order.py`; the documentation test suite runs every gallery script.
