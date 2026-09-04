# Parametric integration with algebraic endpoints

Integrate over x²≤a while CAD controls the parameter-dependent root branches.

## Problem and interpretation

The fiber

\[
\{x:x^2\le a\}
\]

is empty for \(a<0\), collapses to a point at \(a=0\), and has algebraic
endpoints \(\pm\sqrt a\) for \(a>0\).

The implementation constructs a CAD in parameter-plus-fiber space so root
number and ordering are stable on each parameter cell. It then integrates
between exact root-function endpoints. The coalescence point \(a=0\) is part
of the exact stratification rather than a numerically special-cased value.

## Executable example

```python
import sympy as sp

from semialg import integrate_over_region

x, a = sp.symbols("x a", real=True)

result = integrate_over_region(
    1,
    x**2 <= a,
    [x],
    parameters=[a],
    return_stratified=True,
)

print(result)
assert result.certified is True
assert result.method == "parametric_algebraic_root_cad_integration"
assert result.select({a: -1}) == 0
assert result.select({a: 0}) == 0
assert result.select({a: 4}) == 4
assert result.select({a: 9}) == 6
```

The matching executable file is `examples/gallery/11_parametric_algebraic_endpoints.py`; the documentation
test suite runs every gallery script.
