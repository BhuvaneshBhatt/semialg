"""Executable companion to the semialg documentation example."""

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
