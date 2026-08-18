"""FindInstance-style examples across supported domains."""

import sympy as sp

from semialg import find_instance

x, y, p, q = sp.symbols("x y p q", real=True)

print("Real algebraic instances:")
print(find_instance(sp.Eq(x**2, 2), [x], count=2).instances)

print("Complex instances:")
print(find_instance(sp.Eq(x**2 + 1, 0), [x], domain="complexes", count=2).instances)

print("Integer instances:")
print(find_instance(sp.Eq(x**2, 4), [x], domain="integers", count=2).instances)

print("Boolean instances:")
print(find_instance(p | q, [p, q], domain="booleans", count=3).instances)
