"""FindInstance-style examples across supported domains."""

import sympy as sp

from semialg import find_instance

x, p, q = sp.symbols("x p q", real=True)

print("Real algebraic instances:")
print(find_instance(sp.Eq(x**2, 2), [x], count=2))

print("Complex instances:")
print(find_instance(sp.Eq(x**2 + 1, 0), [x], domain="complexes", count=2))

print("Integer instances:")
print(find_instance(sp.Eq(x**2, 4), [x], domain="integers", count=2))

print("Boolean instances:")
print(find_instance(p | q, [p, q], domain="booleans", count=3))


# Sign assumptions on public real symbols are part of the witness contract.
positive_x = sp.Symbol("positive_x", positive=True)
positive_witness = find_instance(sp.true, [positive_x])
assert positive_witness is not None
assert positive_witness[positive_x] > 0
