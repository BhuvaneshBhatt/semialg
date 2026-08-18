import sympy as sp

from semialg import find_instance

x = sp.Symbol("x", real=True)
result = find_instance(sp.Eq(x**2, 2), [x], count=2)
print("exact:", result.instances)
print("approx:", result.approximate)
