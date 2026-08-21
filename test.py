import sympy as sp

from semialg import equivalent, implies, is_satisfiable

x, y = sp.symbols("x y", real=True)

assert is_satisfiable(
    sp.And(x**2 + y**2 <= 1, x > 0, y > 0),
    [x, y],
)

assert implies(
    x > 1,
    x**2 > 1,
    [x],
)

assert equivalent(
    x**2 <= 1,
    sp.And(x >= -1, x <= 1),
    [x],
)

print("SemiAlg smoke test PASSED")
