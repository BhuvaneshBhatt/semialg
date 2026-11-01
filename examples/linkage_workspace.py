"""Use rotational symmetry to derive and verify the workspace of two unit links."""

import sympy as sp

from semialg.decision import equivalent, is_satisfiable
from semialg.qe import quantifier_eliminate

x, y, u, v, r = sp.symbols("x y u v r", real=True)
radial = sp.And(r >= 0, sp.Eq(u * u + v * v, 1), sp.Eq((r - u) ** 2 + v * v, 1))
projection = quantifier_eliminate(radial, [("exists", u), ("exists", v)], variables=(r, u, v))
assert equivalent(projection, sp.And(r >= 0, r <= 2), (r,))
# Rotating any endpoint to (r,0) preserves both squared link lengths.
height = sp.sqrt(1 - r * r / 4)
assert sp.simplify((r / 2) ** 2 + height**2 - 1) == 0
assert sp.simplify((r - r / 2) ** 2 + height**2 - 1) == 0
constraints = sp.And(sp.Eq(u * u + v * v, 1), sp.Eq((x - u) ** 2 + (y - v) ** 2, 1))
assert is_satisfiable(constraints.subs({x: 1, y: 0}), (u, v))
assert not is_satisfiable(constraints.subs({x: 3, y: 0}), (u, v))
assert is_satisfiable(constraints.subs({x: 0, y: 0}), (u, v))
print("The radial projection is 0<=r<=2; rotational symmetry gives the closed disk x^2+y^2<=4.")
