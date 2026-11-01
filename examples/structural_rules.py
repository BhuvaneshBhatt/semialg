"""Exact ellipsoid support, polynomial limits and Cartesian boundaries."""

import sympy as sp

from semialg import region_boundary, semialgebraic_minimize

x, y, z = sp.symbols("x y z", real=True)

# Cauchy–Schwarz in the positive-definite quadratic metric gives the bound.
closed = (x - 2) ** 2 + 4 * (y + 1) ** 2 <= 9
result = semialgebraic_minimize(3 * x + 4 * y + 5, closed, (x, y), return_result=True)
assert sp.simplify(result.value - (7 - 3 * sp.sqrt(13))) == 0
assert result.attained and result.certified
assert sp.simplify(closed.subs(result.point)) is sp.true

opened = (x - 2) ** 2 + 4 * (y + 1) ** 2 < 9
limit = semialgebraic_minimize(3 * x + 4 * y + 5, opened, (x, y), return_result=True)
assert limit.value == result.value
assert not limit.attained and limit.points == () and limit.certified

punctured = semialgebraic_minimize(x * x, sp.Ne(x, 0), (x,), return_result=True)
assert punctured.value == 0 and not punctured.attained and punctured.certified
unbounded = semialgebraic_minimize(-(x**3), sp.true, (x,), return_result=True)
assert unbounded.value == -sp.oo and not unbounded.attained and unbounded.certified

box = sp.And(x > 0, x <= 1, y >= -1, y < 2, z >= 0, z <= 3)
boundary = region_boundary(box, (x, y, z))
assert boundary.subs({x: 0, y: 0, z: 1}) is sp.true
assert boundary.subs({x: sp.Rational(1, 2), y: 0, z: 1}) is sp.false
assert boundary.subs({x: 0, y: 4, z: 1}) is sp.false
