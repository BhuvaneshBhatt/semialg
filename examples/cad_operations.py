"""Exact closure, interior and boundary of a closed interval."""

import sympy as sp

from semialg import cad
from semialg.decision import equivalent

x = sp.Symbol("x", real=True)
region = sp.And(x >= 0, x <= 1)
interior = cad(region, [x], operation="interior")
closure = cad(region, [x], operation="closure")
boundary = cad(region, [x], operation="boundary")
assert equivalent(interior, sp.And(x > 0, x < 1), (x,))
assert equivalent(closure, region, (x,))
assert equivalent(boundary, sp.Or(sp.Eq(x, 0), sp.Eq(x, 1)), (x,))
print(interior, closure, boundary)
