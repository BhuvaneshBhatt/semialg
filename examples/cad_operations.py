"""Closure, interior, and boundary operations with cad."""

import sympy as sp

from semialg import cad

x = sp.Symbol("x", real=True)
region = (x >= 0) & (x <= 1)

print("region:", cad(region, [x]))
print("interior:", cad(region, [x], operation="interior"))
print("closure:", cad(region, [x], operation="closure"))
print("boundary:", cad(region, [x], operation="boundary"))
