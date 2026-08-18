"""Closure, interior, and boundary operations with cad."""

import sympy as sp

from semialg import cad

x = sp.Symbol("x", real=True)
region = (x >= 0) & (x <= 1)

print("region:", cad(region, [x]).formula)
print("interior:", cad(region, [x], operation="interior").formula)
print("closure:", cad(region, [x], operation="closure").formula)
print("boundary:", cad(region, [x], operation="boundary").formula)
