import sympy as sp

from semialg import component_instances

x = sp.Symbol("x", real=True)
res = component_instances((x < -1) | (x > 1), [x])
for comp in res.components:
    print(comp.id, comp.dimension, comp.sample_exact)
