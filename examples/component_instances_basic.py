"""Inspect connected-component metadata from the expert decomposition API."""

import sympy as sp

from semialg.decomposition import component_instances

x = sp.Symbol("x", real=True)
result = component_instances((x < -1) | (x > 1), [x], return_result=True)
for component in result.components:
    print(component.id, component.dimension, component.sample_exact)
