"""Use specialist algorithms from their owning namespaces."""

import sympy as sp

from semialg.map_degree import parametric_map_degree
from semialg.parameters import root_count_conditions

x, a = sp.symbols("x a", real=True)

counts = root_count_conditions(x**2 - a, x, (a,))
degree = parametric_map_degree((x**2,), (x,))

assert sp.simplify(counts[2] ^ (a > 0)) is sp.false
assert degree.degree == 2
