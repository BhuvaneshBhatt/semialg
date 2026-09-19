"""Executable companion to the semialg documentation example."""

import sympy as sp

from semialg import cad, find_instance, semialgebraic_minimize
from semialg.solve.zero_dimensional import solve_zero_dimensional_system

x = sp.Symbol("x", real=True)

region = cad((x >= -1) & (x <= 2), [x])
minimum = semialgebraic_minimize((x - 2) ** 2, variables=[x])
instances = find_instance(sp.Eq(x**2, 2), [x], count=2)
points = solve_zero_dimensional_system([x**2 - 1], variables=[x])

print(region)
print(minimum)
print(instances)
print(points)

assert region == ((x >= -1) & (x <= 2))
assert minimum == [0, [{x: 2}]]
assert len(instances) == 2
assert set(points) == {(-1,), (1,)}
