"""Symbolic region conditions and active-boundary strata."""

import sympy as sp

from semialg import ParametricRegion, region_active_boundary_strata
from semialg.symbolic_regions import (
    RegionElement,
    SemialgebraicRegion,
    as_semialgebraic_region,
    region_element_conditions,
)

x, y, t = sp.symbols("x y t", real=True)
a = sp.Symbol("a")

mapped = as_semialgebraic_region(ParametricRegion((t,), ((t, 0, 1),), (t, t**2)), (x, y))
membership = RegionElement((x, y), mapped)
symbolic = region_element_conditions(membership)
eliminated = region_element_conditions(membership, eliminate=True)
assert symbolic == membership
assert eliminated != membership

parameter_region = SemialgebraicRegion(sp.And(x >= 0, x <= a), (x,))
parameter_atom = RegionElement((x,), parameter_region)
explicit_domains = region_element_conditions(parameter_atom, real_parameters=True)
assert explicit_domains.has(sp.Contains(a, sp.S.Reals, evaluate=False))

square = sp.And(x >= 0, x <= 1, y >= 0, y <= 1)
strata = region_active_boundary_strata(square, (x, y))
assert sum(stratum.active_count == 1 for stratum in strata) == 4
assert sum(stratum.active_count == 2 for stratum in strata) == 4

print("symbolic membership:", symbolic)
print("with explicit real parameters:", explicit_domains)
print("eliminated membership:", eliminated)
print("active boundary strata:", len(strata))
