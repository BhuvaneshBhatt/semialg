"""Local strata combine branch geometry with incident global strata."""

import sympy as sp

from semialg import local_algebraic_strata

x, y = sp.symbols("x y", real=True)
result = local_algebraic_strata((x * y,), {x: 0, y: 0}, (x, y))
assert result.complete
assert result.branch_geometry.singular
assert len(result.branch_geometry.branches) == 2
assert result.branch_geometry.intersection is not None
