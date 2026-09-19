from __future__ import annotations

import sympy as sp

from ..model import QEResult
from ..simplify.boolean import simplify_boolean
from .cell_formulas import qe_cells_and_vars
from .closure import region_closure
from .interior import region_interior


def region_boundary(cells_with_truth, variables):
    closure = region_closure(cells_with_truth, variables)
    interior = region_interior(cells_with_truth, variables)
    return simplify_boolean(sp.And(closure, sp.Not(interior)))


def qe_boundary(qe_result: QEResult):
    cells, variables = qe_cells_and_vars(qe_result)
    return region_boundary(cells, variables)
