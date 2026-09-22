from __future__ import annotations

import sympy as sp

from ..model import QEResult
from ..simplify.boolean import simplify_boolean
from .cell_formulas import qe_cells_and_vars
from .closure import cells_closure_formula
from .interior import cells_interior_formula


def cells_boundary_formula(cells_with_truth, variables):
    closure = cells_closure_formula(cells_with_truth, variables)
    interior = cells_interior_formula(cells_with_truth, variables)
    return simplify_boolean(sp.And(closure, sp.Not(interior)))


def qe_boundary(qe_result: QEResult):
    cells, variables = qe_cells_and_vars(qe_result)
    return cells_boundary_formula(cells, variables)
