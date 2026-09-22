from __future__ import annotations

from ..model import QEResult
from .cell_formulas import qe_cells_and_vars, region_formula


def cells_closure_formula(cells_with_truth, variables):
    return region_formula(cells_with_truth, variables, left_closed=True, right_closed=True)


def qe_closure(qe_result: QEResult):
    cells, variables = qe_cells_and_vars(qe_result)
    return cells_closure_formula(cells, variables)
