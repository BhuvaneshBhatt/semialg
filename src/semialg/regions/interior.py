from __future__ import annotations

from ..model import QEResult
from .cell_formulas import qe_cells_and_vars, region_formula


def cells_interior_formula(cells_with_truth, variables):
    return region_formula(cells_with_truth, variables, left_closed=False, right_closed=False)


def qe_interior(qe_result: QEResult):
    cells, variables = qe_cells_and_vars(qe_result)
    return cells_interior_formula(cells, variables)
