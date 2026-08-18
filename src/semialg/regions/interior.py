from __future__ import annotations

from ..model import QEResult
from .common import qe_cells_and_vars, region_formula


def region_interior(cells_with_truth, variables):
    return region_formula(cells_with_truth, variables, left_closed=False, right_closed=False)


def qe_interior(qe_result: QEResult):
    cells, variables = qe_cells_and_vars(qe_result)
    return region_interior(cells, variables)
