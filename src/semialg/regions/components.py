from __future__ import annotations

import sympy as sp

from ..model import QEResult
from .common import cell_to_formula, qe_cells_and_vars


def region_components(cells_with_truth, variables):
    pieces = []
    for cell, truth in cells_with_truth:
        if not truth:
            continue
        piece = cell_to_formula(cell, variables, left_closed=True, right_closed=True)
        if piece is not sp.false:
            pieces.append(sp.simplify(piece))
    components = []
    seen = set()
    for piece in pieces:
        key = sp.srepr(piece)
        if key not in seen:
            seen.add(key)
            components.append(piece)
    return components


def qe_components(qe_result: QEResult):
    cells, variables = qe_cells_and_vars(qe_result)
    return region_components(cells, variables)
