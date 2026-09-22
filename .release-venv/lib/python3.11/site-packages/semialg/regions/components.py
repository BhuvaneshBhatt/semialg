from __future__ import annotations

import sympy as sp

from ..model import QEResult
from .cell_formulas import cell_to_formula, qe_cells_and_vars


def components_from_cells(cells_with_truth, variables):
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
        if piece not in seen:
            seen.add(piece)
            components.append(piece)
    return components


def qe_components(qe_result: QEResult):
    cells, variables = qe_cells_and_vars(qe_result)
    return components_from_cells(cells, variables)
