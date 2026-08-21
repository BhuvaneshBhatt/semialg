from __future__ import annotations

import sympy as sp

from .exact_arithmetic import compare_exact_reals
from .model import CADResult


class CADInvariantError(ValueError):
    pass


def validate_cad_result(cad: CADResult) -> list[str]:
    issues: list[str] = []
    for level, cells in cad.cells_by_level.items():
        for cell in cells:
            if cell.level != level:
                issues.append(
                    f"Cell {cell.index} stored at level {level} but has level {cell.level}"
                )
            if len(cell.sample) != cell.level:
                issues.append(f"Cell {cell.index} has mismatched sample length")
            if len(cell.intervals) != cell.level:
                issues.append(f"Cell {cell.index} has mismatched interval length")
            for left, right in cell.intervals:
                if left is not None and right is not None:
                    if sp.simplify(left - right) != 0:
                        try:
                            if compare_exact_reals(left, right) >= 0:
                                issues.append(
                                    f"Cell {cell.index} has non-increasing interval "
                                    f"({left}, {right})"
                                )
                        except (TypeError, ValueError, NotImplementedError):
                            pass
    for parent_index, children in cad.children_by_parent.items():
        for child in children:
            if child.parent_index != parent_index:
                issues.append(f"Child {child.index} does not point back to parent {parent_index}")
    return issues
