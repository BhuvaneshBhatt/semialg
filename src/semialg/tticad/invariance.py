from __future__ import annotations

from dataclasses import dataclass

from ..formula import Formula
from ..model import CADResult
from ..qe import eval_quantifier_free


@dataclass(frozen=True)
class InvarianceCheck:
    checked_cells: int
    mismatches: tuple[tuple[int, ...], ...]

    @property
    def truth_invariant(self) -> bool:
        return not self.mismatches


def check_truth_invariance(formula: Formula, cad: CADResult) -> InvarianceCheck:
    mismatches = []
    checked = 0
    for cell in cad.cells_by_level.get(len(cad.vars), []):
        checked += 1
        expected = eval_quantifier_free(formula, dict(zip(cad.vars, cell.sample, strict=True)))
        known = cad.truth_by_cell.get(cell.index, expected)
        if bool(expected) != bool(known):
            mismatches.append(cell.index)
    return InvarianceCheck(checked_cells=checked, mismatches=tuple(mismatches))


__all__ = ["InvarianceCheck", "check_truth_invariance"]
