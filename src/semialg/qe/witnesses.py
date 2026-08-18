from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass

import sympy as sp

from ..algebraic.samples import sample_to_expr
from ..cad.lifting.stack import CADCell
from ..formula import Formula
from .complete import CompleteQEResult, evaluate_formula_on_cell


@dataclass(frozen=True)
class Witness:
    assignment: Mapping[sp.Symbol, sp.Expr]
    source_cell: CADCell | None = None


def witness_from_qe_result(
    qe_result: CompleteQEResult, free_point: Mapping[sp.Symbol, sp.Expr] | None = None
) -> Witness | None:
    if qe_result.is_sentence:
        if qe_result.truth_value:
            return Witness(assignment=dict(free_point or {}), source_cell=None)
        return None
    if qe_result.cell_union is None:
        return None
    free_point = dict(free_point or {})
    for cell in qe_result.cell_union.cells:
        ok = True
        for var, _sample in zip(qe_result.free_variables[: cell.level], cell.sample, strict=True):
            if var not in free_point:
                continue
            # Exact containment tests are delegated to the reconstructed formula
            # when possible; this helper only provides a convenient sample.
            if not sp.simplify(qe_result.formula.subs(free_point)):
                ok = False
                break
        if ok:
            assignment = dict(free_point)
            for var, sample in zip(
                qe_result.free_variables[: cell.level], cell.sample, strict=True
            ):
                assignment.setdefault(var, sample_to_expr(sample))
            return Witness(assignment=assignment, source_cell=cell)
    return None


def existential_wit_assign(
    formula: Formula,
    vars_: Sequence[sp.Symbol],
    free_assignment: Mapping[sp.Symbol, sp.Expr],
    candidate_cells: Sequence[CADCell],
) -> Witness | None:
    for cell in candidate_cells:
        point = dict(free_assignment)
        point.update(
            {
                var: sample_to_expr(sample)
                for var, sample in zip(vars_[: cell.level], cell.sample, strict=True)
            }
        )
        if evaluate_formula_on_cell(formula, cell, vars_):
            return Witness(assignment=point, source_cell=cell)
    return None


__all__ = ["Witness", "witness_from_qe_result", "existential_wit_assign"]
