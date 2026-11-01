from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass

import sympy as sp

from ..algebraic.samples import sample_to_expr
from ..cad_algorithms.lifting.stack import CADCell
from ..errors import ExactComputationUnavailable
from ..exact_arithmetic import exact_truth
from ..formula import Formula, to_sympy
from .complete import CompleteQEResult


@dataclass(frozen=True)
class Witness:
    """A satisfying assignment with its source cell when containment is known."""

    assignment: Mapping[sp.Symbol, sp.Expr]
    source_cell: CADCell | None = None


def witness_from_qe_result(
    qe_result: CompleteQEResult, free_point: Mapping[sp.Symbol, sp.Expr] | None = None
) -> Witness | None:
    """Sample the free-coordinate solution of a completed QE result.

    Supplied coordinates are fixed. Missing free coordinates are completed by
    solving the specialized quantifier-free result. The source cell is located
    at the returned point, so its provenance remains valid away from CAD samples.
    This helper does not reconstruct quantified variables or Skolem functions.
    """
    if qe_result.status != "complete":
        raise ExactComputationUnavailable("cannot sample an incomplete QE result")
    point = {var: sp.sympify(value) for var, value in (free_point or {}).items()}
    free = tuple(qe_result.free_variables)
    if set(point) - set(free):
        raise ValueError("free_point contains coordinates outside the QE free variables")
    if any(value.free_symbols for value in point.values()):
        raise ValueError("free_point coordinates must be resolved real values")
    if qe_result.is_sentence:
        return Witness(point) if qe_result.truth_value else None
    union = qe_result.cell_union
    if union is None or not union.cells:
        return None
    if not point:
        cell = union.cells[0]
        return Witness(
            {var: sample_to_expr(sample) for var, sample in zip(free, cell.sample, strict=True)},
            cell,
        )
    remaining = tuple(var for var in free if var not in point)
    if remaining:
        from ..solve.find_instance import find_instance

        completion = find_instance(qe_result.formula.subs(point), remaining, strict=True)
        if completion is None:
            return None
        point.update(completion)
        if set(free) - set(point):
            raise ExactComputationUnavailable("witness completion left free coordinates unresolved")
    cell = qe_result.cad._point_locator.locate(tuple(point[var] for var in free))
    if cell is None or cell.index not in union.cell_indices:
        return None
    return Witness(point, cell)


def existential_wit_assign(
    formula: Formula,
    vars_: Sequence[sp.Symbol],
    free_assignment: Mapping[sp.Symbol, sp.Expr],
    candidate_cells: Sequence[CADCell],
) -> Witness | None:
    """Select a satisfying candidate while preserving fixed coordinates.

    Candidate samples are checked after the fixed assignment is applied. No
    fiber is relifted; failure to find a candidate does not prove infeasibility.
    Source-cell provenance is retained only for an unchanged sample point.
    """
    expression = to_sympy(formula)
    for cell in candidate_cells:
        sample = {
            var: sample_to_expr(value)
            for var, value in zip(vars_[: cell.level], cell.sample, strict=True)
        }
        point = {**sample, **free_assignment}
        try:
            truth = exact_truth(expression.subs(point))
        except ValueError:
            continue
        if truth:
            return Witness(point, cell if point == sample else None)
    return None


__all__ = ["Witness", "witness_from_qe_result", "existential_wit_assign"]
