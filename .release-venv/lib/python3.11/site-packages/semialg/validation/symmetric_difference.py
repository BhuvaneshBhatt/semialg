from __future__ import annotations

import itertools
from collections.abc import Iterable, Sequence
from dataclasses import dataclass

import sympy as sp

from .solution_checking import form_sat_by_assign

_DEFAULT_COORDS = (
    sp.Integer(-3),
    sp.Integer(-2),
    sp.Integer(-1),
    sp.Rational(-1, 2),
    sp.Integer(0),
    sp.Rational(1, 2),
    sp.Integer(1),
    sp.Integer(2),
    sp.Integer(3),
)


@dataclass(frozen=True)
class SymDiffWit:
    assignment: dict[sp.Symbol, sp.Expr]
    left_value: bool
    right_value: bool


@dataclass(frozen=True)
class SymmetricDifferenceCheck:
    equivalent_on_grid: bool
    checked_points: int
    witness: SymDiffWit | None = None


def find_grid_witness(
    left: sp.Expr,
    right: sp.Expr,
    variables: Sequence[sp.Symbol],
    *,
    coordinates: Iterable[sp.Expr] = _DEFAULT_COORDS,
) -> SymmetricDifferenceCheck:
    coords = tuple(coordinates)
    checked = 0
    for values in itertools.product(coords, repeat=len(variables)):
        assignment = dict(zip(variables, values, strict=True))
        checked += 1
        left_value = form_sat_by_assign(left, assignment, check_numeric_equalities=True)
        right_value = form_sat_by_assign(right, assignment, check_numeric_equalities=True)
        if left_value != right_value:
            return SymmetricDifferenceCheck(
                equivalent_on_grid=False,
                checked_points=checked,
                witness=SymDiffWit(assignment, left_value, right_value),
            )
    return SymmetricDifferenceCheck(equivalent_on_grid=True, checked_points=checked)


__all__ = ["SymDiffWit", "SymmetricDifferenceCheck", "find_grid_witness"]
