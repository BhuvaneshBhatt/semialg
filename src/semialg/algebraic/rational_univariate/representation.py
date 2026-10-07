from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass

import sympy as sp
from algroots.rational_univariate import RationalUnivariatePoint as _AlgrootsRationalUnivariatePoint
from algroots.rational_univariate import RationalUnivariateRepresentation

from ...dimension_validation import assignments_from_points, zip_equal
from ...exceptions import AlgebraicSolvingError
from ...status import SolverStatus


class RationalUnivariateError(AlgebraicSolvingError):
    """Semialg-facing error for rational-univariate backend failures."""


@dataclass(frozen=True)
class RationalUnivariatePoint(_AlgrootsRationalUnivariatePoint):
    """An algroots RUR point enriched with Semialg sign/Thom operations."""

    @property
    def assignment(self) -> Mapping[sp.Symbol, sp.Expr]:
        return dict(zip_equal(self.variables, self.coordinates, context="RUR point coordinates"))

    @property
    def thom_encoding(self):
        """Return the Thom certificate identifying the parameter root."""

        from .thom import rur_parameter_thom_encoding

        return rur_parameter_thom_encoding(self)

    def sign_of(self, polynomial: sp.Poly | sp.Expr) -> int:
        """Return the exact sign of a polynomial at this represented point."""

        from .thom import sign_at_rur_point

        return sign_at_rur_point(polynomial, self)


@dataclass(frozen=True)
class RationalUnivariateFormulaResult:
    """RUR-backed solutions for a Boolean formula with finite equality branches."""

    variables: tuple[sp.Symbol, ...]
    assignments: tuple[Mapping[sp.Symbol, sp.Expr], ...]
    backend: str = "rational-univariate-formula-solver"
    status: SolverStatus | str = SolverStatus.SAT
    solved_branches: int = 0
    skipped_branches: int = 0
    notes: tuple[str, ...] = ()

    @property
    def points(self) -> tuple[tuple[sp.Expr, ...], ...]:
        return tuple(
            tuple(assignment[var] for var in self.variables) for assignment in self.assignments
        )

    @property
    def satisfiable(self) -> bool:
        return bool(self.assignments)

    @property
    def complete(self) -> bool:
        return self.skipped_branches == 0

    @property
    def partial(self) -> bool:
        return self.skipped_branches > 0


@dataclass(frozen=True)
class FilteredRationalUnivariateSolutions:
    """Solutions of a zero-dimensional equality system filtered by constraints."""

    variables: tuple[sp.Symbol, ...]
    representation: RationalUnivariateRepresentation
    points: tuple[tuple[sp.Expr, ...], ...]

    @property
    def assignments(self) -> tuple[Mapping[sp.Symbol, sp.Expr], ...]:
        return assignments_from_points(self.variables, self.points, context="RUR solution point")

    @property
    def satisfiable(self) -> bool:
        return bool(self.points)


__all__ = [
    "RationalUnivariateError",
    "RationalUnivariateRepresentation",
    "RationalUnivariatePoint",
    "RationalUnivariateFormulaResult",
    "FilteredRationalUnivariateSolutions",
]
