from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass

import sympy as sp

from ...dimension_validation import zip_equal
from ...exceptions import AlgebraicSolvingError
from ...status import SolverStatus


class RationalUnivariateError(AlgebraicSolvingError):
    """Raised when a rational univariate representation cannot be computed."""


@dataclass(frozen=True)
class RationalUnivariateRepresentation:
    """Rational univariate representation of a zero-dimensional system.

    If ``parameter`` is named ``t``, each solution is represented by
    ``defining_polynomial(t) == 0`` and
    ``variable_i == coordinate_numerators[i](t) / coordinate_denominator(t)``.
    The implementation targets rational-coefficient, zero-dimensional polynomial
    systems and returns distinct algebraic solution branches.
    """

    variables: tuple[sp.Symbol, ...]
    parameter: sp.Symbol
    defining_polynomial: sp.Poly
    coordinate_denominator: sp.Poly
    coordinate_numerators: tuple[sp.Poly, ...]
    separating_linear_form: sp.Expr
    standard_exponents: tuple[tuple[int, ...], ...]
    quotient_dimension: int | None = None
    geometric_solution_count: int | None = None

    @property
    def dimension(self) -> int:
        return (
            self.quotient_dimension
            if self.quotient_dimension is not None
            else len(self.standard_exponents)
        )

    @property
    def solution_count(self) -> int:
        return (
            self.geometric_solution_count
            if self.geometric_solution_count is not None
            else self.defining_polynomial.degree()
        )

    @property
    def is_empty(self) -> bool:
        return self.defining_polynomial.degree() <= 0

    def coordinate_expressions(self) -> tuple[sp.Expr, ...]:
        denominator = self.coordinate_denominator.as_expr()
        return tuple(
            sp.cancel(numer.as_expr() / denominator) for numer in self.coordinate_numerators
        )

    def normalized_coordinate_polynomials(self) -> tuple[sp.Poly, ...]:
        """Return denominator-free coordinate polynomials modulo the RUR polynomial."""

        if self.is_empty:
            return tuple(sp.Poly(0, self.parameter, domain=sp.QQ) for _ in self.variables)
        try:
            inverse_denominator = sp.invert(self.coordinate_denominator, self.defining_polynomial)
        except (sp.polys.polyerrors.NotInvertible, sp.PolynomialError, ValueError) as exc:
            raise RationalUnivariateError(
                "coordinate denominator is not invertible modulo the RUR polynomial"
            ) from exc
        return tuple(
            sp.Poly(
                (inverse_denominator * numerator).rem(self.defining_polynomial).as_expr(),
                self.parameter,
                domain=sp.QQ,
            )
            for numerator in self.coordinate_numerators
        )


@dataclass(frozen=True)
class RationalUnivariatePoint:
    """A point represented by a RUR parameter root.

    ``root`` is an exact real root of ``representation.defining_polynomial``.
    Coordinates are obtained by evaluating the denominator-free coordinate
    polynomials modulo the RUR defining polynomial. This representation lets
    sign queries reduce multivariate expressions to a single exact algebraic
    sign computation in the RUR parameter.
    """

    representation: RationalUnivariateRepresentation
    root: sp.Expr

    @property
    def variables(self) -> tuple[sp.Symbol, ...]:
        return self.representation.variables

    @property
    def coordinates(self) -> tuple[sp.Expr, ...]:
        t = self.representation.parameter
        return tuple(
            sp.cancel(poly.as_expr().subs(t, self.root))
            for poly in self.representation.normalized_coordinate_polynomials()
        )

    @property
    def assignment(self) -> Mapping[sp.Symbol, sp.Expr]:
        return dict(zip_equal(self.variables, self.coordinates, context="RUR point coordinates"))


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
    """Solutions of a zero-dimensional equality system filtered by constraints.

    ``points`` stores tuples in the same order as ``variables``. ``assignments``
    stores the same solutions as symbol-to-expression dictionaries, which is
    convenient for witness-generation callers.
    """

    variables: tuple[sp.Symbol, ...]
    representation: RationalUnivariateRepresentation
    points: tuple[tuple[sp.Expr, ...], ...]

    @property
    def assignments(self) -> tuple[Mapping[sp.Symbol, sp.Expr], ...]:
        return tuple(
            dict(zip_equal(self.variables, point, context="RUR solution point"))
            for point in self.points
        )

    @property
    def satisfiable(self) -> bool:
        return bool(self.points)
