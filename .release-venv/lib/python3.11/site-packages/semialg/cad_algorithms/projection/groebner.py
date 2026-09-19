"""Groebner-structured projection data for zero-dimensional equality varieties.

This module does not replace complete sign-invariant CAD.  It constructs a
proof-carrying *variety projection* when a formula has common polynomial
equational constraints whose ideal is zero-dimensional.  The resulting
triangular lexicographic basis is aligned with the CAD variable order so that
lifting can follow only algebraic sections on which all common equalities may
vanish.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass

import sympy as sp

from ...algebraic.equality_ideal import EqualityIdealContext
from ...algebraic.rational_univariate import RationalUnivariateError
from ...formula import Formula, equational_constraints, to_sympy
from .collins import ProjectionLevel, ProjectionPolynomial, ProjectionTower, normalize_poly


@dataclass(frozen=True)
class GroebnerVarietyProjection:
    """Triangular projection data for a finite equality variety."""

    equality_context: EqualityIdealContext
    tower: ProjectionTower
    residual_formula: sp.Expr
    triangular_basis: tuple[sp.Expr, ...]
    _coordinate_polynomials: tuple[sp.Expr, ...] | None = None

    @property
    def coordinate_polynomials(self) -> tuple[sp.Expr, ...]:
        """Materialize coordinate eliminants lazily for diagnostics or expert use."""

        if self._coordinate_polynomials is not None:
            return self._coordinate_polynomials
        return self.equality_context.fglm_coordinate_polynomials()

    @property
    def quotient_dimension(self) -> int:
        return int(self.equality_context.quotient_dimension or 0)


def _poly_level(expr: sp.Expr, variables: tuple[sp.Symbol, ...]) -> int:
    symbols = expr.free_symbols.intersection(variables)
    if not symbols:
        return 0
    return max(variables.index(symbol) + 1 for symbol in symbols)


def _triangular_basis(
    context: EqualityIdealContext,
    variables: tuple[sp.Symbol, ...],
) -> tuple[sp.Expr, ...]:
    """Return the context-cached FGLM basis triangular in CAD lifting order."""

    return context.triangular_basis(variables)


def _build_tower(
    triangular_basis: Sequence[sp.Expr],
    variables: tuple[sp.Symbol, ...],
    *,
    quotient_dimension: int,
) -> ProjectionTower:
    """Build a projection tower from a triangular Groebner basis without duplicating level polynomials."""
    by_level: dict[int, list[sp.Poly]] = {level: [] for level in range(1, len(variables) + 1)}
    entries: dict[int, list[ProjectionPolynomial]] = {
        level: [] for level in range(1, len(variables) + 1)
    }
    for expr in triangular_basis:
        level = _poly_level(expr, variables)
        if level == 0:
            continue
        poly = normalize_poly(sp.Poly(expr, *variables, extension=True))
        if any(existing.as_expr() == poly.as_expr() for existing in by_level[level]):
            continue
        by_level[level].append(poly)
        entries[level].append(
            ProjectionPolynomial(
                poly=poly,
                level=level,
                source="groebner-triangular",
                operation_variable=variables[level - 1],
                expression=expr,
            )
        )
    if any(not by_level[level] for level in range(1, len(variables) + 1)):
        raise RationalUnivariateError(
            "triangular FGLM basis does not constrain every CAD lifting level"
        )
    levels = tuple(
        ProjectionLevel(
            level,
            variables[level - 1],
            tuple(by_level[level]),
            tuple(entries[level]),
        )
        for level in range(1, len(variables) + 1)
    )
    original = tuple(poly for level in levels for poly in level.polynomials)
    return ProjectionTower(
        variables=variables,
        levels=levels,
        original_polynomials=original,
        metadata={
            "projection": "groebner-variety",
            "variety_only": True,
            "complete_space_decomposition": False,
            "equality_dimension": 0,
            "quotient_dimension": quotient_dimension,
            "coordinate_polynomials": "lazy",
            "triangular_basis": tuple(triangular_basis),
        },
    )


def _empty_variety_tower(variables: tuple[sp.Symbol, ...]) -> ProjectionTower:
    levels = tuple(
        ProjectionLevel(level, variable, tuple(), tuple())
        for level, variable in enumerate(variables, start=1)
    )
    return ProjectionTower(
        variables=variables,
        levels=levels,
        original_polynomials=tuple(),
        metadata={
            "projection": "groebner-variety",
            "variety_only": True,
            "complete_space_decomposition": False,
            "equality_dimension": -1,
            "quotient_dimension": 0,
            "coordinate_polynomials": tuple(),
            "triangular_basis": (sp.Integer(1),),
            "inconsistent_equality_ideal": True,
        },
    )


def build_groebner_variety_projection(
    formula: Formula,
    variables: Sequence[sp.Symbol],
) -> GroebnerVarietyProjection | None:
    """Build finite-variety projection data when exact prerequisites hold.

    The equational constraints must be common to every Boolean branch.  If the
    exact equality ideal is absent, unsupported, inconsistent, or
    positive-dimensional, ``None`` is returned and callers should use their
    generic CAD backend.
    """

    vars_tuple = tuple(variables)
    equations = tuple(equational_constraints(formula))
    if not equations or not vars_tuple:
        return None
    try:
        context = EqualityIdealContext(equations, vars_tuple)
    except (
        RationalUnivariateError,
        sp.PolynomialError,
        sp.polys.polyerrors.CoercionFailed,
        TypeError,
        ValueError,
        NotImplementedError,
    ):
        return None
    if context.inconsistent:
        return GroebnerVarietyProjection(
            equality_context=context,
            tower=_empty_variety_tower(vars_tuple),
            residual_formula=sp.false,
            triangular_basis=(sp.Integer(1),),
            _coordinate_polynomials=tuple(),
        )
    if not context.zero_dimensional:
        return None
    try:
        triangular = _triangular_basis(context, vars_tuple)
        residual = context.simplify_constraints(to_sympy(formula))
        tower = _build_tower(
            triangular,
            vars_tuple,
            quotient_dimension=int(context.quotient_dimension or 0),
        )
    except (
        RationalUnivariateError,
        sp.PolynomialError,
        sp.polys.polyerrors.CoercionFailed,
        TypeError,
        ValueError,
        NotImplementedError,
    ):
        return None
    return GroebnerVarietyProjection(
        equality_context=context,
        tower=tower,
        residual_formula=sp.sympify(residual),
        triangular_basis=triangular,
        _coordinate_polynomials=None,
    )


__all__ = ["GroebnerVarietyProjection", "build_groebner_variety_projection"]
