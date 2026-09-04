"""Solver-facing access to exact polynomial equality-ideal analysis."""

from __future__ import annotations

from collections.abc import Iterable, Sequence

import sympy as sp

from ..algebraic.equality_ideal import (
    EqualityIdealAnalysis,
    EqualityIdealContext,
    LinearVariableElimination,
    ZeroDimensionalFilterResult,
    analyze_equality_ideal,
)
from ..algebraic.rational_univariate import RationalUnivariateError


def build_equality_ideal_context(
    equations: Iterable[sp.Expr | sp.Equality | bool],
    variables: Sequence[sp.Symbol],
) -> EqualityIdealContext | None:
    """Build an exact equality-ideal context, or return ``None`` if unsupported."""

    try:
        return EqualityIdealContext(equations, variables)
    except (
        RationalUnivariateError,
        sp.PolynomialError,
        sp.polys.polyerrors.CoercionFailed,
        TypeError,
        ValueError,
        NotImplementedError,
    ):
        return None


__all__ = [
    "EqualityIdealAnalysis",
    "EqualityIdealContext",
    "LinearVariableElimination",
    "ZeroDimensionalFilterResult",
    "analyze_equality_ideal",
    "build_equality_ideal_context",
]
