from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass

import sympy as sp


@dataclass(frozen=True)
class RegionIntegralResult:
    """Integral over a supported semialgebraic region.

    Integrals use ambient Lebesgue measure by default. Lower-dimensional
    equality-only subsets therefore contribute zero unless the caller requests
    ``measure_dimension="intrinsic"`` or an explicit Hausdorff dimension.

    ``exact`` is true for symbolic evaluation and false for numerical
    evaluation. ``evaluated`` records whether the returned value has been
    evaluated rather than left as an explicit ``Integral`` expression.
    """

    value: sp.Expr
    integrand: sp.Expr
    condition: sp.Expr
    variables: tuple[sp.Symbol, ...]
    method: str
    diagnostics: Mapping[str, object] | None = None
    exact: bool = True
    evaluated: bool = True
    error_estimate: sp.Expr | None = None


@dataclass(frozen=True)
class RegionIntegralPiece:
    """One iterated-integral piece for a semialgebraic region integral.

    ``limits`` uses SymPy's integral-limit convention, for example
    ``((x, a, b), (y, lower(x), upper(x)))``. Pieces may be signed: Boolean
    differences such as annuli are represented by using a negative integrand on
    the subtracted part.
    """

    integrand: sp.Expr
    limits: tuple[tuple[sp.Symbol, sp.Expr, sp.Expr], ...]
    method: str
    diagnostics: Mapping[str, object] | None = None

    def as_integral(self) -> sp.Integral:
        """Return this piece as an unevaluated SymPy ``Integral``."""

        return sp.Integral(self.integrand, *self.limits)


@dataclass(frozen=True)
class ReducedRegionIntegral:
    """Reduction of a region integral to explicit iterated-integral pieces."""

    integrand: sp.Expr
    condition: sp.Expr
    variables: tuple[sp.Symbol, ...]
    pieces: tuple[RegionIntegralPiece, ...]
    method: str
    diagnostics: Mapping[str, object] | None = None

    def as_integrals(self) -> tuple[sp.Integral, ...]:
        """Return each reduced piece as an unevaluated SymPy ``Integral``."""

        return tuple(piece.as_integral() for piece in self.pieces)

    def unevaluated_sum(self) -> sp.Expr:
        """Return the formal sum of unevaluated piece integrals."""

        if not self.pieces:
            return sp.Integer(0)
        return sp.Add(*self.as_integrals())


__all__ = ["RegionIntegralResult", "RegionIntegralPiece", "ReducedRegionIntegral"]
