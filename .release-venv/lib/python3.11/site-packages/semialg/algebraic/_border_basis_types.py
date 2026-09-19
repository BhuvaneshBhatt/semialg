"""Representation and diagnostics for exact border-basis construction."""

from __future__ import annotations

from dataclasses import dataclass

import sympy as sp

from .rational_univariate.representation import RationalUnivariateError


class BorderBasisError(RationalUnivariateError):
    """Raised when an exact border basis cannot be constructed."""


@dataclass(frozen=True)
class BorderBasisDiagnostics:
    """Diagnostics for exact border-basis construction.

    The diagnostics are exact/symbolic. Numerical border-basis
    algorithms should report residual norms separately; this class records the
    structural checks used by the current exact implementation.
    """

    success: bool = True
    messages: tuple[str, ...] = tuple()
    quotient_basis_rank: int | None = None
    border_rank: int | None = None
    commutators_zero: bool | None = None
    failed_border_monomial: sp.Expr | None = None
    failed_expression: sp.Expr | None = None

    def add_message(self, message: str) -> BorderBasisDiagnostics:
        return BorderBasisDiagnostics(
            success=self.success,
            messages=self.messages + (message,),
            quotient_basis_rank=self.quotient_basis_rank,
            border_rank=self.border_rank,
            commutators_zero=self.commutators_zero,
            failed_border_monomial=self.failed_border_monomial,
            failed_expression=self.failed_expression,
        )
