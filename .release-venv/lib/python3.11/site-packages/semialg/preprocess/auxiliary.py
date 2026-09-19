from __future__ import annotations

from dataclasses import dataclass
from typing import Literal

import sympy as sp

AuxKind = Literal["abs", "rational_power"]


@dataclass(frozen=True)
class AuxiliaryDef:
    """Provenance for an auxiliary introduced by semialgebraic preprocessing."""

    symbol: sp.Symbol
    expression: sp.Expr
    kind: AuxKind
    constraints: tuple[sp.Expr, ...]
    branch_points: tuple[sp.Expr, ...] = ()


class AuxiliaryFactory:
    """Deterministic fresh-symbol factory avoiding user variables."""

    def __init__(self, existing: set[sp.Symbol] | None = None, prefix: str = "_sa") -> None:
        self.existing = set(existing or set())
        self.prefix = prefix
        self.counter = 0

    def fresh(self, stem: str) -> sp.Symbol:
        while True:
            sym = sp.Symbol(f"{self.prefix}_{stem}_{self.counter}", real=True)
            self.counter += 1
            if sym not in self.existing:
                self.existing.add(sym)
                return sym


__all__ = ["AuxKind", "AuxiliaryDef", "AuxiliaryFactory"]
