"""Dependency-aware exact specialization of cylindrical CAD samples."""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field

import sympy as sp

from ..algebraic.samples import sample_to_expr
from ..reconstruct.root_functions import root_of


@dataclass
class CylindricalSampleContext:
    """Close a CAD sample tower strictly in cylindrical dependency order.

    Coordinate ``i`` may depend only on coordinates ``< i``.  Prefix closure is
    memoized so incidence, sign evaluation and topology share one specialization
    rule instead of open-coding nested ``root_of`` substitution.
    """

    variables: tuple[sp.Symbol, ...]
    raw_values: tuple[sp.Expr, ...]
    _prefixes: dict[int, tuple[tuple[sp.Symbol, sp.Expr], ...]] = field(default_factory=dict)

    @classmethod
    def from_cell(cls, cell, variables: Sequence[sp.Symbol]) -> CylindricalSampleContext:
        vars_ = tuple(variables)
        return cls(vars_, tuple(sample_to_expr(s) for s in cell.sample))

    @classmethod
    def from_mapping(cls, assignments: Mapping[sp.Symbol, sp.Expr]) -> CylindricalSampleContext:
        return cls(tuple(assignments), tuple(map(sp.sympify, assignments.values())))

    def prefix_assignment(self, level: int | None = None) -> dict[sp.Symbol, sp.Expr]:
        """Return closed coordinates strictly below ``level`` (all if omitted)."""
        stop = len(self.variables) if level is None else int(level)
        if stop < 0 or stop > len(self.variables):
            raise ValueError("cylindrical sample prefix level is out of range")
        cached = self._prefixes.get(stop)
        if cached is not None:
            return dict(cached)
        resolved: dict[sp.Symbol, sp.Expr] = {}
        # Local import avoids a module cycle; specialization itself consumes only
        # an already-closed lower prefix.
        from .incidence import _specialize_root_functions

        for var, value in zip(self.variables[:stop], self.raw_values[:stop], strict=True):
            closed = sp.sympify(value).subs(resolved)
            if closed.has(root_of):
                closed = _specialize_root_functions(closed, resolved)
                closed = sp.sympify(closed).subs(resolved)
            unresolved = closed.free_symbols.intersection(self.variables[stop:])
            if unresolved:
                names = ", ".join(sorted(s.name for s in unresolved))
                raise ValueError(f"cylindrical sample depends on a higher coordinate: {names}")
            resolved[var] = closed
        frozen = tuple(resolved.items())
        self._prefixes[stop] = frozen
        return dict(frozen)

    def assignment_before(self, fiber: sp.Symbol) -> dict[sp.Symbol, sp.Expr]:
        try:
            level = self.variables.index(fiber)
        except ValueError:
            # A recursive specialization receives only an already-resolved lower
            # prefix; in that context every supplied coordinate precedes fiber.
            return self.prefix_assignment()
        return self.prefix_assignment(level)

    def full_assignment(self) -> dict[sp.Symbol, sp.Expr]:
        return self.prefix_assignment()
