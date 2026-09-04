from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass

import sympy as sp

from ..decomposition.cylindrical import CADResult, CellSet
from ..normalization import normalize_formula, normalize_problem_variables
from ..structural_keys import symbol_identity_key
from .bounds import (
    AlgebraicRootFunction,
    CADBound,
    DelineabilityCertificate,
    RootOrderCertificate,
    as_cad_bound,
)
from .structured_cells import (
    StructuredCADCell,
    StructuredCADCellDecomposition,
    extract_structured_cad_cells,
)


@dataclass(frozen=True)
class CylindricalCoordinateConstraint:
    """One coordinate constraint in a cylindrical solution cell.

    The constraint is interpreted relative to the preceding coordinates in the
    same cell. For example, in variables ``(x, y, z)`` the third level may have
    bounds involving both ``x`` and ``y``. Sections represent equalities
    ``variable == lower == upper``; sectors represent open intervals by default
    and may be rendered with weak inequalities when ``closed=True`` is requested.
    """

    variable: sp.Symbol
    level: int
    kind: str
    lower: sp.Expr
    upper: sp.Expr
    sample: sp.Expr
    index: tuple[int, ...]
    lower_bound: CADBound | None = None
    upper_bound: CADBound | None = None
    lower_closed: bool = False
    upper_closed: bool = False
    delineability: DelineabilityCertificate | None = None
    root_order: RootOrderCertificate | None = None

    @property
    def typed_lower(self) -> CADBound:
        return self.lower_bound or as_cad_bound(self.lower, closed=self.lower_closed)

    @property
    def typed_upper(self) -> CADBound:
        return self.upper_bound or as_cad_bound(self.upper, closed=self.upper_closed)

    @property
    def is_section(self) -> bool:
        return self.kind == "section"

    @property
    def is_sector(self) -> bool:
        return self.kind == "sector"

    @property
    def dimension(self) -> int:
        return 0 if self.is_section else 1

    def as_formula(self, *, closed: bool = False) -> sp.Expr:
        if self.is_section:
            return sp.Eq(self.variable, self.lower)
        parts: list[sp.Expr] = []
        if self.lower != -sp.oo:
            lower_closed = closed or self.lower_closed
            parts.append(
                self.variable >= self.lower if lower_closed else self.variable > self.lower
            )
        if self.upper != sp.oo:
            upper_closed = closed or self.upper_closed
            parts.append(
                self.variable <= self.upper if upper_closed else self.variable < self.upper
            )
        return sp.And(*parts) if parts else sp.true

    def as_limit(self) -> tuple[sp.Symbol, sp.Expr, sp.Expr]:
        """Return ``(variable, lower, upper)`` for iterated-integral style use."""

        return (self.variable, self.lower, self.upper)


@dataclass(frozen=True)
class CylindricalSolutionCell:
    """One nested cylindrical solution cell.

    A cell is an ordered path through a CAD: first a constraint on ``x1``, then
    a constraint on ``x2`` whose bounds may depend on ``x1``, and so on. This
    is the user-facing representation for solutions such as ``x`` in a base
    cell and ``y`` between algebraic functions over that base cell.
    """

    variables: tuple[sp.Symbol, ...]
    levels: tuple[CylindricalCoordinateConstraint, ...]
    sample: Mapping[sp.Symbol, sp.Expr]
    index: tuple[int, ...]
    selected: bool = True
    source_cell: StructuredCADCell | None = None

    @property
    def dimension(self) -> int:
        return sum(level.dimension for level in self.levels)

    @property
    def codimension(self) -> int:
        return len(self.variables) - self.dimension

    @property
    def bounded(self) -> bool:
        return all(level.lower != -sp.oo and level.upper != sp.oo for level in self.levels)

    @property
    def is_full_dimensional(self) -> bool:
        return self.dimension == len(self.variables)

    def as_formula(self, *, closed: bool = False) -> sp.Expr:
        parts = [level.as_formula(closed=closed) for level in self.levels]
        return sp.And(*parts) if parts else sp.true

    def sample_point(self) -> dict[sp.Symbol, sp.Expr]:
        return {sym: sp.simplify(value) for sym, value in self.sample.items()}

    def iterated_limits(self) -> tuple[tuple[sp.Symbol, sp.Expr, sp.Expr], ...]:
        """Return nested expression limits in CAD variable order."""

        return tuple(level.as_limit() for level in self.levels)

    def cylindrical_bounds(self) -> tuple[tuple[sp.Symbol, CADBound, CADBound], ...]:
        """Return typed nested bounds preserving algebraic root functions."""

        return tuple(
            (level.variable, level.typed_lower, level.typed_upper) for level in self.levels
        )

    def verify_bounds(self):
        from .bounds import verify_cad_cell_bounds

        return verify_cad_cell_bounds(self)


@dataclass(frozen=True)
class CylindricalDecompositionCertificate:
    """Certificate that solution cells form a complete disjoint decomposition."""

    coverage_verified: bool
    pairwise_disjoint: bool
    cells_verified: bool
    source: str = "unknown"
    notes: tuple[str, ...] = ()

    @property
    def certified(self) -> bool:
        return self.coverage_verified and self.pairwise_disjoint and self.cells_verified

    def verify(self) -> bool:
        return self.certified


@dataclass(frozen=True)
class CylindricalSolution:
    """Full cylindrical solution representation for a semialgebraic formula."""

    variables: tuple[sp.Symbol, ...]
    cells: tuple[CylindricalSolutionCell, ...]
    formula: sp.Expr
    source_decomposition: StructuredCADCellDecomposition | None = None
    decomposition_cert: CylindricalDecompositionCertificate | None = None

    @property
    def dimension(self) -> int | None:
        if not self.cells:
            return None
        return max(cell.dimension for cell in self.cells)

    @property
    def bounded(self) -> bool | None:
        if not self.cells:
            return True
        return all(cell.bounded for cell in self.cells)

    @property
    def full_dimensional_cells(self) -> tuple[CylindricalSolutionCell, ...]:
        n = len(self.variables)
        return tuple(cell for cell in self.cells if cell.dimension == n)

    @property
    def sample_points(self) -> tuple[Mapping[sp.Symbol, sp.Expr], ...]:
        return tuple(cell.sample_point() for cell in self.cells)

    def cells_by_dimension(self, dimension: int) -> tuple[CylindricalSolutionCell, ...]:
        return tuple(cell for cell in self.cells if cell.dimension == dimension)

    def as_formula(self, *, closed: bool = False) -> sp.Expr:
        if not self.cells:
            return sp.false
        pieces = [cell.as_formula(closed=closed) for cell in self.cells]
        return sp.Or(*pieces) if len(pieces) > 1 else pieces[0]


def _cylindrical_cell_from_structured(cell: StructuredCADCell) -> CylindricalSolutionCell:
    levels = tuple(
        CylindricalCoordinateConstraint(
            variable=level.variable,
            level=level.level,
            kind=level.kind,
            lower=level.lower,
            upper=level.upper,
            sample=level.sample,
            index=level.index,
            lower_bound=level.lower_bound,
            upper_bound=level.upper_bound,
            lower_closed=level.lower_closed,
            upper_closed=level.upper_closed,
            delineability=level.delineability,
            root_order=level.root_order,
        )
        for level in cell.levels
    )
    return CylindricalSolutionCell(
        variables=cell.variables,
        levels=levels,
        sample=cell.sample,
        index=cell.index,
        selected=cell.selected,
        source_cell=cell,
    )


def cylindrical_solution_from_structured(
    decomposition: StructuredCADCellDecomposition,
) -> CylindricalSolution:
    """Convert structured CAD cells to the public cylindrical solution form."""

    cells = tuple(_cylindrical_cell_from_structured(cell) for cell in decomposition.cells)
    from .bounds import verify_cad_cell_bounds

    cert = CylindricalDecompositionCertificate(
        coverage_verified=True,
        pairwise_disjoint=True,
        cells_verified=all(verify_cad_cell_bounds(cell).verify() for cell in cells),
        source="cad",
        notes=("CAD leaf cells form a cylindrical partition of the selected formula",),
    )
    return CylindricalSolution(
        variables=decomposition.variables,
        cells=cells,
        formula=decomposition.formula,
        source_decomposition=decomposition,
        decomposition_cert=cert,
    )


def extract_cylindrical_solution(
    condition_or_cad: object,
    variables: Sequence[sp.Symbol | str] | None = None,
    *,
    selected_only: bool = True,
) -> CylindricalSolution:
    """Return a nested cylindrical representation of the selected solution cells.

    This is the high-level public form of the CAD cell extraction layer. It can
    represent arbitrary CAD paths in any dimension: each cell stores coordinate
    bounds in variable order, and bounds at level ``k`` may involve variables
    from levels ``< k``.
    """

    if not isinstance(condition_or_cad, (CADResult, CellSet)):
        formula = normalize_formula(condition_or_cad)
        vars_ = (
            normalize_problem_variables(variables, formula)
            if variables is not None
            else tuple(sorted(formula.free_symbols, key=symbol_identity_key))
        )
        from .explicit_cells import extract_explicit_cylindrical_solution

        explicit = extract_explicit_cylindrical_solution(formula, vars_)
        if explicit is not None and any(
            isinstance(bound, AlgebraicRootFunction)
            for cell in explicit.cells
            for level in cell.levels
            for bound in (level.lower_bound, level.upper_bound)
            if bound is not None
        ):
            return explicit
    decomp = extract_structured_cad_cells(condition_or_cad, variables, selected_only=selected_only)
    return cylindrical_solution_from_structured(decomp)


__all__ = [
    "CylindricalCoordinateConstraint",
    "CylindricalSolutionCell",
    "CylindricalDecompositionCertificate",
    "CylindricalSolution",
    "cylindrical_solution_from_structured",
    "extract_cylindrical_solution",
]
