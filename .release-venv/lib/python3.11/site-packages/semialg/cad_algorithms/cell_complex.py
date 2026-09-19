"""Exact incidence complexes extracted from selected CAD cells."""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass

import sympy as sp

from ..decision import is_satisfiable
from ..normalization import normalize_formula
from .cylindrical_solution import (
    CylindricalSolution,
    CylindricalSolutionCell,
    extract_cylindrical_solution,
)


@dataclass(frozen=True)
class CADIncidence:
    lower: int
    higher: int
    lower_index: tuple[int, ...]
    higher_index: tuple[int, ...]
    codimension: int = 1


@dataclass(frozen=True)
class CADCellComplex:
    """Finite exact cell complex induced by the selected cells of one CAD.

    ``incidences`` records codimension-one closure incidence.  Transitive
    incidence is obtained by following the boundary/coboundary graph.
    """

    variables: tuple[sp.Symbol, ...]
    formula: sp.Expr
    cells: tuple[CylindricalSolutionCell, ...]
    incidences: tuple[CADIncidence, ...]
    solution: CylindricalSolution

    def cells_by_dimension(self, dimension: int) -> tuple[CylindricalSolutionCell, ...]:
        return tuple(cell for cell in self.cells if cell.dimension == dimension)

    @property
    def dimensions(self) -> tuple[int, ...]:
        return tuple(sorted({cell.dimension for cell in self.cells}))

    @property
    def f_vector(self) -> tuple[int, ...]:
        ambient = len(self.variables)
        return tuple(len(self.cells_by_dimension(d)) for d in range(ambient + 1))

    def boundary_indices(self, cell: int | tuple[int, ...]) -> tuple[int, ...]:
        pos = self._position(cell)
        return tuple(sorted(i.lower for i in self.incidences if i.higher == pos))

    def coboundary_indices(self, cell: int | tuple[int, ...]) -> tuple[int, ...]:
        pos = self._position(cell)
        return tuple(sorted(i.higher for i in self.incidences if i.lower == pos))

    def boundary_of(self, cell: int | tuple[int, ...]) -> tuple[CylindricalSolutionCell, ...]:
        return tuple(self.cells[i] for i in self.boundary_indices(cell))

    def coboundary_of(self, cell: int | tuple[int, ...]) -> tuple[CylindricalSolutionCell, ...]:
        return tuple(self.cells[i] for i in self.coboundary_indices(cell))

    def incidence_matrix(self, dimension: int) -> sp.ImmutableMatrix:
        """Unsigned codimension-one incidence matrix C_d -> C_{d-1}."""
        high = [i for i, c in enumerate(self.cells) if c.dimension == dimension]
        low = [i for i, c in enumerate(self.cells) if c.dimension == dimension - 1]
        hpos = {idx: j for j, idx in enumerate(high)}
        lpos = {idx: i for i, idx in enumerate(low)}
        data = [[0] * len(high) for _ in low]
        for inc in self.incidences:
            if inc.lower in lpos and inc.higher in hpos:
                data[lpos[inc.lower]][hpos[inc.higher]] = 1
        return sp.ImmutableMatrix(data)

    def euler_characteristic(self) -> int:
        return int(sum((-1) ** cell.dimension for cell in self.cells))

    def _position(self, cell: int | tuple[int, ...]) -> int:
        if isinstance(cell, int):
            if cell < 0 or cell >= len(self.cells):
                raise IndexError(cell)
            return cell
        target = tuple(cell)
        for i, item in enumerate(self.cells):
            if item.index == target:
                return i
        raise KeyError(target)


def _formula_subset(lhs: sp.Expr, rhs: sp.Expr, variables: Sequence[sp.Symbol]) -> bool:
    witness = sp.And(lhs, sp.Not(rhs))
    return not bool(is_satisfiable(witness, tuple(variables), strategy="cad"))


def build_cad_cell_complex(
    region: object,
    variables: Sequence[sp.Symbol] | None = None,
    *,
    solution: CylindricalSolution | None = None,
) -> CADCellComplex:
    """Build exact codimension-one incidence from a selected CAD decomposition."""
    if solution is None:
        if hasattr(region, "result") and hasattr(region.result, "variables"):
            result = region.result
            formula = normalize_formula(result.formula)
            vars_ = tuple(result.variables)
            solution = extract_cylindrical_solution(result, selected_only=True)
        else:
            formula = normalize_formula(getattr(region, "formula", region))
            vars_ = tuple(getattr(region, "variables", ()) if variables is None else variables)
            solution = extract_cylindrical_solution(formula, vars_, selected_only=True)
    else:
        formula = normalize_formula(solution.formula)
        vars_ = tuple(solution.variables)

    cells = tuple(solution.cells)
    incidences: list[CADIncidence] = []
    by_dim: dict[int, list[int]] = {}
    for i, cell in enumerate(cells):
        by_dim.setdefault(cell.dimension, []).append(i)

    # Prefer the source CAD's exact recursive closure test.  Reconstructing a
    # closed cylindrical formula can introduce radicals such as Abs(x), which
    # are mathematically harmless but no longer polynomial input for CAD/QE.
    source_result = None
    if solution.source_decomposition is not None:
        source_result = solution.source_decomposition.cad_result
    source_cells = {}
    if source_result is not None:
        source_cells = {
            cell.index: cell for cell in source_result.cad.cells_by_level.get(len(vars_), ())
        }

    for d in range(1, len(vars_) + 1):
        for hi in by_dim.get(d, ()):
            for lo in by_dim.get(d - 1, ()):
                incident = False
                lower_source = source_cells.get(cells[lo].index)
                higher_source = source_cells.get(cells[hi].index)
                if lower_source is not None and higher_source is not None:
                    from ..topology.incidence import is_cell_in_closure

                    incident = is_cell_in_closure(
                        lower_source, higher_source, vars_, source_result.cad.cells_by_level
                    )
                else:
                    high_closed = cells[hi].as_formula(closed=True)
                    incident = _formula_subset(cells[lo].as_formula(), high_closed, vars_)
                if incident:
                    incidences.append(CADIncidence(lo, hi, cells[lo].index, cells[hi].index))
    return CADCellComplex(vars_, formula, cells, tuple(incidences), solution)


__all__ = ["CADIncidence", "CADCellComplex", "build_cad_cell_complex"]
