from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass

import sympy as sp

from ..cad.lifting.stack import CADCell
from .cylindrical import _sector_bound_expr, level_cell_condition, section_value_expr


@dataclass(frozen=True)
class FormulaCompressionStats:
    """Small diagnostic record for cylindrical formula compression."""

    selected_cells: int
    emitted_blocks: int
    fallback_used: bool = False


@dataclass(frozen=True)
class NestedFormulaResult:
    """Result of recursive cylindrical formula reconstruction."""

    formula: sp.Expr
    stats: FormulaCompressionStats


def _children_by_parent(
    cells_by_level: Mapping[int, Sequence[CADCell]],
) -> dict[tuple[int, ...] | None, tuple[CADCell, ...]]:
    grouped: dict[tuple[int, ...] | None, list[CADCell]] = {}
    for level_cells in cells_by_level.values():
        for cell in level_cells:
            grouped.setdefault(cell.parent_index, []).append(cell)
    return {
        parent: tuple(sorted(children, key=lambda cell: cell.stack_position))
        for parent, children in grouped.items()
    }


def _leaf_indices_under(
    parent: tuple[int, ...], leaf_indices: set[tuple[int, ...]]
) -> set[tuple[int, ...]]:
    plen = len(parent)
    return {idx for idx in leaf_indices if idx[:plen] == parent}


def _all_leaf_indices_under(
    parent: tuple[int, ...],
    children: Mapping[tuple[int, ...] | None, tuple[CADCell, ...]],
    full_level: int,
) -> set[tuple[int, ...]]:
    if len(parent) == full_level:
        return {parent}
    out: set[tuple[int, ...]] = set()
    for child in children.get(parent, ()):
        out.update(_all_leaf_indices_under(child.index, children, full_level))
    return out


def _is_full_subtree(
    cell: CADCell,
    leaf_indices: set[tuple[int, ...]],
    children: Mapping[tuple[int, ...] | None, tuple[CADCell, ...]],
    full_level: int,
) -> bool:
    return _all_leaf_indices_under(cell.index, children, full_level).issubset(leaf_indices)


def _blocks(cells: Sequence[CADCell]) -> list[list[CADCell]]:
    if not cells:
        return []
    sorted_cells = sorted(cells, key=lambda cell: cell.stack_position)
    out: list[list[CADCell]] = [[sorted_cells[0]]]
    for cell in sorted_cells[1:]:
        if cell.stack_position == out[-1][-1].stack_position + 1:
            out[-1].append(cell)
        else:
            out.append([cell])
    return out


def _block_condition(
    block: Sequence[CADCell],
    variable: sp.Symbol,
    cells_by_level: Mapping[int, Sequence[CADCell]],
    *,
    closed: bool,
) -> sp.Expr:
    """Return a compact condition for a consecutive stack block.

    A block such as sector-section-sector is emitted as one closed or half-open
    interval instead of three disjuncts. Non-consecutive cells are handled by
    the caller as separate blocks.
    """

    if not block:
        return sp.false
    if len(block) == 1:
        return level_cell_condition(block[0], variable, cells_by_level, closed=closed)
    first = block[0]
    last = block[-1]
    left_sample = first.interval[0] if first.interval is not None else None
    right_sample = last.interval[1] if last.interval is not None else None
    pieces: list[sp.Expr] = []
    if left_sample is not None:
        if first.kind == "section":
            pieces.append(variable >= section_value_expr(first, variable))
        else:
            lower = _sector_bound_expr(first, variable, cells_by_level, side="left")
            if lower is not None:
                pieces.append(variable > lower)
    if right_sample is not None:
        if last.kind == "section":
            pieces.append(variable <= section_value_expr(last, variable))
        else:
            upper = _sector_bound_expr(last, variable, cells_by_level, side="right")
            if upper is not None:
                pieces.append(variable < upper)
    return sp.And(*pieces) if pieces else sp.true


def nested_formula_from_cells(
    cells: Sequence[CADCell],
    variables: Sequence[sp.Symbol],
    cells_by_level: Mapping[int, Sequence[CADCell]],
    *,
    max_terms: int = 512,
    closed: bool = False,
) -> NestedFormulaResult:
    """Compress a union of CAD leaf cells into a nested cylindrical formula.

    The semantic input is the cell set; the emitted expression is only a display
    and downstream simplification form. The compressor shares common prefixes,
    collapses complete subtrees to their parent condition, and merges contiguous
    stack positions into interval blocks where possible.
    """

    full_level = len(tuple(variables))
    if not cells:
        return NestedFormulaResult(sp.false, FormulaCompressionStats(0, 0))
    leaf_indices = {cell.index for cell in cells}
    children = _children_by_parent(cells_by_level)
    emitted = 0

    def rec(parent: tuple[int, ...] | None, level: int) -> sp.Expr:
        nonlocal emitted
        if emitted > max_terms:
            return sp.false
        if level > full_level:
            return sp.true
        candidates = [
            cell
            for cell in children.get(parent, ())
            if _leaf_indices_under(cell.index, leaf_indices)
        ]
        if not candidates:
            return sp.false
        var = variables[level - 1]
        formulas: list[sp.Expr] = []
        full_cells: list[CADCell] = []
        partial_cells: list[CADCell] = []
        for cell in candidates:
            if _is_full_subtree(cell, leaf_indices, children, full_level):
                full_cells.append(cell)
            else:
                partial_cells.append(cell)
        # Consecutive full subtrees can be emitted as a single interval block.
        for block in _blocks(full_cells):
            emitted += 1
            formulas.append(_block_condition(block, var, cells_by_level, closed=closed))
        for cell in partial_cells:
            emitted += 1
            head = level_cell_condition(cell, var, cells_by_level, closed=closed)
            tail = rec(cell.index, level + 1)
            if tail is sp.false or tail == sp.false:
                continue
            formulas.append(sp.And(head, tail) if head is not sp.true and head != sp.true else tail)
        if not formulas:
            return sp.false
        if len(formulas) == 1:
            return formulas[0]
        return sp.Or(*formulas)

    formula = rec(None, 1)
    fallback = emitted > max_terms
    if fallback:
        # The caller can fall back to DNF if needed; returning false silently
        # would be unsound, so expose the fallback condition in the stats.
        formula = sp.false
    return NestedFormulaResult(
        formula,
        FormulaCompressionStats(
            selected_cells=len(cells),
            emitted_blocks=min(emitted, len(cells)),
            fallback_used=fallback,
        ),
    )


__all__ = ["FormulaCompressionStats", "NestedFormulaResult", "nested_formula_from_cells"]
