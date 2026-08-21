from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass

import sympy as sp

from ..cache_utils import BoundedLRU
from ..normalization import conjuncts
from .solution import IntervalComponent

_RECOVERABLE_ERRORS = (
    ArithmeticError,
    TypeError,
    ValueError,
    NotImplementedError,
    RuntimeError,
    sp.PolynomialError,
)


@dataclass(frozen=True)
class MetadataRequest:
    """Control which potentially expensive solution metadata is computed."""

    cells: bool = False
    cylindrical: bool = False
    connectivity: bool = False

    @property
    def structural(self) -> bool:
        return self.cells or self.cylindrical or self.connectivity


_METADATA_CACHE: BoundedLRU[dict[str, object]] = BoundedLRU(128, "decision.metadata")


def _empty_metadata() -> dict[str, object]:
    return {
        "dimension": None,
        "bounded": None,
        "closed": None,
        "compact": None,
        "components": (),
        "cells": (),
        "cylindrical_solution": None,
        "connectivity": None,
    }


def interval_components_from_set(set_expr: sp.Set, var: sp.Symbol) -> tuple[IntervalComponent, ...]:
    if set_expr is sp.S.EmptySet or set_expr == sp.S.EmptySet:
        return ()
    pieces = set_expr.args if isinstance(set_expr, sp.Union) else (set_expr,)
    components: list[IntervalComponent] = []
    for piece in pieces:
        if isinstance(piece, sp.Interval):
            components.append(
                IntervalComponent(
                    var,
                    piece.start,
                    piece.end,
                    not bool(piece.left_open),
                    not bool(piece.right_open),
                )
            )
        elif isinstance(piece, sp.FiniteSet):
            for point in sorted(piece, key=sp.default_sort_key):
                components.append(IntervalComponent(var, point, point, True, True))
        else:
            raise NotImplementedError(f"unsupported one-dimensional solution set piece: {piece!r}")
    components.sort(
        key=lambda comp: (sp.default_sort_key(comp.lower), sp.default_sort_key(comp.upper))
    )
    return tuple(components)


def one_dim_components(expr: sp.Expr, var: sp.Symbol) -> tuple[IntervalComponent, ...] | None:
    if expr is sp.false or expr == sp.false:
        return ()
    if expr is sp.true or expr == sp.true:
        return (IntervalComponent(var, -sp.oo, sp.oo, False, False),)
    try:
        set_expr = expr.as_set()
        if isinstance(set_expr, sp.ConditionSet):
            return None
        return interval_components_from_set(set_expr, var)
    except (TypeError, ValueError, NotImplementedError):
        return None


def components_formula(components: Sequence[IntervalComponent]) -> sp.Expr:
    if not components:
        return sp.false
    formulas = [component.as_formula() for component in components]
    return sp.Or(*formulas) if len(formulas) > 1 else formulas[0]


def has_strict_atom(expr: sp.Expr) -> bool:
    if isinstance(expr, (sp.StrictLessThan, sp.StrictGreaterThan)):
        return True
    if isinstance(expr, (sp.And, sp.Or)):
        return any(has_strict_atom(arg) for arg in expr.args)
    return False


def one_dim_bounds(expr: sp.Expr, var: sp.Symbol) -> tuple[sp.Expr, sp.Expr] | None:
    lo: sp.Expr = -sp.oo
    hi: sp.Expr = sp.oo
    try:
        reduced = sp.reduce_inequalities(list(conjuncts(expr)), var)
    except (TypeError, ValueError, NotImplementedError):
        reduced = expr
    for atom in conjuncts(reduced):
        if not isinstance(atom, sp.core.relational.Relational):
            continue
        lhs, rhs = atom.lhs, atom.rhs
        if isinstance(atom, (sp.LessThan, sp.StrictLessThan)):
            if lhs == var and not rhs.has(var):
                hi = sp.Min(hi, rhs) if hi != sp.oo else rhs
            elif rhs == var and not lhs.has(var):
                lo = sp.Max(lo, lhs) if lo != -sp.oo else lhs
        elif isinstance(atom, (sp.GreaterThan, sp.StrictGreaterThan)):
            if lhs == var and not rhs.has(var):
                lo = sp.Max(lo, rhs) if lo != -sp.oo else rhs
            elif rhs == var and not lhs.has(var):
                hi = sp.Min(hi, lhs) if hi != sp.oo else lhs
        elif isinstance(atom, sp.Equality):
            if lhs == var and not rhs.has(var):
                lo = hi = rhs
            elif rhs == var and not lhs.has(var):
                lo = hi = lhs
    return (lo, hi)


def _update_from_components(
    metadata: dict[str, object], formula: sp.Expr, variables: tuple[sp.Symbol, ...]
) -> bool:
    if len(variables) != 1:
        return False
    components = one_dim_components(formula, variables[0])
    if components is not None:
        if not components:
            metadata.update(
                dimension=None,
                bounded=True,
                closed=True,
                compact=True,
                components=(),
            )
        else:
            dimension = max(component.dimension for component in components)
            bounded = all(component.bounded for component in components)
            closed = all(component.closed for component in components)
            metadata.update(
                dimension=dimension,
                bounded=bounded,
                closed=closed,
                compact=bounded and closed,
                components=components,
            )
        return True
    bounds = one_dim_bounds(formula, variables[0])
    if bounds is not None:
        lo, hi = bounds
        finite = lo != -sp.oo and hi != sp.oo
        closed = not has_strict_atom(formula)
        metadata.update(dimension=1, bounded=finite, closed=closed, compact=finite and closed)
    return False


def _update_from_box(
    metadata: dict[str, object], formula: sp.Expr, variables: tuple[sp.Symbol, ...]
) -> None:
    try:
        from ..implicit_geometry import extract_symbolic_box_bounds

        box = extract_symbolic_box_bounds(formula, variables)
    except _RECOVERABLE_ERRORS:
        return
    if box is None:
        return
    finite = all(lo != -sp.oo and hi != sp.oo for _, lo, hi in box.limits)
    closed = not has_strict_atom(formula)
    metadata.update(
        dimension=len(variables),
        bounded=finite,
        closed=closed,
        compact=finite and closed,
    )


def _finite_vertical_cells(cells: Sequence[object]) -> bool:
    return all(
        cell.x_interval[0] != -sp.oo
        and cell.x_interval[1] != sp.oo
        and all(lower != -sp.oo and upper != sp.oo for lower, upper in cell.y_bounds)
        for cell in cells
    )


def _update_from_vertical_cells(
    metadata: dict[str, object], formula: sp.Expr, variables: tuple[sp.Symbol, ...]
) -> None:
    if len(variables) != 2:
        return
    cells_all: tuple[object, ...] = ()
    try:
        from ..implicit_geometry import decompose_cylindrical_formula_to_vertical_bounds_2d

        cells_all = tuple(decompose_cylindrical_formula_to_vertical_bounds_2d(formula, variables))
    except _RECOVERABLE_ERRORS:
        try:
            from ..cad.cells import extract_vertical_bounds_from_cad_2d

            cells_all = tuple(
                extract_vertical_bounds_from_cad_2d(formula, variables, full_dimensional_only=False)
            )
        except _RECOVERABLE_ERRORS:
            return
    cells = tuple(cell for cell in cells_all if getattr(cell, "dimension", 2) == 2) or cells_all
    if not cells:
        return
    metadata["cells"] = cells
    metadata["dimension"] = max(getattr(cell, "dimension", 2) for cell in cells)
    finite = _finite_vertical_cells(cells)
    closed = not has_strict_atom(formula)
    metadata.update(bounded=finite, closed=closed, compact=finite and closed)


def _has_nonlinear_second_level(
    formula: sp.Expr, variables: tuple[sp.Symbol, ...], metadata: dict[str, object]
) -> bool:
    if len(variables) < 2 or not metadata.get("cells"):
        return False
    y_var = variables[1]
    for atom in getattr(formula, "args", (formula,)):
        if not getattr(atom, "is_Relational", False) or y_var not in getattr(
            atom, "free_symbols", set()
        ):
            continue
        try:
            if sp.Poly(atom.lhs - atom.rhs, y_var).degree() > 1:
                return True
        except _RECOVERABLE_ERRORS:
            continue
    return False


def _update_from_cylindrical(
    metadata: dict[str, object],
    formula: sp.Expr,
    variables: tuple[sp.Symbol, ...],
    *,
    include_connectivity: bool,
) -> None:
    if len(variables) < 2 or _has_nonlinear_second_level(formula, variables, metadata):
        return
    try:
        from ..cad.cells import extract_cylindrical_solution, extract_explicit_cylindrical_solution

        cyl = extract_explicit_cylindrical_solution(formula, variables)
        if cyl is None:
            cyl = extract_cylindrical_solution(formula, variables, selected_only=True)
    except _RECOVERABLE_ERRORS:
        return
    metadata["cylindrical_solution"] = cyl
    if not getattr(cyl, "cells", ()):
        return
    metadata["dimension"] = cyl.dimension
    metadata["bounded"] = cyl.bounded
    if not metadata.get("cells"):
        metadata["cells"] = cyl.cells
    if not include_connectivity:
        return
    try:
        from ..connectivity import build_cad_adjacency_graph

        connectivity = build_cad_adjacency_graph(cyl, formula=formula)
    except _RECOVERABLE_ERRORS:
        return
    metadata["connectivity"] = connectivity
    if connectivity.components:
        metadata["components"] = connectivity.components


def _infer_dimension(
    metadata: dict[str, object], formula: sp.Expr, variables: tuple[sp.Symbol, ...]
) -> None:
    if metadata["dimension"] is not None:
        return
    equalities = [atom for atom in conjuncts(formula) if isinstance(atom, sp.Equality)]
    metadata["dimension"] = (
        max(0, len(variables) - len(equalities)) if equalities else len(variables)
    )


def collect_solution_metadata(
    formula: sp.Expr,
    variables: tuple[sp.Symbol, ...],
    *,
    request: MetadataRequest | None = None,
) -> dict[str, object]:
    """Collect cached solution metadata, computing structural CAD data only on demand."""

    request = request or MetadataRequest()
    key = (formula, variables, request.cells, request.cylindrical, request.connectivity)
    cached = _METADATA_CACHE.get(key)
    if cached is not None:
        return dict(cached)

    metadata = _empty_metadata()
    if formula is sp.false or formula == sp.false:
        metadata.update(dimension=None, bounded=True, closed=True, compact=True)
    elif formula is sp.true or formula == sp.true:
        metadata.update(
            dimension=len(variables),
            bounded=len(variables) == 0,
            closed=True,
            compact=len(variables) == 0,
        )
    else:
        finished_1d = _update_from_components(metadata, formula, variables)
        if not finished_1d:
            _update_from_box(metadata, formula, variables)
        if request.cells or request.cylindrical or request.connectivity:
            _update_from_vertical_cells(metadata, formula, variables)
        if request.cylindrical or request.connectivity:
            _update_from_cylindrical(
                metadata,
                formula,
                variables,
                include_connectivity=request.connectivity,
            )
        _infer_dimension(metadata, formula, variables)

    _METADATA_CACHE.put(key, dict(metadata))
    return metadata


def metadata_request_for_output(
    output: str | None,
    sample_mode: str,
) -> MetadataRequest:
    """Resolve structural metadata needs from the requested public view."""

    if output is None:
        # Preserve the structured-result contract: cells and a cylindrical view
        # remain populated when supported, while connectivity is opt-in.
        return MetadataRequest(
            cells=True,
            cylindrical=True,
            connectivity=sample_mode == "per_component",
        )
    key = output.lower().replace("-", "_")
    connectivity_outputs = {
        "connectivity",
        "adjacency",
        "roadmap",
        "roadmap_graph",
        "components_graph",
    }
    cylindrical_outputs = {"cylindrical", "cylindrical_solution", "cylindrical_cells", "cad_cells"}
    cell_outputs = {"cells", "cell"}
    if key in connectivity_outputs or sample_mode == "per_component":
        return MetadataRequest(cells=True, cylindrical=True, connectivity=True)
    if key in cylindrical_outputs:
        return MetadataRequest(cells=True, cylindrical=True)
    if key in cell_outputs or sample_mode == "per_cell":
        return MetadataRequest(cells=True)
    return MetadataRequest()


def clear_solution_metadata_cache() -> None:
    _METADATA_CACHE.clear()


__all__ = [
    "MetadataRequest",
    "clear_solution_metadata_cache",
    "collect_solution_metadata",
    "components_formula",
    "has_strict_atom",
    "interval_components_from_set",
    "metadata_request_for_output",
    "one_dim_bounds",
    "one_dim_components",
]
