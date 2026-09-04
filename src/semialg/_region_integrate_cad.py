from __future__ import annotations

from collections.abc import Mapping, Sequence
from itertools import permutations

import sympy as sp

from ._errors import EXACT_OPERATION_ERRORS as _RECOVERABLE_ERRORS
from ._region_integrate_geometry import _radial_radii_squared, _vertical_slice_data
from .implicit_geometry import decompose_cylindrical_formula_to_vertical_bounds_2d
from .interval_decomposition import finite_real_roots as _finite_real_roots  # noqa: F401
from .region_integral_results import (
    RegionIntegralPiece,
)


def _reduce_radial_vertical_pieces(
    integrand: sp.Expr,
    condition: sp.Expr,
    x: sp.Symbol,
    y: sp.Symbol,
) -> tuple[RegionIntegralPiece, ...] | None:
    radii = _radial_radii_squared(condition, x, y)
    if radii is None:
        return None
    lower_sq, upper_sq = radii
    if lower_sq == upper_sq:
        return ()
    if upper_sq == sp.oo:
        raise NotImplementedError(
            "unbounded radial regions are outside finite iterated-piece reduction"
        )

    pieces: list[RegionIntegralPiece] = []

    def add_disk_piece(radius_sq: sp.Expr, signed_integrand: sp.Expr, label: str) -> None:
        if sp.simplify(radius_sq) == 0:
            return
        radius = sp.sqrt(radius_sq)
        height = sp.sqrt(radius_sq - x**2)
        pieces.append(
            RegionIntegralPiece(
                integrand=signed_integrand,
                limits=((y, -height, height), (x, -radius, radius)),
                method=label,
                diagnostics={"radius_squared": radius_sq},
            )
        )

    add_disk_piece(upper_sq, integrand, "radial_region_outer_disk_vertical_slice")
    if sp.simplify(lower_sq) != 0:
        add_disk_piece(lower_sq, -integrand, "radial_region_inner_disk_subtraction")
    return tuple(pieces)


def _reduce_vertical_slice(
    integrand: sp.Expr,
    condition: sp.Expr,
    x: sp.Symbol,
    y: sp.Symbol,
    bounds: Mapping[sp.Symbol, tuple[sp.Expr, sp.Expr]],
) -> tuple[RegionIntegralPiece, ...] | None:
    data = _vertical_slice_data(condition, x, y, bounds)
    if data is None:
        return None
    lower, upper, intervals = data
    return tuple(
        RegionIntegralPiece(
            integrand=integrand,
            limits=((y, lower, upper), (x, lo, hi)),
            method="vertical_slice_iterated_integral",
            diagnostics={"lower": lower, "upper": upper},
        )
        for lo, hi in intervals
    )


def _reduce_cylindrical_vertical_bounds_2d(
    integrand: sp.Expr,
    condition: sp.Expr,
    x: sp.Symbol,
    y: sp.Symbol,
) -> tuple[RegionIntegralPiece, ...] | None:
    """Reduce supported CAD-like 2D cylindrical formulas to vertical pieces."""

    try:
        cells = decompose_cylindrical_formula_to_vertical_bounds_2d(condition, (x, y))
    except NotImplementedError:
        return None
    pieces: list[RegionIntegralPiece] = []
    for cell in cells:
        lo, hi = cell.x_interval
        if lo == hi:
            continue
        if lo == -sp.oo or hi == sp.oo:
            return None
        for lower, upper in cell.y_bounds:
            if lower == upper:
                continue
            if lower == -sp.oo or upper == sp.oo:
                return None
            pieces.append(
                RegionIntegralPiece(
                    integrand=integrand,
                    limits=((y, lower, upper), (x, lo, hi)),
                    method="cylindrical_formula_vertical_bounds_2d",
                    diagnostics={
                        "x_interval": (lo, hi),
                        "y_bounds": (lower, upper),
                        "source_formula": cell.source_formula,
                    },
                )
            )
    return tuple(pieces)


def _reduce_cad_extracted_vertical_bounds_2d(
    integrand: sp.Expr,
    condition: sp.Expr,
    x: sp.Symbol,
    y: sp.Symbol,
) -> tuple[RegionIntegralPiece, ...] | None:
    """Reduce arbitrary complete-CAD 2D cells to vertical integral pieces."""

    try:
        from .cad_algorithms.cells import extract_vertical_bounds_from_cad_2d

        cells = extract_vertical_bounds_from_cad_2d(condition, (x, y), full_dimensional_only=True)
    except _RECOVERABLE_ERRORS:
        return None
    pieces: list[RegionIntegralPiece] = []
    for cell in cells:
        lo, hi = cell.x_interval
        if lo == hi or lo == -sp.oo or hi == sp.oo:
            continue
        for lower, upper in cell.y_bounds:
            if lower == upper or lower == -sp.oo or upper == sp.oo:
                continue
            pieces.append(
                RegionIntegralPiece(
                    integrand=integrand,
                    limits=((y, lower, upper), (x, lo, hi)),
                    method="complete_cad_vertical_bounds_2d",
                    diagnostics={
                        "x_interval": (lo, hi),
                        "y_bounds": (lower, upper),
                        "source_formula": cell.source_formula,
                    },
                )
            )
    return tuple(pieces) if pieces else None


def _reduce_cylindrical_solution_cells_nd(
    integrand: sp.Expr,
    condition: sp.Expr,
    variables: Sequence[sp.Symbol],
) -> tuple[RegionIntegralPiece, ...] | None:
    """Reduce arbitrary-dimensional full CAD cells to nested iterated integrals.

    This exploits the public cylindrical-solution representation. It is exact
    for full-dimensional cells whose coordinate bounds can be expressed by the
    current CAD extractor. Lower-dimensional cells are ignored for ambient
    Lebesgue integration.
    """

    try:
        from .cad_algorithms.cells import (
            extract_cylindrical_solution,
            extract_explicit_cylindrical_solution,
        )
        from .cad_algorithms.integration import full_dimensional_cell_integral

        cyl = extract_explicit_cylindrical_solution(condition, variables)
        if cyl is None:
            cyl = extract_cylindrical_solution(condition, variables, selected_only=True)
    except _RECOVERABLE_ERRORS:
        return None
    pieces: list[RegionIntegralPiece] = []
    n = len(tuple(variables))
    for cell in getattr(cyl, "cells", ()):
        if getattr(cell, "dimension", None) != n:
            continue
        try:
            adapted = full_dimensional_cell_integral(cell, integrand, require_verified=True)
        except _RECOVERABLE_ERRORS:
            continue
        limits = adapted.limits
        if len(limits) != n:
            continue
        pieces.append(
            RegionIntegralPiece(
                integrand=integrand,
                limits=limits,
                method="cylindrical_solution_cell_iterated_integral",
                diagnostics={
                    "cell_index": getattr(cell, "index", None),
                    "cell_dimension": getattr(cell, "dimension", None),
                    "typed_bounds_verified": adapted.certified_bounds,
                    "integration_variable_order": tuple(variables),
                },
            )
        )
    return tuple(pieces) if pieces else None


def _radical_penalty(expr: sp.Expr) -> int:
    """Cheap score favoring polynomial/rational CAD integration bounds."""

    penalty = 0
    for power in sp.preorder_traversal(sp.sympify(expr)):
        if isinstance(power, sp.Pow) and power.exp.is_Rational and power.exp.q != 1:
            penalty += 1
        elif isinstance(power, sp.RootOf):
            penalty += 4
    return penalty


def _integration_piece_score(
    pieces: Sequence[RegionIntegralPiece],
    order: Sequence[sp.Symbol],
) -> tuple[int, int, int, tuple[str, ...]]:
    """Rank exact decompositions by cell count and symbolic bound complexity."""

    radicals = 0
    operations = 0
    for piece in pieces:
        for _var, lower, upper in piece.limits:
            radicals += _radical_penalty(lower) + _radical_penalty(upper)
            operations += int(sp.count_ops(lower)) + int(sp.count_ops(upper))
    return (len(pieces), radicals, operations, tuple(v.name for v in order))


def _explicit_integration_orders(
    variables: tuple[sp.Symbol, ...],
) -> tuple[tuple[sp.Symbol, ...], ...]:
    """Return bounded coordinate permutations worth trying without a CAD build."""

    if len(variables) <= 1:
        return (variables,)
    if len(variables) <= 3:
        return tuple(permutations(variables))
    return (variables,)


def _reduce_with_best_explicit_cylindrical_order(
    integrand: sp.Expr,
    condition: sp.Expr,
    variables: tuple[sp.Symbol, ...],
) -> tuple[tuple[RegionIntegralPiece, ...], tuple[sp.Symbol, ...]] | None:
    """Choose the simplest exact explicit cylindrical coordinate order.

    This is deliberately tried before launching a CAD.  A region written in a
    poor caller order can already be cylindrical after permuting coordinates;
    selecting that order avoids introducing algebraic root bounds solely as an
    artifact of orientation.
    """

    candidates: list[
        tuple[
            tuple[int, int, int, tuple[str, ...]],
            tuple[sp.Symbol, ...],
            tuple[RegionIntegralPiece, ...],
        ]
    ] = []
    for order in _explicit_integration_orders(variables):
        try:
            from .cad_algorithms.cells import extract_explicit_cylindrical_solution
            from .cad_algorithms.integration import full_dimensional_cell_integral

            solution = extract_explicit_cylindrical_solution(condition, order)
        except _RECOVERABLE_ERRORS:
            continue
        if solution is None:
            continue
        pieces: list[RegionIntegralPiece] = []
        for cell in solution.full_dimensional_cells:
            try:
                adapted = full_dimensional_cell_integral(cell, integrand, require_verified=True)
            except _RECOVERABLE_ERRORS:
                pieces = []
                break
            pieces.append(
                RegionIntegralPiece(
                    integrand=integrand,
                    limits=adapted.limits,
                    method="coordinate_permuted_cylindrical_integral",
                    diagnostics={
                        "cell_index": getattr(cell, "index", None),
                        "typed_bounds_verified": adapted.certified_bounds,
                        "integration_variable_order": order,
                        "order_source": "explicit_cylindrical_search",
                    },
                )
            )
        if pieces:
            pieces_tuple = tuple(pieces)
            candidates.append((_integration_piece_score(pieces_tuple, order), order, pieces_tuple))
    if not candidates:
        return None
    _score, order, pieces = min(candidates, key=lambda item: item[0])
    return pieces, order


def _cad_integration_orders(
    condition: sp.Expr,
    variables: tuple[sp.Symbol, ...],
) -> tuple[tuple[sp.Symbol, ...], ...]:
    """Return a small deterministic set of CAD lifting orders for integration."""

    orders: list[tuple[sp.Symbol, ...]] = [variables]
    if len(variables) == 2:
        orders.append(tuple(reversed(variables)))
    try:
        from .formula import formula_polynomials, parse_formula
        from .heuristics import suggest_cad_variable_order

        parsed = parse_formula(condition)
        polys = tuple(formula_polynomials(parsed))
        if polys:
            suggested = suggest_cad_variable_order(polys, variables, strategy="auto").order
            orders.append(tuple(suggested))
    except _RECOVERABLE_ERRORS:
        pass
    unique: list[tuple[sp.Symbol, ...]] = []
    for order in orders:
        if order not in unique:
            unique.append(order)
    return tuple(unique)


def _reduce_with_best_cad_order(
    integrand: sp.Expr,
    condition: sp.Expr,
    variables: tuple[sp.Symbol, ...],
) -> tuple[tuple[RegionIntegralPiece, ...], tuple[sp.Symbol, ...]] | None:
    """Build certified CAD cell integrals using the best bounded order search."""

    candidates: list[
        tuple[
            tuple[int, int, int, tuple[str, ...]],
            tuple[sp.Symbol, ...],
            tuple[RegionIntegralPiece, ...],
        ]
    ] = []
    for order in _cad_integration_orders(condition, variables):
        pieces = _reduce_cylindrical_solution_cells_nd(integrand, condition, order)
        if pieces is None:
            continue
        candidates.append((_integration_piece_score(pieces, order), order, pieces))
    if not candidates:
        return None
    _score, order, pieces = min(candidates, key=lambda item: item[0])
    annotated = tuple(
        RegionIntegralPiece(
            integrand=piece.integrand,
            limits=piece.limits,
            method="cad_variable_order_cell_integral",
            diagnostics={
                **(piece.diagnostics or {}),
                "integration_variable_order": order,
                "order_source": "cad_order_search",
            },
        )
        for piece in pieces
    )
    return annotated, order
