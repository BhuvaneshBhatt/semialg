"""Reusable exact point-location and sign-vector queries on a complete CAD."""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass

import sympy as sp

from ..exact_arithmetic import compare_exact_reals


@dataclass(frozen=True)
class CADSignEntry:
    level: int
    polynomial: sp.Expr
    sign: int


@dataclass(frozen=True)
class CADPointLocation:
    variables: tuple[sp.Symbol, ...]
    point: Mapping[sp.Symbol, sp.Expr]
    cell_index: tuple[int, ...]
    selected: bool
    dimension: int
    signs: tuple[CADSignEntry, ...]

    @property
    def sign_vector(self) -> tuple[int, ...]:
        return tuple(entry.sign for entry in self.signs)


def _point_mapping(variables, point) -> dict[sp.Symbol, sp.Expr]:
    vars_ = tuple(variables)
    if isinstance(point, Mapping):
        missing = tuple(v for v in vars_ if v not in point)
        if missing:
            names = ", ".join(sp.sstr(v) for v in missing)
            raise ValueError(f"missing coordinate(s): {names}")
        values = {v: sp.sympify(point[v]) for v in vars_}
    else:
        seq = tuple(point)
        if len(seq) != len(vars_):
            raise ValueError(f"expected {len(vars_)} coordinate(s), got {len(seq)}")
        values = {v: sp.sympify(x) for v, x in zip(vars_, seq, strict=True)}
    if any(x.free_symbols for x in values.values()):
        raise ValueError("CAD point location requires exact resolved coordinates")
    return values


def _exact_sign(expr: sp.Expr) -> int:
    value = sp.simplify(expr)
    if value == 0:
        return 0
    if value.is_positive is True:
        return 1
    if value.is_negative is True:
        return -1
    return compare_exact_reals(value, sp.Integer(0))


def locate_cad_point(region: object, point) -> CADPointLocation:
    """Locate ``point`` once in a reusable CAD and return its tower sign vector."""
    from ..cad_region import _leaf_dimension, as_cad_region

    reg = as_cad_region(region)
    values = _point_mapping(reg.variables, point)
    function = reg.result.as_function()
    leaf = function.locate_cell(values)
    if leaf is None:
        raise ValueError("point could not be located in the complete CAD")
    node = function.tree_by_index().get(leaf.index)
    if node is None or node.truth is None:
        raise ValueError("located CAD leaf is missing a certified truth value")
    selected = node.truth
    entries: list[CADSignEntry] = []
    for level in sorted(reg.cad.tower.by_level):
        for poly in sorted(
            reg.cad.tower.by_level[level],
            key=lambda p: sp.srepr(p.as_expr()),
        ):
            sign = _exact_sign(poly.as_expr().subs(values))
            entries.append(CADSignEntry(level, poly.as_expr(), sign))
    return CADPointLocation(
        tuple(reg.variables),
        values,
        tuple(leaf.index),
        selected,
        _leaf_dimension(leaf, reg.result),
        tuple(entries),
    )


__all__ = ["CADSignEntry", "CADPointLocation", "locate_cad_point"]
