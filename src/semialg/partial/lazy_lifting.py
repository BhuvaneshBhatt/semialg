from __future__ import annotations

from collections.abc import Sequence
from functools import cmp_to_key

import sympy as sp

from ..cad import _real_root_list, _roots_over_cell, _stack_from_roots
from ..exact_arithmetic import compare_exact_reals
from ..model import CADResult, Cell, ProjectionConfig


def ensure_children(
    cad: CADResult,
    parent: Cell,
    vars_: Sequence[sp.Symbol],
    projection_config: ProjectionConfig,
    collins_projection_sets=None,
):
    vars_ = tuple(vars_)
    if parent.index in cad.children_by_parent:
        return cad.children_by_parent[parent.index]
    child_level = parent.level + 1
    if child_level > len(vars_):
        cad.children_by_parent[parent.index] = []
        return []
    if child_level == 1:
        roots = []
        x1 = vars_[0]
        for poly in cad.projection_sets[1]:
            if x1 in poly.free_symbols:
                roots.extend(_real_root_list(poly.as_expr(), x1, diagnostics=cad.diagnostics))
        unique = []
        for root in roots:
            if not any(sp.simplify(root - other) == 0 for other in unique):
                unique.append(root)
        unique.sort(key=cmp_to_key(compare_exact_reals))
        stack = _stack_from_roots(parent, unique)
    else:
        target_var = vars_[child_level - 1]
        fallback_level = (
            collins_projection_sets[child_level] if collins_projection_sets is not None else None
        )
        roots, events, used_fallback, supplemented, certificate = _roots_over_cell(
            cad.projection_sets[child_level],
            vars_,
            parent.sample,
            target_var,
            parent_cell=parent,
            config=projection_config,
            fallback_level_polys=fallback_level,
            projection_info=cad.projection_metadata.get(child_level),
            diagnostics=cad.diagnostics,
        )
        cad.lifting_certificates = tuple(list(cad.lifting_certificates) + [certificate])
        cad.lifting_supplemented = cad.lifting_supplemented or supplemented
        if events:
            cad.nullification_events = tuple(list(cad.nullification_events) + list(events))
            cad.well_oriented = False
        if used_fallback:
            cad.used_fallback_projection = True
        stack = _stack_from_roots(parent, roots)
    cad.children_by_parent[parent.index] = stack
    cad.cells_by_level.setdefault(child_level, []).extend(stack)
    cad.lazy_cells_by_level.setdefault(child_level, []).extend(stack)
    return stack


__all__ = ["ensure_children"]
