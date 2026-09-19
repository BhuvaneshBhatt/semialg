from __future__ import annotations

from collections.abc import Mapping

import sympy as sp


def qe_by_complete_cad(*args, **kwargs):
    """Lazy exact-QE seam shared by optimization submodules."""

    from .qe import qe_by_complete_cad as impl

    return impl(*args, **kwargs)


def optimization_is_feasible(condition: sp.Expr, point: Mapping[sp.Symbol, sp.Expr]) -> bool:
    """Lazy feasibility seam avoiding an optimization import cycle."""

    from .optimization import _is_feasible

    return _is_feasible(condition, point)


__all__ = ["optimization_is_feasible", "qe_by_complete_cad"]
