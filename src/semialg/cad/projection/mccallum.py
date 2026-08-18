from __future__ import annotations

from collections.abc import Sequence

import sympy as sp

from .reduced import ReducedProjectionTower, build_reduced_proj_tower


def build_mccallum_proj_set(
    polys: Sequence[sp.Expr | sp.Poly],
    variables: Sequence[sp.Symbol],
    *,
    equational_constraints: Sequence[sp.Expr] = (),
    certify: bool = False,
) -> ReducedProjectionTower:
    """Build a McCallum-style reduced projection tower for diagnostics.

    The returned validity object must be consulted before using the tower as a
    complete backend. The The safe projection driver safe driver falls back to Collins unless
    the tower is explicitly certified.
    """

    return build_reduced_proj_tower(
        polys,
        variables,
        theory="mccallum",
        equational_constraints=equational_constraints,
        certify=certify,
    )
