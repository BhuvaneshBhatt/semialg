from __future__ import annotations

from collections.abc import Sequence

import sympy as sp

from .reduced import ReducedProjectionTower, build_reduced_proj_tower


def build_lazard_proj_set(
    polys: Sequence[sp.Expr | sp.Poly],
    variables: Sequence[sp.Symbol],
    *,
    equational_constraints: Sequence[sp.Expr] = (),
    certify: bool = False,
) -> ReducedProjectionTower:
    """Build a Lazard-style reduced projection tower for diagnostics.

    Full Lazard valuation-aware lifting belongs in future work. This constructor
    therefore records a non-certified reduced tower unless certification is
    explicitly possible.
    """

    return build_reduced_proj_tower(
        polys,
        variables,
        theory="lazard",
        equational_constraints=equational_constraints,
        certify=certify,
    )
