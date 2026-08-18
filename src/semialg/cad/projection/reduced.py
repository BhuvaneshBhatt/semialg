from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field
from typing import Literal

import sympy as sp

from .collins import ProjectionLevel, ProjectionTower, collins_projection_step, squarefree_basis

ReducedTheory = Literal["mccallum", "lazard", "tticad"]


@dataclass(frozen=True)
class ProjectionValidity:
    """Auditable status for a reduced projection attempt.

    Reduced CAD theories are useful only when their side conditions have been
    checked. The reduced projection path makes this status proof-carrying: safe drivers may
    actively lift a reduced tower, but they accept it only after an explicit
    certificate proves sign, truth, or truth-table invariance against a complete
    Collins refinement.
    """

    theory: ReducedTheory
    valid: bool
    complete_if_used: bool
    reason: str
    checked_conditions: tuple[str, ...] = ()
    failed_conditions: tuple[str, ...] = ()
    fallback_backend: str | None = "collins-complete"
    details: Mapping[str, object] = field(default_factory=dict)


@dataclass(frozen=True)
class ReducedProjectionTower:
    """Reduced projection data plus a conservative validity report."""

    requested_theory: ReducedTheory
    tower: ProjectionTower
    validity: ProjectionValidity


def _as_polys(
    polys: Sequence[sp.Expr | sp.Poly], variables: Sequence[sp.Symbol]
) -> tuple[sp.Poly, ...]:
    return squarefree_basis(
        poly if isinstance(poly, sp.Poly) else sp.Poly(sp.expand(poly), *variables)
        for poly in polys
    )


def _degree(poly: sp.Poly, var: sp.Symbol) -> int:
    try:
        return int(poly.degree(var))
    except Exception:
        return 0


def _lower_poly(expr: sp.Expr, lower_gens: Sequence[sp.Symbol]) -> sp.Poly | None:
    expr = sp.expand(expr)
    if expr == 0 or not lower_gens:
        return None
    try:
        poly = sp.Poly(expr, *lower_gens)
    except Exception:
        return None
    basis = squarefree_basis([poly])
    return basis[0] if basis else None


def _coefficients(poly: sp.Poly, var: sp.Symbol, lower_gens: Sequence[sp.Symbol]) -> list[sp.Poly]:
    if not lower_gens:
        return []
    out: list[sp.Poly] = []
    try:
        coeffs = sp.Poly(poly.as_expr(), var).all_coeffs()
    except Exception:
        return []
    for coeff in coeffs:
        item = _lower_poly(coeff, lower_gens)
        if item is not None:
            out.append(item)
    return out


def _leading_coefficient(
    poly: sp.Poly, var: sp.Symbol, lower_gens: Sequence[sp.Symbol]
) -> list[sp.Poly]:
    if not lower_gens:
        return []
    try:
        return [
            item
            for item in [_lower_poly(sp.Poly(poly.as_expr(), var).LC(), lower_gens)]
            if item is not None
        ]
    except Exception:
        return []


def _discriminant(poly: sp.Poly, var: sp.Symbol, lower_gens: Sequence[sp.Symbol]) -> list[sp.Poly]:
    if not lower_gens or _degree(poly, var) <= 1:
        return []
    return [
        item
        for item in [_lower_poly(sp.discriminant(poly.as_expr(), var), lower_gens)]
        if item is not None
    ]


def _resultant(
    left: sp.Poly, right: sp.Poly, var: sp.Symbol, lower_gens: Sequence[sp.Symbol]
) -> list[sp.Poly]:
    if not lower_gens:
        return []
    return [
        item
        for item in [_lower_poly(sp.resultant(left.as_expr(), right.as_expr(), var), lower_gens)]
        if item is not None
    ]


def subres_discrim(poly: sp.Poly, var: sp.Symbol, lower_gens: Sequence[sp.Symbol]) -> list[sp.Poly]:
    if not lower_gens or _degree(poly, var) <= 1:
        return []
    deriv = sp.diff(poly.as_expr(), var)
    return [
        item
        for item in [_lower_poly(sp.resultant(poly.as_expr(), deriv, var), lower_gens)]
        if item is not None
    ]


def _candidate_ec(
    active: Sequence[sp.Poly], equational_constraints: Sequence[sp.Expr]
) -> sp.Poly | None:
    expanded_ecs = [sp.expand(item) for item in equational_constraints]
    for poly in active:
        expr = sp.expand(poly.as_expr())
        if any(sp.expand(expr - ec) == 0 or sp.expand(expr + ec) == 0 for ec in expanded_ecs):
            return poly
    return None


def reduced_projection_step(
    polys: Sequence[sp.Poly],
    var: sp.Symbol,
    lower_gens: Sequence[sp.Symbol],
    *,
    theory: ReducedTheory,
    equational_constraints: Sequence[sp.Expr] = (),
) -> tuple[sp.Poly, ...]:
    """Construct a diagnostic reduced-projection step.

    This is deliberately not advertised as complete by itself. It implements
    the standard EC-shaped subset used by McCallum/Lazard/TTICAD-style paths:
    leading coefficient/discriminant of a designated EC plus resultants between
    that EC and the other active polynomials. Without a designated EC it falls
    back to the full Collins projection step for that level.
    """

    basis = squarefree_basis(polys)
    active = [poly for poly in basis if _degree(poly, var) > 0]
    inactive = [poly for poly in basis if _degree(poly, var) == 0]
    projected: list[sp.Poly] = []
    if lower_gens:
        projected.extend(sp.Poly(poly.as_expr(), *lower_gens) for poly in inactive)
    designated = _candidate_ec(active, equational_constraints)
    if designated is None:
        return collins_projection_step(basis, var, lower_gens)
    projected.extend(_leading_coefficient(designated, var, lower_gens))
    projected.extend(_discriminant(designated, var, lower_gens))
    if theory in {"mccallum", "lazard", "tticad"}:
        projected.extend(subres_discrim(designated, var, lower_gens))
    for other in active:
        if sp.expand(other.as_expr() - designated.as_expr()) != 0:
            projected.extend(_resultant(designated, other, var, lower_gens))
    return squarefree_basis(projected)


def build_reduced_proj_tower(
    polys: Sequence[sp.Expr | sp.Poly],
    variables: Sequence[sp.Symbol],
    *,
    theory: ReducedTheory,
    equational_constraints: Sequence[sp.Expr] = (),
    certify: bool = False,
) -> ReducedProjectionTower:
    """Build a reduced projection tower with explicit safe-use metadata.

    The tower itself is only a projection object. The active safe drivers in
    :mod:`semialg.cad.reduced` and :mod:`semialg.tticad.safe` decide whether it
    can be used by constructing a reduced CAD and certifying it against a full
    Collins refinement.
    """

    vars_tuple = tuple(variables)
    current = _as_polys(polys, vars_tuple)
    levels: dict[int, tuple[sp.Poly, ...]] = {len(vars_tuple): current}
    reduced_levels: list[int] = []
    for level in range(len(vars_tuple), 1, -1):
        reduced = reduced_projection_step(
            current,
            vars_tuple[level - 1],
            vars_tuple[: level - 1],
            theory=theory,
            equational_constraints=equational_constraints,
        )
        full = collins_projection_step(current, vars_tuple[level - 1], vars_tuple[: level - 1])
        if len(reduced) < len(full):
            reduced_levels.append(level)
        current = reduced
        levels[level - 1] = current
    tower = ProjectionTower(
        variables=vars_tuple,
        levels=tuple(
            ProjectionLevel(i, vars_tuple[i - 1], levels.get(i, ()))
            for i in range(1, len(vars_tuple) + 1)
        ),
        original_polynomials=levels[len(vars_tuple)],
        metadata={
            "projection": theory,
            "reduced_levels": tuple(reduced_levels),
            "complete": bool(certify and not reduced_levels),
        },
    )
    valid = bool(certify and not reduced_levels)
    validity = ProjectionValidity(
        theory=theory,
        valid=valid,
        complete_if_used=valid,
        reason=(
            "no reduced steps were used; equivalent to full Collins projection"
            if valid
            else "projection tower constructed; active safe driver must certify reduced lifting before use"
        ),
        checked_conditions=("projection construction",),
        failed_conditions=() if valid else ("reduced lifting certificate pending",),
        fallback_backend=None if valid else "collins-complete",
        details={
            "reduced_levels": tuple(reduced_levels),
            "poly_count_by_level": {k: len(v) for k, v in levels.items()},
        },
    )
    return ReducedProjectionTower(requested_theory=theory, tower=tower, validity=validity)
