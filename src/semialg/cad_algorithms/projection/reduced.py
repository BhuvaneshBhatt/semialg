from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field
from typing import Literal

import sympy as sp

from ...ec.scoring import rank_eq_cons
from ..polynomial_utils import (
    _cad_input_polynomials,
    derivative_resultant_expression,
    discriminant_expression,
    lower_polynomial,
    resultant_expression,
    univariate_coefficients,
    univariate_leading_coefficient,
)
from .collins import ProjectionLevel, ProjectionTower, collins_projection_step, squarefree_basis
from .ec_hypotheses import ECLevelHypotheses, propagated_equations

ReducedTheory = Literal["mccallum", "lazard", "tticad"]


@dataclass(frozen=True)
class ProjectionValidity:
    """Auditable status for a reduced projection attempt.

    Reduced CAD theories are useful only when their side conditions have been
    checked. Safe drivers accept a reduced tower only after degree/order
    induction or a supported reference certificate proves the requested sign,
    truth or truth-table invariant.
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
    except (sp.PolynomialError, TypeError, ValueError, NotImplementedError, ArithmeticError):
        return 0


def _leading_coefficient(
    poly: sp.Poly, var: sp.Symbol, lower_gens: Sequence[sp.Symbol]
) -> list[sp.Poly]:
    if not lower_gens:
        return []
    try:
        item = lower_polynomial(univariate_leading_coefficient(poly, var), lower_gens)
    except (sp.PolynomialError, ValueError, TypeError):
        return []
    return [] if item is None else [item]


def _discriminant(poly: sp.Poly, var: sp.Symbol, lower_gens: Sequence[sp.Symbol]) -> list[sp.Poly]:
    if not lower_gens or _degree(poly, var) <= 1:
        return []
    item = lower_polynomial(discriminant_expression(poly, var), lower_gens)
    return [] if item is None else [item]


def _resultant(
    left: sp.Poly, right: sp.Poly, var: sp.Symbol, lower_gens: Sequence[sp.Symbol]
) -> list[sp.Poly]:
    if not lower_gens:
        return []
    item = lower_polynomial(resultant_expression(left, right, var), lower_gens)
    return [] if item is None else [item]


def subres_discrim(poly: sp.Poly, var: sp.Symbol, lower_gens: Sequence[sp.Symbol]) -> list[sp.Poly]:
    if not lower_gens or _degree(poly, var) <= 1:
        return []
    item = lower_polynomial(derivative_resultant_expression(poly, var), lower_gens)
    return [] if item is None else [item]


def _candidate_ec(
    active: Sequence[sp.Poly], equational_constraints: Sequence[sp.Expr]
) -> sp.Poly | None:
    """Choose the simplest active polynomial matching a declared EC.

    Ranking by degree and sparsity makes the choice deterministic and usually
    reduces discriminant/resultant growth compared with basis-order selection.
    """
    matches: list[tuple[sp.Expr, sp.Poly]] = []
    expanded_ecs = tuple(sp.expand(item) for item in equational_constraints)
    for poly in active:
        expr = sp.expand(poly.as_expr())
        for ec in expanded_ecs:
            if expr == ec or expr == -ec:
                matches.append((ec, poly))
                break
    if not matches:
        return None
    ranked = rank_eq_cons([expr for expr, _ in matches], tuple(active[0].gens) if active else ())
    best = ranked[0].expr
    for expr, poly in matches:
        if expr == best or expr == -best:
            return poly
    return matches[0][1]


def reduced_projection_step(
    polys: Sequence[sp.Poly],
    var: sp.Symbol,
    lower_gens: Sequence[sp.Symbol],
    *,
    theory: ReducedTheory,
    equational_constraints: Sequence[sp.Expr] = (),
    _normalized_basis: tuple[sp.Poly, ...] | None = None,
    preserve_order: bool = False,
) -> tuple[sp.Poly, ...]:
    """Construct a diagnostic reduced-projection step.

    References: [McCallum1988], [Brown2001], [Lazard1994], and [BDEMW2016]
    in ``docs/references.md``. The implementation is conservative:
    reduced towers are accelerators whose safe drivers certify or fall back to
    complete Collins-style projection/lifting.

    This is not advertised as complete by itself. It implements
    the standard EC-shaped subset used by McCallum/Lazard/TTICAD-style paths:
    leading coefficient/discriminant of a designated EC plus resultants between
    that EC and the other active polynomials. Without a designated EC it falls
    back to full coefficient projection for that level. ``preserve_order``
    adds all active leading coefficients/discriminants, pairwise resultants
    and designated EC coefficients for intermediate recursive EC lifts.
    """

    basis = _normalized_basis if _normalized_basis is not None else squarefree_basis(polys)
    active = [poly for poly in basis if _degree(poly, var) > 0]
    inactive = [poly for poly in basis if _degree(poly, var) == 0]
    projected: list[sp.Poly] = []
    if lower_gens:
        projected.extend(
            item
            for poly in inactive
            if (item := lower_polynomial(poly.as_expr(), lower_gens)) is not None
        )
    designated = _candidate_ec(active, equational_constraints)
    if designated is None:
        return collins_projection_step(basis, var, lower_gens)
    projected.extend(_leading_coefficient(designated, var, lower_gens))
    projected.extend(_discriminant(designated, var, lower_gens))
    if theory in {"mccallum", "lazard", "tticad"}:
        projected.extend(subres_discrim(designated, var, lower_gens))
    for other in active:
        if sp.expand(other.as_expr()) != sp.expand(designated.as_expr()):
            projected.extend(_resultant(designated, other, var, lower_gens))
    if preserve_order:
        # Subsequent lifts need order invariance of non-EC polynomials on
        # EC sections. Their discriminants and degree guards are essential;
        # repeating the single-EC sign operator does not supply them.
        for other in active:
            projected.extend(_leading_coefficient(other, var, lower_gens))
            projected.extend(_discriminant(other, var, lower_gens))
        others = [poly for poly in active if poly != designated]
        for index, left in enumerate(others):
            for right in others[index + 1 :]:
                projected.extend(_resultant(left, right, var, lower_gens))
        for coefficient in univariate_coefficients(designated, var):
            item = lower_polynomial(coefficient, lower_gens)
            if item is not None:
                projected.append(item)
    return squarefree_basis(projected)


def build_reduced_proj_tower(
    polys: Sequence[sp.Expr | sp.Poly],
    variables: Sequence[sp.Symbol],
    *,
    theory: ReducedTheory,
    equational_constraints: Sequence[sp.Expr] = (),
    certify: bool = False,
    recursive_ec: bool = False,
) -> ReducedProjectionTower:
    """Build a reduced projection tower with explicit safe-use metadata.

    The tower itself is only a projection object. Formula drivers replay EC
    implications and intermediate degree/order conditions; other safe drivers
    use their supported reference certificates. Failed hypotheses select an
    ambient sign-invariant fallback.
    """

    vars_tuple = tuple(variables)
    current = _as_polys(_cad_input_polynomials(polys, vars_tuple), vars_tuple)
    implications = propagated_equations(equational_constraints, vars_tuple) if recursive_ec else ()
    declared = (
        tuple(item.polynomial for item in implications)
        if recursive_ec
        else tuple(equational_constraints)
    )
    plans = []
    levels: dict[int, tuple[sp.Poly, ...]] = {len(vars_tuple): current}
    reduced_levels: list[int] = []
    for level in range(len(vars_tuple), 1, -1):
        var = vars_tuple[level - 1]
        normalized = squarefree_basis(current)
        active = [poly for poly in normalized if _degree(poly, var) > 0]
        designated = _candidate_ec(active, declared)
        uses_reduced_operator = designated is not None
        preserve_order = recursive_ec and level < len(vars_tuple)
        reduced = reduced_projection_step(
            current,
            var,
            vars_tuple[: level - 1],
            theory=theory,
            equational_constraints=declared,
            _normalized_basis=normalized,
            preserve_order=preserve_order,
        )
        if recursive_ec:
            degree_polys = (
                active if preserve_order else ([] if designated is None else [designated])
            )
            plans.append(
                ECLevelHypotheses(
                    level,
                    None if designated is None else designated.as_expr(),
                    "full" if designated is None else "ec-order" if preserve_order else "ec-sign",
                    tuple(univariate_leading_coefficient(poly, var) for poly in degree_polys)
                    if designated is not None
                    else (),
                    tuple(poly.as_expr() for poly in reduced),
                )
            )
        if uses_reduced_operator:
            reduced_levels.append(level)
        current = reduced
        levels[level - 1] = current
    if recursive_ec and vars_tuple:
        designated = _candidate_ec(
            [poly for poly in current if poly.degree(vars_tuple[0]) > 0], declared
        )
        plans.append(
            ECLevelHypotheses(1, None if designated is None else designated.as_expr(), "base")
        )
    has_reduction = bool(reduced_levels) or (
        recursive_ec and any(plan.equational_constraint is not None for plan in plans)
    )
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
            "complete": bool(certify and not has_reduction),
            **(
                {"ec_lifting": True, "ec_levels": tuple(plans), "ec_implications": implications}
                if recursive_ec
                else {}
            ),
        },
    )
    valid = bool(certify and not has_reduction)
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
