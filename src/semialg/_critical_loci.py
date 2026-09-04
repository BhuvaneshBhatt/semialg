"""Shared exact construction of projected positive-dimensional critical loci."""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass

import sympy as sp
from sympy.polys.polyerrors import PolynomialError

from .algebraic.groebner_utils import compute_groebner_basis
from .optimization_active_sets import jacobian_rank_equations, kkt_system, pruned_active_subsets
from .optimization_geometry import polynomial_locus_dimension

LocusKey = tuple[str, tuple[sp.Expr, ...]]


@dataclass(frozen=True)
class ProjectedCriticalLocus:
    """A distinct positive-dimensional locus in the original problem variables."""

    kind: str
    equations: tuple[sp.Expr, ...]
    dimension: int
    key: LocusKey


def project_kkt_locus(
    equations: Sequence[sp.Expr],
    multipliers: Sequence[sp.Symbol],
    variables: Sequence[sp.Symbol],
) -> tuple[sp.Expr, ...]:
    """Eliminate KKT multipliers, returning equations on original variables."""

    if not multipliers:
        expanded = tuple(sp.expand(eq) for eq in equations)
        return tuple(dict.fromkeys(eq for eq in expanded if eq != 0))
    all_vars = (*multipliers, *variables)
    try:
        basis = compute_groebner_basis(equations, all_vars, order="lex", domain=sp.QQ)
    except (PolynomialError, ValueError, TypeError):
        return ()
    multiplier_set = set(multipliers)
    projected = [
        sp.expand(poly.as_expr())
        for poly in basis.polys
        if not (poly.as_expr().free_symbols & multiplier_set)
    ]
    return tuple(dict.fromkeys(expr for expr in projected if expr != 0))


def _locus_key(kind: str, equations: Sequence[sp.Expr]) -> LocusKey:
    canonical = tuple(
        sorted(
            {sp.factor(eq) for eq in equations},
            key=sp.default_sort_key,
        )
    )
    return kind, canonical


def iter_positive_dimensional_critical_loci(
    objective: sp.Expr,
    condition: sp.Expr,
    variables: tuple[sp.Symbol, ...],
    equalities: tuple[sp.Expr, ...],
    inequalities: tuple[sp.Expr, ...],
    *,
    include_singular: bool = True,
) -> tuple[ProjectedCriticalLocus, ...]:
    """Return distinct projected positive-dimensional KKT and singular loci.

    Active-set construction, multiplier projection, dimension checks, and
    structural deduplication are centralized here so optimization and geometry
    queries cannot drift to different critical-locus semantics.
    """

    out: list[ProjectedCriticalLocus] = []
    seen: set[LocusKey] = set()
    for subset in pruned_active_subsets(equalities, inequalities, variables, condition):
        active = tuple(equalities) + tuple(subset)
        equations, multipliers = kkt_system(objective, variables, active)
        full_dimension = polynomial_locus_dimension(equations, (*variables, *multipliers))
        if full_dimension is not None and full_dimension > 0:
            projected = project_kkt_locus(equations, multipliers, variables)
            if projected:
                dimension = polynomial_locus_dimension(projected, variables)
                if dimension is not None and 0 < dimension < len(variables):
                    key = _locus_key("kkt", projected)
                    if key not in seen:
                        seen.add(key)
                        out.append(ProjectedCriticalLocus("kkt", projected, dimension, key))

        if include_singular and active:
            minors = jacobian_rank_equations(active, variables)
            if minors:
                singular_equations = tuple(active) + tuple(minors)
                dimension = polynomial_locus_dimension(singular_equations, variables)
                if dimension is not None and 0 < dimension < len(variables):
                    key = _locus_key("singular", singular_equations)
                    if key not in seen:
                        seen.add(key)
                        out.append(
                            ProjectedCriticalLocus(
                                "singular", tuple(singular_equations), dimension, key
                            )
                        )
    return tuple(out)


__all__ = [
    "LocusKey",
    "ProjectedCriticalLocus",
    "iter_positive_dimensional_critical_loci",
    "project_kkt_locus",
]
