"""Exact probability calculations for polynomial densities on semialgebraic supports."""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field

import sympy as sp

from ..decision import implies
from ..exact_arithmetic import compare_exact_reals
from ..normalization import normalize_bounds, normalize_formula, normalize_variables
from ..region_integrate import integrate_over_region

FormulaLike = sp.Expr | sp.logic.boolalg.Boolean | bool
_DEFAULT_DENSITY = sp.Integer(1)


@dataclass(frozen=True)
class PolynomialProbabilityResult:
    """Exact normalized probability of a semialgebraic event."""

    event: sp.Expr
    support: sp.Expr
    density: sp.Expr
    variables: tuple[sp.Symbol, ...]
    normalizing_mass: sp.Expr
    event_mass: sp.Expr
    probability: sp.Expr
    density_nonnegative: bool
    certified: bool = True
    method: str = "exact_semialgebraic_probability"
    diagnostics: Mapping[str, object] = field(default_factory=dict)


def polynomial_probability(
    event: FormulaLike,
    variables: Sequence[sp.Symbol | str],
    *,
    support: FormulaLike | None = None,
    density: sp.Expr = _DEFAULT_DENSITY,
    bounds: Sequence[tuple[sp.Symbol | str, object, object]]
    | Mapping[sp.Symbol | str, tuple[object, object]]
    | None = None,
) -> PolynomialProbabilityResult:
    """Compute an exact probability for a polynomial density.

    ``density`` need not be normalized.  semialg certifies that it is
    nonnegative on ``support``, integrates it over the support to obtain the
    normalizing mass, integrates over ``support & event`` for the event mass,
    and returns their exact ratio.  A positive finite normalizing mass is
    required.
    """

    event_formula = normalize_formula(event)
    support_formula = normalize_formula(support) if support is not None else sp.true
    density_expr = sp.expand(sp.sympify(density))
    context = sp.Tuple(event_formula, support_formula, density_expr)
    vars_ = normalize_variables(variables, context, append_context_symbols=False)
    if not vars_:
        raise ValueError("at least one random variable is required")
    bound_map = normalize_bounds(bounds, vars_)

    extras = (
        event_formula.free_symbols | support_formula.free_symbols | density_expr.free_symbols
    ) - set(vars_)
    if extras:
        names = ", ".join(sorted(symbol.name for symbol in extras))
        raise ValueError(
            "undeclared symbolic parameters are not supported by polynomial probability: " + names
        )
    try:
        sp.Poly(density_expr, *vars_)
    except sp.PolynomialError as exc:
        raise ValueError("probability density must be polynomial in the random variables") from exc

    bound_condition = (
        sp.And(*(sp.And(var >= lower, var <= upper) for var, (lower, upper) in bound_map.items()))
        if bound_map
        else sp.true
    )
    effective_support = sp.And(support_formula, bound_condition)

    nonnegative_check = implies(effective_support, density_expr >= 0, vars_, return_result=True)
    if not bool(nonnegative_check):
        raise ValueError("probability density must be nonnegative throughout the support")

    total_mass = integrate_over_region(
        density_expr,
        support_formula,
        vars_,
        bounds=bound_map,
    )
    if total_mass in (sp.oo, -sp.oo):
        raise ValueError("probability support must have a finite positive normalizing mass")
    relation = compare_exact_reals(sp.sympify(total_mass), sp.Integer(0))
    if relation <= 0:
        raise ValueError("probability support must have a finite positive normalizing mass")

    event_mass = integrate_over_region(
        density_expr,
        sp.And(support_formula, event_formula),
        vars_,
        bounds=bound_map,
    )
    probability = sp.simplify(event_mass / total_mass)
    return PolynomialProbabilityResult(
        event=event_formula,
        support=support_formula,
        density=density_expr,
        variables=vars_,
        normalizing_mass=sp.simplify(total_mass),
        event_mass=sp.simplify(event_mass),
        probability=probability,
        density_nonnegative=True,
        diagnostics={
            "bounds": bound_map,
            "density_certificate_method": getattr(nonnegative_check, "method", "implication"),
        },
    )


def geometric_probability(
    event: FormulaLike,
    variables: Sequence[sp.Symbol | str],
    *,
    support: FormulaLike,
    bounds: Sequence[tuple[sp.Symbol | str, object, object]]
    | Mapping[sp.Symbol | str, tuple[object, object]]
    | None = None,
) -> PolynomialProbabilityResult:
    """Return exact uniform geometric probability on a semialgebraic support."""

    return polynomial_probability(
        event,
        variables,
        support=support,
        density=sp.Integer(1),
        bounds=bounds,
    )


__all__ = ["PolynomialProbabilityResult", "geometric_probability", "polynomial_probability"]
