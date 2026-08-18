from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass

import sympy as sp

from .formula import to_sympy
from .region_integrate import integrate_over_region


@dataclass(frozen=True)
class RegionMomentResult:
    """Moment integral over a semialgebraic region.

    ``value`` is the raw moment, not normalized by the region measure. Use
    ``region_centroid`` or ``region_covariance`` for normalized first and
    second central moments.
    """

    value: sp.Expr
    integrand: sp.Expr
    condition: sp.Expr
    variables: tuple[sp.Symbol, ...]
    powers: tuple[int, ...] | None
    method: str
    diagnostics: Mapping[str, object] | None = None
    exact: bool = True


@dataclass(frozen=True)
class RegionCentroidResult:
    """Centroid of a semialgebraic region."""

    centroid: Mapping[sp.Symbol, sp.Expr]
    measure: sp.Expr
    condition: sp.Expr
    variables: tuple[sp.Symbol, ...]
    method: str
    diagnostics: Mapping[str, object] | None = None
    exact: bool = True


@dataclass(frozen=True)
class RegionCovarianceResult:
    """Covariance matrix for the normalized uniform measure on a region."""

    covariance: sp.Matrix
    centroid: Mapping[sp.Symbol, sp.Expr]
    measure: sp.Expr
    condition: sp.Expr
    variables: tuple[sp.Symbol, ...]
    method: str
    diagnostics: Mapping[str, object] | None = None
    exact: bool = True


def _normalize_variables(variables: Sequence[sp.Symbol | str]) -> tuple[sp.Symbol, ...]:
    out: list[sp.Symbol] = []
    seen: set[sp.Symbol] = set()
    for var in variables:
        sym = sp.Symbol(var, real=True) if isinstance(var, str) else var
        if sym not in seen:
            out.append(sym)
            seen.add(sym)
    return tuple(out)


def _normalize_formula(condition: object) -> sp.Expr:
    if isinstance(condition, (sp.Basic, sp.logic.boolalg.Boolean)):
        return condition  # type: ignore[return-value]
    return to_sympy(condition)  # type: ignore[arg-type]


def _monomial_from_powers(variables: Sequence[sp.Symbol], powers: Sequence[int]) -> sp.Expr:
    if len(powers) != len(variables):
        raise ValueError("powers must have the same length as variables")
    if any(power < 0 for power in powers):
        raise ValueError("moment powers must be nonnegative integers")
    monomial = sp.Integer(1)
    for var, power in zip(variables, powers, strict=True):
        monomial *= var ** int(power)
    return sp.expand(monomial)


def _moment_integrand(
    variables: Sequence[sp.Symbol],
    powers: Sequence[int] | sp.Expr | None,
    integrand: object | None,
) -> tuple[sp.Expr, tuple[int, ...] | None]:
    if powers is not None and integrand is not None:
        raise ValueError("pass either powers or integrand, not both")
    if integrand is not None:
        return sp.sympify(integrand), None
    if powers is None:
        powers_tuple = tuple(0 for _ in variables)
        return sp.Integer(1), powers_tuple
    if isinstance(powers, sp.Expr):
        return sp.sympify(powers), None
    powers_tuple = tuple(int(power) for power in powers)
    return _monomial_from_powers(variables, powers_tuple), powers_tuple


def _integral_value_and_exact(value: sp.Expr | object) -> tuple[sp.Expr, bool]:
    expr = sp.sympify(value)
    return expr, not bool(expr.has(sp.Float))


def _validate_finite_nonzero_measure(measure: sp.Expr) -> None:
    if sp.simplify(measure) == 0:
        raise ValueError("region has zero measure, so normalized moments are undefined")
    if measure in (sp.oo, -sp.oo, sp.zoo) or measure.has(sp.oo, -sp.oo, sp.zoo):
        raise ValueError(
            "region has infinite or indeterminate measure, so normalized moments are undefined"
        )


def region_moment(
    condition: object,
    variables: Sequence[sp.Symbol | str],
    powers: Sequence[int] | sp.Expr | None = None,
    *,
    integrand: object | None = None,
    bounds: Sequence[tuple[sp.Symbol | str, object, object]]
    | Mapping[sp.Symbol | str, tuple[object, object]]
    | None = None,
    method: str = "symbolic",
    precision: int = 50,
    measure_dimension: object = "ambient",
    return_result: bool = False,
) -> sp.Expr | RegionMomentResult:
    """Return a raw moment integral over a semialgebraic region.

    ``powers`` gives a monomial moment. For variables ``(x, y)`` and powers
    ``(2, 1)``, the integrated moment is ``x**2*y``. Alternatively, callers
    may pass an explicit ``integrand`` for a general raw moment.
    """

    vars_ = _normalize_variables(variables)
    formula = _normalize_formula(condition)
    moment_integrand, power_tuple = _moment_integrand(vars_, powers, integrand)
    value = integrate_over_region(
        moment_integrand,
        formula,
        vars_,
        bounds=bounds,
        method=method,
        precision=precision,
        measure_dimension=measure_dimension,
    )
    expr, exact = _integral_value_and_exact(value)
    result = RegionMomentResult(
        value=sp.simplify(expr) if exact else expr,
        integrand=moment_integrand,
        condition=formula,
        variables=vars_,
        powers=power_tuple,
        method=method,
        diagnostics={
            "bounds": bounds,
            "precision": precision,
            "measure_dimension": measure_dimension,
        },
        exact=exact,
    )
    return result if return_result else result.value


def region_centroid(
    condition: object,
    variables: Sequence[sp.Symbol | str],
    *,
    bounds: Sequence[tuple[sp.Symbol | str, object, object]]
    | Mapping[sp.Symbol | str, tuple[object, object]]
    | None = None,
    method: str = "symbolic",
    precision: int = 50,
    measure_dimension: object = "ambient",
    return_result: bool = False,
) -> Mapping[sp.Symbol, sp.Expr] | RegionCentroidResult:
    """Return the centroid of a finite-measure semialgebraic region."""

    vars_ = _normalize_variables(variables)
    formula = _normalize_formula(condition)
    measure = integrate_over_region(
        sp.Integer(1),
        formula,
        vars_,
        bounds=bounds,
        method=method,
        precision=precision,
        measure_dimension=measure_dimension,
    )
    measure_expr, measure_exact = _integral_value_and_exact(measure)
    _validate_finite_nonzero_measure(measure_expr)

    centroid: dict[sp.Symbol, sp.Expr] = {}
    exact = measure_exact
    for var in vars_:
        raw = integrate_over_region(
            var,
            formula,
            vars_,
            bounds=bounds,
            method=method,
            precision=precision,
            measure_dimension=measure_dimension,
        )
        raw_expr, raw_exact = _integral_value_and_exact(raw)
        exact = exact and raw_exact
        centroid[var] = (
            sp.simplify(raw_expr / measure_expr)
            if exact
            else sp.N(raw_expr / measure_expr, precision)
        )

    result = RegionCentroidResult(
        centroid=centroid,
        measure=sp.simplify(measure_expr) if measure_exact else measure_expr,
        condition=formula,
        variables=vars_,
        method=method,
        diagnostics={
            "bounds": bounds,
            "precision": precision,
            "measure_dimension": measure_dimension,
        },
        exact=exact,
    )
    return result if return_result else result.centroid


def region_covariance(
    condition: object,
    variables: Sequence[sp.Symbol | str],
    *,
    bounds: Sequence[tuple[sp.Symbol | str, object, object]]
    | Mapping[sp.Symbol | str, tuple[object, object]]
    | None = None,
    method: str = "symbolic",
    precision: int = 50,
    measure_dimension: object = "ambient",
    return_result: bool = False,
) -> sp.Matrix | RegionCovarianceResult:
    """Return the covariance matrix of the uniform measure on a region."""

    vars_ = _normalize_variables(variables)
    formula = _normalize_formula(condition)
    centroid_result = region_centroid(
        formula,
        vars_,
        bounds=bounds,
        method=method,
        precision=precision,
        measure_dimension=measure_dimension,
        return_result=True,
    )
    assert isinstance(centroid_result, RegionCentroidResult)
    measure = centroid_result.measure
    _validate_finite_nonzero_measure(measure)

    entries: list[list[sp.Expr]] = []
    exact = centroid_result.exact
    for _i, vi in enumerate(vars_):
        row: list[sp.Expr] = []
        for _j, vj in enumerate(vars_):
            raw = integrate_over_region(
                vi * vj,
                formula,
                vars_,
                bounds=bounds,
                method=method,
                precision=precision,
                measure_dimension=measure_dimension,
            )
            raw_expr, raw_exact = _integral_value_and_exact(raw)
            exact = exact and raw_exact
            value = raw_expr / measure - centroid_result.centroid[vi] * centroid_result.centroid[vj]
            row.append(sp.simplify(value) if exact else sp.N(value, precision))
        entries.append(row)
    matrix = sp.Matrix(entries)
    result = RegionCovarianceResult(
        covariance=matrix,
        centroid=centroid_result.centroid,
        measure=measure,
        condition=formula,
        variables=vars_,
        method=method,
        diagnostics={
            "bounds": bounds,
            "precision": precision,
            "measure_dimension": measure_dimension,
        },
        exact=exact,
    )
    return result if return_result else result.covariance


__all__ = [
    "RegionMomentResult",
    "RegionCentroidResult",
    "RegionCovarianceResult",
    "region_moment",
    "region_centroid",
    "region_covariance",
]
