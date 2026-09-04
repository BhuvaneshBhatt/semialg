from __future__ import annotations

import random as _random
from collections.abc import Iterable, Mapping, Sequence
from itertools import product

import sympy as sp
from sympy.logic.boolalg import Boolean

from .algebraic import (
    AlgebraicRoot,
    RationalSample,
    RationalUnivariatePoint,
    Sample,
    sample_to_expr,
    sign_at_sample,
    sign_of_algebraic_expression,
)
from .algebraic.rational_univariate import solve_formula_with_rur
from .dimension_validation import require_point_dimension, require_same_length, zip_equal
from .errors import DimensionMismatchError, SemialgStrategyFailure
from .formula import to_sympy
from .instances.real_fallbacks import satisfies_formula
from .normalization import normalize_sampling_variables
from .solve import find_instance
from .structural_keys import point_key

_RECOVERABLE_ERRORS = (
    ArithmeticError,
    TypeError,
    ValueError,
    NotImplementedError,
    sp.PolynomialError,
    SemialgStrategyFailure,
)


def _as_expr_value(value: object) -> sp.Expr:
    if isinstance(value, (RationalSample, AlgebraicRoot)):
        return sample_to_expr(value)
    return sp.sympify(value)


def _is_rur_point(point: object) -> bool:
    return isinstance(point, RationalUnivariatePoint)


def _sign_at_rur_point(
    expr: sp.Expr, point: RationalUnivariatePoint, variables: Sequence[sp.Symbol]
) -> int:
    """Evaluate ``expr`` at a RUR point via univariate reduction.

    Directly substituting algebraic coordinate expressions can be expensive and
    can lose the useful certificate carried by the RUR. Instead we substitute
    the coordinate polynomials in the RUR parameter, reduce modulo the defining
    polynomial, and decide the sign of the resulting algebraic expression at the
    represented parameter root.
    """

    representation = point.representation
    if tuple(variables) != tuple(representation.variables):
        missing = [var for var in variables if var not in representation.variables]
        if missing:
            raise DimensionMismatchError(f"RUR point does not contain assignments for {missing!r}")
    t = representation.parameter
    coord_polys = representation.normalized_coordinate_polynomials()
    coord_map = dict(
        zip_equal(
            representation.variables,
            (poly.as_expr() for poly in coord_polys),
            context="RUR coordinate polynomials",
        )
    )
    numerator, denominator = sp.fraction(sp.cancel(sp.sympify(expr).subs(coord_map)))
    q = representation.defining_polynomial
    numerator_poly = sp.Poly(numerator, t, domain=sp.QQ).rem(q)
    denominator_poly = sp.Poly(denominator, t, domain=sp.QQ).rem(q)
    numerator_sign = sign_of_algebraic_expression(numerator_poly.as_expr().subs(t, point.root))
    if numerator_sign == 0:
        return 0
    denominator_sign = sign_of_algebraic_expression(denominator_poly.as_expr().subs(t, point.root))
    if denominator_sign == 0:
        raise ValueError(f"expression {sp.sstr(expr)} is undefined at the supplied RUR point")
    return numerator_sign * denominator_sign


def _ordered_values(
    point: Mapping[sp.Symbol, object] | Sequence[object],
    variables: Sequence[sp.Symbol],
) -> tuple[sp.Expr, ...]:
    if isinstance(point, Mapping):
        missing = [var for var in variables if var not in point]
        if missing:
            raise DimensionMismatchError(f"point is missing assignments for {missing!r}")
        return tuple(_as_expr_value(point[var]) for var in variables)
    values = tuple(_as_expr_value(value) for value in point)
    require_point_dimension(values, variables, context="sampling point")
    return values


def sign_at(
    poly: sp.Poly | sp.Expr,
    point: Mapping[sp.Symbol, object] | Sequence[object] | RationalUnivariatePoint,
    *,
    variables: Sequence[sp.Symbol | str] | None = None,
    exact: bool = True,
    numeric_precision: int = 120,
) -> int:
    """Return the sign of a polynomial/expression at a point.

    The exact path handles ordinary rational SymPy values, semialg
    ``RationalSample``/``AlgebraicRoot`` objects, exact SymPy algebraic
    expressions such as ``sqrt(2)``/``RootOf``, and ``RationalUnivariatePoint``
    objects. Numeric fallback is available only when ``exact=False``.
    """

    if numeric_precision < 1:
        raise ValueError("numeric_precision must be positive")
    expr = poly.as_expr() if isinstance(poly, sp.Poly) else sp.sympify(poly)
    vars_ = normalize_sampling_variables(variables, sp.sympify(expr))

    if _is_rur_point(point):
        return _sign_at_rur_point(expr, point, vars_)

    values = _ordered_values(point, vars_)
    samples: list[Sample] = []
    all_semialg_samples = True
    for value in values:
        if isinstance(value, (RationalSample, AlgebraicRoot)):
            samples.append(value)
        elif getattr(value, "is_Rational", False):
            samples.append(RationalSample(sp.Rational(value)))
        else:
            all_semialg_samples = False
            break
    if all_semialg_samples and len(samples) == len(values):
        try:
            return sign_at_sample(expr, samples)
        except _RECOVERABLE_ERRORS:
            pass

    substituted = sp.cancel(
        expr.subs(dict(zip_equal(vars_, values, context="sampling substitution")))
    )
    try:
        return sign_of_algebraic_expression(substituted)
    except _RECOVERABLE_ERRORS:
        if exact:
            raise ValueError(
                f"could not determine exact sign of {sp.sstr(expr)} at {point!r}"
            ) from None

    numeric = sp.N(substituted, numeric_precision)
    if numeric == 0:
        return 0
    if numeric > 0:
        return 1
    if numeric < 0:
        return -1
    raise ValueError(f"could not determine sign of {sp.sstr(expr)} at {point!r}")


def sign_vector(
    polys: Iterable[sp.Poly | sp.Expr],
    point: Mapping[sp.Symbol, object] | Sequence[object] | RationalUnivariatePoint,
    *,
    variables: Sequence[sp.Symbol | str] | None = None,
    exact: bool = True,
    as_dict: bool = False,
    numeric_precision: int = 120,
) -> tuple[int, ...] | dict[sp.Expr, int]:
    """Return the signs of ``polys`` at ``point`` in input order.

    Set ``as_dict=True`` to receive a mapping from each input expression to its
    sign. This returns a tuple by default while supporting
    more inspectable CAD/debugging workflows.
    """

    exprs = tuple(
        poly.as_expr() if isinstance(poly, sp.Poly) else sp.sympify(poly) for poly in polys
    )
    vars_ = normalize_sampling_variables(variables, *(sp.sympify(expr) for expr in exprs))
    signs = tuple(
        sign_at(expr, point, variables=vars_, exact=exact, numeric_precision=numeric_precision)
        for expr in exprs
    )
    if as_dict:
        return dict(zip_equal(exprs, signs, context="sign vector"))
    return signs


def _small_rational_grid(variables: Sequence[sp.Symbol], radius: int = 3):
    values = [sp.Integer(0)]
    for k in range(1, radius + 1):
        values.extend(
            [sp.Integer(k), sp.Integer(-k), sp.Rational(1, k + 1), -sp.Rational(1, k + 1)]
        )
    for coords in product(values, repeat=len(variables)):
        yield dict(zip_equal(variables, coords, context="sampling grid point"))


def _bounds_for_variables(
    variables: Sequence[sp.Symbol],
    bounds: Sequence[tuple[object, object]] | Mapping[sp.Symbol, tuple[object, object]] | None,
    *,
    default_radius: int = 3,
) -> tuple[tuple[sp.Expr, sp.Expr], ...]:
    if bounds is None:
        return tuple((sp.Integer(-default_radius), sp.Integer(default_radius)) for _ in variables)
    if isinstance(bounds, Mapping):
        out = []
        for var in variables:
            if var not in bounds:
                out.append((sp.Integer(-default_radius), sp.Integer(default_radius)))
            else:
                lo, hi = bounds[var]
                out.append((sp.sympify(lo), sp.sympify(hi)))
        return tuple(out)
    raw = tuple(bounds)
    require_same_length(raw, variables, context="sampling bounds", names=("bounds", "variables"))
    return tuple((sp.sympify(lo), sp.sympify(hi)) for lo, hi in raw)


def _grid_values_for_interval(
    lo: sp.Expr, hi: sp.Expr, *, resolution: int, exact: bool, numeric_precision: int
) -> tuple[sp.Expr, ...]:
    def finalize(value: sp.Expr) -> sp.Expr:
        return sp.simplify(value) if exact else sp.Float(value, numeric_precision)

    if resolution <= 1:
        return (finalize((lo + hi) / 2),)
    return tuple(
        finalize(lo + (hi - lo) * sp.Rational(k, resolution - 1)) for k in range(resolution)
    )


def _bounded_rational_grid(
    variables: Sequence[sp.Symbol],
    bounds: Sequence[tuple[object, object]] | Mapping[sp.Symbol, tuple[object, object]] | None,
    *,
    resolution: int,
    exact: bool,
    default_radius: int,
    numeric_precision: int,
):
    intervals = _bounds_for_variables(variables, bounds, default_radius=default_radius)
    value_lists = [
        _grid_values_for_interval(
            lo, hi, resolution=resolution, exact=exact, numeric_precision=numeric_precision
        )
        for lo, hi in intervals
    ]
    for coords in product(*value_lists):
        yield dict(zip_equal(variables, coords, context="sampling grid point"))


def _random_points(
    variables: Sequence[sp.Symbol],
    bounds: Sequence[tuple[object, object]] | Mapping[sp.Symbol, tuple[object, object]] | None,
    *,
    attempts: int,
    seed: int | None,
    exact: bool,
    default_radius: int,
    max_denominator: int,
    numeric_precision: int,
):
    rng = _random.Random(seed)
    intervals = _bounds_for_variables(variables, bounds, default_radius=default_radius)
    for _ in range(max(0, attempts)):
        point: dict[sp.Symbol, sp.Expr] = {}
        for var, (lo, hi) in zip_equal(variables, intervals, context="sampling intervals"):
            lo_f = float(sp.N(lo))
            hi_f = float(sp.N(hi))
            u = rng.random()
            if exact:
                # Use a finite-denominator rational approximation so exact
                # validation remains possible for polynomial inequalities.
                value = sp.Rational(str(lo_f + (hi_f - lo_f) * u)).limit_denominator(
                    max_denominator
                )
            else:
                value = sp.Float(lo_f + (hi_f - lo_f) * u, numeric_precision)
            point[var] = value
        yield point


def _dedupe_points(
    points: Iterable[Mapping[sp.Symbol, sp.Expr]],
) -> tuple[dict[sp.Symbol, sp.Expr], ...]:
    out: list[dict[sp.Symbol, sp.Expr]] = []
    seen: set[tuple[tuple[sp.Symbol, sp.Expr], ...]] = set()
    for point in points:
        key = point_key(point)
        if key not in seen:
            out.append({var: sp.sympify(value) for var, value in point.items()})
            seen.add(key)
    return tuple(out)


def _validated(
    points: Iterable[Mapping[sp.Symbol, sp.Expr]], formula: sp.Expr, *, strict: bool
) -> tuple[dict[sp.Symbol, sp.Expr], ...]:
    return tuple(
        dict(point)
        for point in _dedupe_points(points)
        if satisfies_formula(formula, point, strict=strict)
    )


def _rur_points(
    formula: sp.Expr, variables: Sequence[sp.Symbol], *, count: int | None
) -> tuple[dict[sp.Symbol, sp.Expr], ...]:
    try:
        rur_result = solve_formula_with_rur(
            formula, tuple(variables), real=True, max_solutions=count
        )
    except _RECOVERABLE_ERRORS:
        return ()
    if rur_result is None or rur_result.partial:
        return ()
    return tuple(dict(point) for point in rur_result.assignments)


def _cad_points(
    formula: sp.Expr,
    variables: Sequence[sp.Symbol],
    *,
    count: int,
    strict: bool,
    per_component: bool,
) -> tuple[dict[sp.Symbol, sp.Expr], ...]:
    try:
        result = find_instance(
            formula,
            variables,
            count=count,
            domain="reals",
            strategy="cad",
            strict=strict,
            return_result=True,
        )
    except _RECOVERABLE_ERRORS:
        return ()
    return tuple(dict(point) for point in getattr(result, "instances", ()))


def sample_points(
    formula: sp.Expr,
    variables: Sequence[sp.Symbol | str] | None = None,
    *,
    count: int = 1,
    per_component: bool = False,
    strategy: str = "auto",
    strict: bool = False,
    exact: bool = True,
    seed: int | None = None,
    bounds: Sequence[tuple[object, object]]
    | Mapping[sp.Symbol, tuple[object, object]]
    | None = None,
    grid_resolution: int = 9,
    random_attempts: int | None = None,
    default_sampling_radius: int = 3,
    random_attempts_min: int = 100,
    random_attempts_each: int = 50,
    max_random_denominator: int = 1_000_000,
    numeric_precision: int = 30,
) -> tuple[dict[sp.Symbol, sp.Expr], ...]:
    """Return satisfying real sample points for a quantifier-free formula.

    ``strategy`` separates the main Sampling meanings:

    ``"representative"``/``"auto"``
        Prefer exact algebraic representatives from finite RUR systems, then
        deterministic rational witnesses, and finally CAD when explicitly asked
        via ``per_component=True`` or ``strategy="cad_cells"``.
    ``"rational"``
        Deterministic exact rational witnesses from a small grid.
    ``"grid"``
        Deterministic rational grid sampling over ``bounds``.
    ``"random"``
        Random sampling over ``bounds``. With ``exact=True`` the generated
        coordinates are rationalized; with ``exact=False`` they are SymPy Floats.
    ``"cad_cells"``
        Request CAD/instance-finder representatives.

    Every returned point is revalidated against ``formula``.
    """

    expr = to_sympy(formula) if not isinstance(formula, (sp.Basic, Boolean)) else formula
    vars_ = normalize_sampling_variables(variables, sp.sympify(expr))
    if default_sampling_radius < 0:
        raise ValueError("default_sampling_radius must be nonnegative")
    if random_attempts_min < 0 or random_attempts_each < 0:
        raise ValueError("random attempt policy values must be nonnegative")
    if max_random_denominator < 1:
        raise ValueError("max_random_denominator must be positive")
    if numeric_precision < 1:
        raise ValueError("numeric_precision must be positive")
    if count <= 0:
        return ()
    key = strategy.lower().replace("-", "_")
    aliases = {
        "auto": "representative",
        "rep": "representative",
        "representatives": "representative",
        "cad": "cad_cells",
        "cad_cell": "cad_cells",
        "cells": "cad_cells",
        "cell": "cad_cells",
        "small_grid": "rational",
        "rational_grid": "rational",
        "numeric": "random",
        "fallback": "representative",
    }
    key = aliases.get(key, key)
    if key not in {"representative", "rational", "grid", "random", "cad_cells", "complete"}:
        raise ValueError(f"unsupported sample strategy: {strategy!r}")

    candidate_points: list[dict[sp.Symbol, sp.Expr]] = []

    if key in {"representative", "complete"}:
        candidate_points.extend(_rur_points(expr, vars_, count=count))
        valid = _validated(candidate_points, expr, strict=strict)
        if len(valid) >= count:
            return valid[:count]

    if key in {"representative", "rational", "complete"}:
        for point in _small_rational_grid(vars_, radius=default_sampling_radius):
            candidate_points.append(point)
            valid = _validated(candidate_points, expr, strict=strict)
            if len(valid) >= count:
                return valid[:count]

    if key == "grid":
        for point in _bounded_rational_grid(
            vars_,
            bounds,
            resolution=grid_resolution,
            exact=exact,
            default_radius=default_sampling_radius,
            numeric_precision=numeric_precision,
        ):
            candidate_points.append(point)
            valid = _validated(candidate_points, expr, strict=strict)
            if len(valid) >= count:
                return valid[:count]

    if key == "random":
        attempts = (
            random_attempts
            if random_attempts is not None
            else max(random_attempts_min, random_attempts_each * count)
        )
        for point in _random_points(
            vars_,
            bounds,
            attempts=attempts,
            seed=seed,
            exact=exact,
            default_radius=default_sampling_radius,
            max_denominator=max_random_denominator,
            numeric_precision=numeric_precision,
        ):
            candidate_points.append(point)
            valid = _validated(candidate_points, expr, strict=strict)
            if len(valid) >= count:
                return valid[:count]

    if key in {"cad_cells", "complete"} or per_component:
        candidate_points.extend(
            _cad_points(expr, vars_, count=count, strict=strict, per_component=per_component)
        )

    valid = _validated(candidate_points, expr, strict=strict)
    return valid[:count]


def sample_point(
    formula: sp.Expr,
    variables: Sequence[sp.Symbol | str] | None = None,
    *,
    strategy: str = "auto",
    strict: bool = False,
    exact: bool = True,
    seed: int | None = None,
    bounds: Sequence[tuple[object, object]]
    | Mapping[sp.Symbol, tuple[object, object]]
    | None = None,
    grid_resolution: int = 9,
    random_attempts: int | None = None,
    default_sampling_radius: int = 3,
    random_attempts_min: int = 100,
    random_attempts_each: int = 50,
    max_random_denominator: int = 1_000_000,
    numeric_precision: int = 30,
) -> dict[sp.Symbol, sp.Expr] | None:
    """Return one satisfying real sample point for ``formula``, or ``None``."""

    points = sample_points(
        formula,
        variables,
        count=1,
        strategy=strategy,
        strict=strict,
        exact=exact,
        seed=seed,
        bounds=bounds,
        grid_resolution=grid_resolution,
        random_attempts=random_attempts,
        default_sampling_radius=default_sampling_radius,
        random_attempts_min=random_attempts_min,
        random_attempts_each=random_attempts_each,
        max_random_denominator=max_random_denominator,
        numeric_precision=numeric_precision,
    )
    return points[0] if points else None


__all__ = ["sample_point", "sample_points", "sign_at", "sign_vector"]
