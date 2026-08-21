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
from .errors import DimensionMismatchError
from .formula import to_sympy
from .instances.real_fallbacks import satisfies_formula
from .normalization import normalize_sampling_variables as _shared_variables
from .solve import find_instance


def _normalize_variables(
    variables: Sequence[sp.Symbol | str] | None,
    exprs: Iterable[sp.Expr],
) -> tuple[sp.Symbol, ...]:
    """Resolve sampling variables against symbols in the sampled formulas."""

    expr_tuple = tuple(sp.sympify(expr) for expr in exprs)
    return _shared_variables(variables, *expr_tuple)


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
) -> int:
    """Return the sign of a polynomial/expression at a point.

    The exact path handles ordinary rational SymPy values, semialg
    ``RationalSample``/``AlgebraicRoot`` objects, exact SymPy algebraic
    expressions such as ``sqrt(2)``/``RootOf``, and ``RationalUnivariatePoint``
    objects. Numeric fallback is available only when ``exact=False``.
    """

    expr = poly.as_expr() if isinstance(poly, sp.Poly) else sp.sympify(poly)
    vars_ = _normalize_variables(variables, [expr])

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
        except Exception:
            pass

    substituted = sp.cancel(
        expr.subs(dict(zip_equal(vars_, values, context="sampling substitution")))
    )
    try:
        return sign_of_algebraic_expression(substituted)
    except Exception:
        if exact:
            raise ValueError(
                f"could not determine exact sign of {sp.sstr(expr)} at {point!r}"
            ) from None

    numeric = sp.N(substituted, 120)
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
) -> tuple[int, ...] | dict[sp.Expr, int]:
    """Return the signs of ``polys`` at ``point`` in input order.

    Set ``as_dict=True`` to receive a mapping from each input expression to its
    sign. This keeps the historical tuple return by default while supporting
    more inspectable CAD/debugging workflows.
    """

    exprs = tuple(
        poly.as_expr() if isinstance(poly, sp.Poly) else sp.sympify(poly) for poly in polys
    )
    vars_ = _normalize_variables(variables, exprs)
    signs = tuple(sign_at(expr, point, variables=vars_, exact=exact) for expr in exprs)
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
    lo: sp.Expr, hi: sp.Expr, *, resolution: int, exact: bool
) -> tuple[sp.Expr, ...]:
    def finalize(value: sp.Expr) -> sp.Expr:
        return sp.simplify(value) if exact else sp.Float(value, 30)

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
):
    intervals = _bounds_for_variables(variables, bounds)
    value_lists = [
        _grid_values_for_interval(lo, hi, resolution=resolution, exact=exact)
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
):
    rng = _random.Random(seed)
    intervals = _bounds_for_variables(variables, bounds)
    for _ in range(max(0, attempts)):
        point: dict[sp.Symbol, sp.Expr] = {}
        for var, (lo, hi) in zip_equal(variables, intervals, context="sampling intervals"):
            lo_f = float(sp.N(lo))
            hi_f = float(sp.N(hi))
            u = rng.random()
            if exact:
                # Use a finite-denominator rational approximation so exact
                # validation remains possible for polynomial inequalities.
                value = sp.Rational(str(lo_f + (hi_f - lo_f) * u)).limit_denominator(10**6)
            else:
                value = sp.Float(lo_f + (hi_f - lo_f) * u, 30)
            point[var] = value
        yield point


def _dedupe_points(
    points: Iterable[Mapping[sp.Symbol, sp.Expr]],
) -> tuple[dict[sp.Symbol, sp.Expr], ...]:
    out: list[dict[sp.Symbol, sp.Expr]] = []
    seen: set[tuple[tuple[str, str], ...]] = set()
    for point in points:
        key = tuple(
            sorted((sp.sstr(var), sp.sstr(sp.simplify(value))) for var, value in point.items())
        )
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
    except Exception:
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
    except Exception:
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
    vars_ = _normalize_variables(variables, [expr])
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
        for point in _small_rational_grid(vars_):
            candidate_points.append(point)
            valid = _validated(candidate_points, expr, strict=strict)
            if len(valid) >= count:
                return valid[:count]

    if key == "grid":
        for point in _bounded_rational_grid(vars_, bounds, resolution=grid_resolution, exact=exact):
            candidate_points.append(point)
            valid = _validated(candidate_points, expr, strict=strict)
            if len(valid) >= count:
                return valid[:count]

    if key == "random":
        attempts = random_attempts if random_attempts is not None else max(100, 50 * count)
        for point in _random_points(vars_, bounds, attempts=attempts, seed=seed, exact=exact):
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
    )
    return points[0] if points else None


__all__ = ["sample_point", "sample_points", "sign_at", "sign_vector"]
