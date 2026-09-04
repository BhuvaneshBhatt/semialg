from __future__ import annotations

from collections.abc import Mapping, Sequence

import sympy as sp
from sympy.core.relational import Equality, Unequality

from ._errors import EXACT_OPERATION_ERRORS as _RECOVERABLE_ERRORS
from ._linear_relations import certified_sign, safe_linear_solution
from .interval_decomposition import finite_real_roots as _finite_real_roots  # noqa: F401
from .interval_decomposition import (
    truth_at as _truth_at,
)
from .normalization import (
    require_conjunction,
)
from .region_integral_results import (
    RegionIntegralPiece,
)
from .region_integrate import _one_dimensional_intervals
from .relations import make_zero_relation
from .relations import split_relation as _relation_parts


def _radial_radii_squared(
    condition: sp.Expr, x: sp.Symbol, y: sp.Symbol
) -> tuple[sp.Expr, sp.Expr] | None:
    """Extract exact inner and outer squared radii for a radial semialgebraic region."""
    try:
        atoms = require_conjunction(
            condition, message="this reconstruction path supports conjunctions only"
        )
    except NotImplementedError:
        return None
    lower = sp.Integer(0)
    upper = sp.oo
    found = False
    for atom in atoms:
        if atom is sp.false:
            return (sp.Integer(0), sp.Integer(0))
        if not getattr(atom, "is_Relational", False):
            return None
        expr, op = _relation_parts(atom)
        try:
            poly = sp.Poly(expr, x, y)
        except _RECOVERABLE_ERRORS:
            return None
        coeff_x2 = poly.coeff_monomial(x**2)
        coeff_y2 = poly.coeff_monomial(y**2)
        if coeff_x2 == 0 or sp.simplify(coeff_x2 - coeff_y2) != 0:
            return None
        if any(monom not in {(2, 0), (0, 2), (0, 0)} for monom in poly.monoms()):
            return None
        constant = poly.coeff_monomial(1)
        radius_sq = sp.simplify(-constant / coeff_x2)
        if op == "==":
            return (sp.Integer(0), sp.Integer(0))
        if op == "!=":
            found = True
            continue
        if sp.simplify(coeff_x2).is_negative:
            op = {"<": ">", "<=": ">=", ">": "<", ">=": "<="}[op]
        if op in ("<", "<="):
            upper = sp.Min(upper, radius_sq) if upper != sp.oo else radius_sq
            found = True
        elif op in (">", ">="):
            lower = sp.Max(lower, radius_sq)
            found = True
    if not found:
        return None
    lower = sp.simplify(lower)
    upper = sp.simplify(upper)
    if upper == sp.oo:
        raise NotImplementedError("unbounded radial integrals are outside the supported fragment")
    if bool(upper <= lower) or bool(upper <= 0):
        return (sp.Integer(0), sp.Integer(0))
    if bool(lower < 0):
        lower = sp.Integer(0)
    return lower, upper


def _angular_monomial_integral(i: int, j: int) -> sp.Expr:
    if i % 2 or j % 2:
        return sp.Integer(0)
    return sp.simplify(
        2
        * sp.gamma(sp.Rational(i + 1, 2))
        * sp.gamma(sp.Rational(j + 1, 2))
        / sp.gamma(sp.Rational(i + j + 2, 2))
    )


def _integrate_radial_polynomial(
    integrand: sp.Expr, condition: sp.Expr, x: sp.Symbol, y: sp.Symbol
) -> sp.Expr | None:
    radii = _radial_radii_squared(condition, x, y)
    if radii is None:
        return None
    lower_sq, upper_sq = radii
    if lower_sq == upper_sq:
        return sp.Integer(0)
    poly = sp.Poly(sp.expand(integrand), x, y)
    if set(poly.as_expr().free_symbols) - {x, y}:
        return None
    total = sp.Integer(0)
    for (i, j), coeff in poly.terms():
        angular = _angular_monomial_integral(i, j)
        if angular == 0:
            continue
        power = sp.Rational(i + j + 2, 2)
        radial = sp.simplify((upper_sq**power - lower_sq**power) / (i + j + 2))
        total += coeff * angular * radial
    return sp.simplify(total)


def _vertical_slice_data(
    condition: sp.Expr,
    x: sp.Symbol,
    y: sp.Symbol,
    bounds: Mapping[sp.Symbol, tuple[sp.Expr, sp.Expr]],
) -> tuple[sp.Expr, sp.Expr, tuple[tuple[sp.Expr, sp.Expr], ...]] | None:
    """Extract exact vertical slice bounds and residual base conditions for integration."""
    try:
        atoms = require_conjunction(
            condition, message="this reconstruction path supports conjunctions only"
        )
    except NotImplementedError:
        return None
    lower_bounds: list[sp.Expr] = []
    upper_bounds: list[sp.Expr] = []
    x_conditions: list[sp.Expr] = []
    for atom in atoms:
        if atom is sp.false:
            return (sp.Integer(0), sp.Integer(0), ())
        expr, op = _relation_parts(atom)
        if op == "==":
            if sp.simplify(expr) != 0:
                return (sp.Integer(0), sp.Integer(0), ())
            continue
        if op == "!=":
            continue
        if y not in expr.free_symbols:
            x_conditions.append(make_zero_relation(expr, op))
            continue
        try:
            poly_y = sp.Poly(expr, y)
        except _RECOVERABLE_ERRORS:
            return None
        if poly_y.degree() != 1:
            return None
        coeff = sp.simplify(poly_y.coeff_monomial(y))
        sign = certified_sign(coeff)
        if sign not in (-1, 1):
            return None
        boundary = safe_linear_solution(expr, y)
        if boundary is None:
            return None
        if sign < 0:
            flipped = {"<": ">", "<=": ">=", ">": "<", ">=": "<="}[op]
        else:
            flipped = op
        if flipped in ("<", "<="):
            upper_bounds.append(boundary)
        elif flipped in (">", ">="):
            lower_bounds.append(boundary)
    if len(lower_bounds) != 1 or len(upper_bounds) != 1:
        return None
    lower = lower_bounds[0]
    upper = upper_bounds[0]
    x_conditions.append(lower <= upper)
    if x in bounds:
        lo, hi = bounds[x]
        x_conditions.extend([x >= lo, x <= hi])
    if y in bounds:
        lo, hi = bounds[y]
        # Keep slice heights explicit rather than introducing piecewise max/min expressions.
        if sp.simplify(lower - lo) != 0 or sp.simplify(upper - hi) != 0:
            x_conditions.extend([lower >= lo, upper <= hi])
    x_condition = sp.And(*x_conditions) if x_conditions else sp.true
    intervals = _one_dimensional_intervals(x_condition, x, None)
    return lower, upper, intervals


def _merge_interval_bound(
    lower: sp.Expr,
    upper: sp.Expr,
    relation_expr: sp.Expr,
    op: str,
    variable: sp.Symbol,
    *,
    parameter_symbols: Sequence[sp.Symbol] = (),
) -> tuple[sp.Expr, sp.Expr] | None:
    """Update a one-variable interval bound from a linear relation."""

    try:
        poly = sp.Poly(relation_expr, variable)
    except _RECOVERABLE_ERRORS:
        return None
    allowed_symbols = {variable, *parameter_symbols}
    if poly.degree() > 1 or relation_expr.free_symbols - allowed_symbols:
        return None
    coeff = sp.simplify(poly.coeff_monomial(variable))
    sign = certified_sign(coeff)
    if sign == 0:
        return (lower, upper) if _truth_at(make_zero_relation(relation_expr, op), {}) else None
    if sign not in (-1, 1):
        return None
    boundary = safe_linear_solution(relation_expr, variable)
    if boundary is None:
        return None
    normalized_op = op
    if sign < 0:
        normalized_op = {"<": ">", "<=": ">=", ">": "<", ">=": "<="}[op]
    if normalized_op in ("<", "<="):
        upper = boundary if upper == sp.oo else sp.Min(upper, boundary)
    elif normalized_op in (">", ">="):
        lower = boundary if lower == -sp.oo else sp.Max(lower, boundary)
    else:
        return None
    return sp.simplify(lower), sp.simplify(upper)


def _box_limits_from_condition(
    condition: sp.Expr,
    variables: Sequence[sp.Symbol],
    bounds: Mapping[sp.Symbol, tuple[sp.Expr, sp.Expr]],
    *,
    parameter_symbols: Sequence[sp.Symbol] = (),
) -> tuple[tuple[sp.Symbol, sp.Expr, sp.Expr], ...] | None:
    """Recognize axis-aligned boxes from independent linear bounds."""

    try:
        atoms = require_conjunction(
            condition, message="this reconstruction path supports conjunctions only"
        )
    except NotImplementedError:
        return None
    lower: dict[sp.Symbol, sp.Expr] = {
        var: bounds.get(var, (-sp.oo, sp.oo))[0] for var in variables
    }
    upper: dict[sp.Symbol, sp.Expr] = {
        var: bounds.get(var, (-sp.oo, sp.oo))[1] for var in variables
    }
    saw_condition = bool(bounds)
    for atom in atoms:
        if atom is sp.false:
            return ()
        if isinstance(atom, (Equality, Unequality)):
            return None
        if not getattr(atom, "is_Relational", False):
            return None
        expr, op = _relation_parts(atom)
        involved = [var for var in variables if var in expr.free_symbols]
        if len(involved) != 1:
            return None
        if expr.free_symbols - set(variables) - set(parameter_symbols):
            return None
        var = involved[0]
        merged = _merge_interval_bound(
            lower[var], upper[var], expr, op, var, parameter_symbols=parameter_symbols
        )
        if merged is None:
            return None
        lower[var], upper[var] = merged
        saw_condition = True
    if not saw_condition:
        return None
    limits: list[tuple[sp.Symbol, sp.Expr, sp.Expr]] = []
    for var in reversed(tuple(variables)):
        lo = sp.simplify(lower[var])
        hi = sp.simplify(upper[var])
        if lo == -sp.oo or hi == sp.oo:
            return None
        width = sp.simplify(hi - lo)
        if width.is_negative is True:
            return ()
        limits.append((var, lo, hi))
    return tuple(limits)


def _reduce_box(
    integrand: sp.Expr,
    condition: sp.Expr,
    variables: Sequence[sp.Symbol],
    bounds: Mapping[sp.Symbol, tuple[sp.Expr, sp.Expr]],
    *,
    parameter_symbols: Sequence[sp.Symbol] = (),
) -> tuple[RegionIntegralPiece, ...] | None:
    limits = _box_limits_from_condition(
        condition, variables, bounds, parameter_symbols=parameter_symbols
    )
    if limits is None:
        return None
    if limits == ():
        return ()
    return (
        RegionIntegralPiece(
            integrand=integrand,
            limits=limits,
            method="axis_aligned_box_iterated_integral",
            diagnostics={"shape": "box"},
        ),
    )


def _canonical_linear_coefficients(
    expr: sp.Expr, variables: Sequence[sp.Symbol]
) -> tuple[tuple[sp.Expr, ...], sp.Expr] | None:
    try:
        poly = sp.Poly(expr, *variables)
    except _RECOVERABLE_ERRORS:
        return None
    if poly.total_degree() > 1:
        return None
    coeffs = tuple(sp.simplify(poly.coeff_monomial(var)) for var in variables)
    constant = sp.simplify(poly.coeff_monomial(1))
    return coeffs, constant


def _reduce_standard_simplex_2d(
    integrand: sp.Expr,
    condition: sp.Expr,
    x: sp.Symbol,
    y: sp.Symbol,
) -> tuple[RegionIntegralPiece, ...] | None:
    """Recognize the 2D unit simplex x >= 0, y >= 0, x + y <= 1."""

    try:
        atoms = require_conjunction(
            condition, message="this reconstruction path supports conjunctions only"
        )
    except NotImplementedError:
        return None
    needed = {"x_nonnegative": False, "y_nonnegative": False, "sum_le_one": False}
    for atom in atoms:
        if not getattr(atom, "is_Relational", False) or isinstance(atom, (Equality, Unequality)):
            return None
        expr, op = _relation_parts(atom)
        linear = _canonical_linear_coefficients(expr, (x, y))
        if linear is None:
            return None
        coeffs, constant = linear
        signature = (*coeffs, constant)
        if (signature == (sp.Integer(1), sp.Integer(0), sp.Integer(0)) and op in (">", ">=")) or (
            signature == (sp.Integer(-1), sp.Integer(0), sp.Integer(0)) and op in ("<", "<=")
        ):
            needed["x_nonnegative"] = True
        elif (signature == (sp.Integer(0), sp.Integer(1), sp.Integer(0)) and op in (">", ">=")) or (
            signature == (sp.Integer(0), sp.Integer(-1), sp.Integer(0)) and op in ("<", "<=")
        ):
            needed["y_nonnegative"] = True
        elif (
            signature == (sp.Integer(1), sp.Integer(1), sp.Integer(-1)) and op in ("<", "<=")
        ) or (signature == (sp.Integer(-1), sp.Integer(-1), sp.Integer(1)) and op in (">", ">=")):
            needed["sum_le_one"] = True
        else:
            return None
    if not all(needed.values()):
        return None
    return (
        RegionIntegralPiece(
            integrand=integrand,
            limits=((y, sp.Integer(0), 1 - x), (x, sp.Integer(0), sp.Integer(1))),
            method="unit_simplex_iterated_integral",
            diagnostics={"shape": "unit_simplex_2d"},
        ),
    )


def _axis_aligned_ellipse_data(
    condition: sp.Expr, x: sp.Symbol, y: sp.Symbol
) -> tuple[sp.Expr, sp.Expr, sp.Expr, sp.Expr] | None:
    """Recognize a filled axis-aligned ellipse from one quadratic inequality."""

    try:
        atoms = require_conjunction(
            condition, message="this reconstruction path supports conjunctions only"
        )
    except NotImplementedError:
        return None
    if len(atoms) != 1 or not getattr(atoms[0], "is_Relational", False):
        return None
    atom = atoms[0]
    if isinstance(atom, (Equality, Unequality)):
        return None
    expr, op = _relation_parts(atom)
    if op not in ("<", "<=", ">", ">="):
        return None
    try:
        poly = sp.Poly(expr, x, y)
    except _RECOVERABLE_ERRORS:
        return None
    if poly.total_degree() != 2 or poly.coeff_monomial(x * y) != 0:
        return None
    ax = sp.simplify(poly.coeff_monomial(x**2))
    ay = sp.simplify(poly.coeff_monomial(y**2))
    bx = sp.simplify(poly.coeff_monomial(x))
    by = sp.simplify(poly.coeff_monomial(y))
    c0 = sp.simplify(poly.coeff_monomial(1))
    if ax == 0 or ay == 0:
        return None
    if op in (">", ">="):
        ax, ay, bx, by, c0 = (-ax, -ay, -bx, -by, -c0)
    if not (ax.is_positive and ay.is_positive):
        return None
    cx = sp.simplify(-bx / (2 * ax))
    cy = sp.simplify(-by / (2 * ay))
    completed_constant = sp.simplify(c0 - ax * cx**2 - ay * cy**2)
    if not completed_constant.is_negative:
        return None
    level = sp.simplify(-completed_constant)
    rx2 = sp.simplify(level / ax)
    ry2 = sp.simplify(level / ay)
    if bool(rx2 <= 0) or bool(ry2 <= 0):
        return None
    return cx, cy, rx2, ry2


def _reduce_axis_aligned_ellipse(
    integrand: sp.Expr,
    condition: sp.Expr,
    x: sp.Symbol,
    y: sp.Symbol,
) -> tuple[RegionIntegralPiece, ...] | None:
    data = _axis_aligned_ellipse_data(condition, x, y)
    if data is None:
        return None
    cx, cy, rx2, ry2 = data
    u = sp.Symbol(f"_{x.name}_unit", real=True)
    v = sp.Symbol(f"_{y.name}_unit", real=True)
    rx = sp.sqrt(rx2)
    ry = sp.sqrt(ry2)
    transformed = sp.simplify(rx * ry * integrand.subs({x: cx + rx * u, y: cy + ry * v}))
    height = sp.sqrt(1 - u**2)
    return (
        RegionIntegralPiece(
            integrand=transformed,
            limits=((v, -height, height), (u, -1, 1)),
            method="axis_aligned_ellipse_affine_unit_disk",
            diagnostics={
                "shape": "axis_aligned_ellipse",
                "center": (cx, cy),
                "radii_squared": (rx2, ry2),
            },
        ),
    )


def _integrate_axis_aligned_ellipse_polynomial(
    integrand: sp.Expr,
    condition: sp.Expr,
    x: sp.Symbol,
    y: sp.Symbol,
) -> sp.Expr | None:
    """Fast exact polynomial moments over a recognized axis-aligned ellipse."""

    data = _axis_aligned_ellipse_data(condition, x, y)
    if data is None:
        return None
    cx, cy, rx2, ry2 = data
    u = sp.Symbol(f"_{x.name}_unit", real=True)
    v = sp.Symbol(f"_{y.name}_unit", real=True)
    rx = sp.sqrt(rx2)
    ry = sp.sqrt(ry2)
    transformed = sp.expand(rx * ry * integrand.subs({x: cx + rx * u, y: cy + ry * v}))
    return _integrate_radial_polynomial(transformed, u**2 + v**2 <= 1, u, v)
