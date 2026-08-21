from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass

import sympy as sp
from sympy.core.relational import Equality, Unequality
from sympy.logic.boolalg import And as SymAnd
from sympy.logic.boolalg import BooleanFalse, BooleanTrue
from sympy.logic.boolalg import Not as SymNot
from sympy.logic.boolalg import Or as SymOr

from .implicit_geometry import decompose_cylindrical_formula_to_vertical_bounds_2d
from .interval_decomposition import (
    finite_real_roots as _finite_real_roots,
)
from .interval_decomposition import (
    one_dimensional_intervals as _shared_intervals,
)
from .interval_decomposition import (
    truth_at as _truth_at,
)
from .normalization import (
    normalize_bounds as _normalize_bounds,
)
from .normalization import (
    normalize_formula as _normalize_formula,
)
from .normalization import (
    normalize_variables as _normalize_variables,
)
from .relations import make_zero_relation
from .relations import split_relation as _relation_parts
from .standard_region_integrate import integrate_over_standard_region
from .standard_regions import StandardRegion

_RECOVERABLE_ERRORS = (
    ArithmeticError,
    TypeError,
    ValueError,
    NotImplementedError,
    RuntimeError,
    sp.PolynomialError,
)


@dataclass(frozen=True)
class RegionIntegralResult:
    """Integral over a supported semialgebraic region.

    Integrals use ambient Lebesgue measure by default. Lower-dimensional
    equality-only subsets therefore contribute zero unless the caller requests
    ``measure_dimension="intrinsic"`` or an explicit Hausdorff dimension.

    ``exact`` is true for symbolic evaluation and false for numerical
    evaluation. ``evaluated`` records whether the returned value has been
    evaluated rather than left as an explicit ``Integral`` expression.
    """

    value: sp.Expr
    integrand: sp.Expr
    condition: sp.Expr
    variables: tuple[sp.Symbol, ...]
    method: str
    diagnostics: Mapping[str, object] | None = None
    exact: bool = True
    evaluated: bool = True
    error_estimate: sp.Expr | None = None


@dataclass(frozen=True)
class RegionIntegralPiece:
    """One iterated-integral piece for a semialgebraic region integral.

    ``limits`` uses SymPy's integral-limit convention, for example
    ``((x, a, b), (y, lower(x), upper(x)))``. Pieces may be signed: Boolean
    differences such as annuli are represented by using a negative integrand on
    the subtracted part.
    """

    integrand: sp.Expr
    limits: tuple[tuple[sp.Symbol, sp.Expr, sp.Expr], ...]
    method: str
    diagnostics: Mapping[str, object] | None = None

    def as_integral(self) -> sp.Integral:
        """Return this piece as an unevaluated SymPy ``Integral``."""

        return sp.Integral(self.integrand, *self.limits)


@dataclass(frozen=True)
class ReducedRegionIntegral:
    """Reduction of a region integral to explicit iterated-integral pieces."""

    integrand: sp.Expr
    condition: sp.Expr
    variables: tuple[sp.Symbol, ...]
    pieces: tuple[RegionIntegralPiece, ...]
    method: str
    diagnostics: Mapping[str, object] | None = None

    def as_integrals(self) -> tuple[sp.Integral, ...]:
        """Return each reduced piece as an unevaluated SymPy ``Integral``."""

        return tuple(piece.as_integral() for piece in self.pieces)

    def unevaluated_sum(self) -> sp.Expr:
        """Return the formal sum of unevaluated piece integrals."""

        if not self.pieces:
            return sp.Integer(0)
        return sp.Add(*self.as_integrals())


def _atoms(condition: sp.Expr) -> tuple[sp.Expr, ...]:
    if condition is sp.true or isinstance(condition, BooleanTrue):
        return ()
    if condition is sp.false or isinstance(condition, BooleanFalse):
        return (sp.false,)
    if isinstance(condition, SymAnd):
        out: list[sp.Expr] = []
        for arg in condition.args:
            out.extend(_atoms(arg))
        return tuple(out)
    if isinstance(condition, (SymOr, SymNot)):
        raise NotImplementedError("this reconstruction path currently supports conjunctions only")
    if getattr(condition, "is_Relational", False):
        return (condition,)
    raise TypeError(f"unsupported formula expression: {condition!r}")


def _one_dimensional_intervals(
    condition: sp.Expr,
    variable: sp.Symbol,
    bound: tuple[sp.Expr, sp.Expr] | None,
) -> tuple[tuple[sp.Expr, sp.Expr], ...]:
    """Return exact true intervals for a one-dimensional integral."""

    return _shared_intervals(
        condition,
        variable,
        bound,
        extra_symbol_error=(
            "1D integration formula contains symbols outside the integration variable"
        ),
    )


def _integrate_1d(
    integrand: sp.Expr,
    condition: sp.Expr,
    variable: sp.Symbol,
    bound: tuple[sp.Expr, sp.Expr] | None,
) -> sp.Expr:
    total = sp.Integer(0)
    for lo, hi in _one_dimensional_intervals(condition, variable, bound):
        if lo == -sp.oo or hi == sp.oo:
            val = sp.integrate(integrand, (variable, lo, hi))
        else:
            val = sp.integrate(integrand, (variable, lo, hi))
        if isinstance(val, sp.Integral):
            raise NotImplementedError("SymPy could not evaluate one of the exact 1D integrals")
        total += val
    return sp.simplify(total)


def _radial_radii_squared(
    condition: sp.Expr, x: sp.Symbol, y: sp.Symbol
) -> tuple[sp.Expr, sp.Expr] | None:
    try:
        atoms = _atoms(condition)
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
            upper = sp.Min(upper, radius_sq) if upper is not sp.oo else radius_sq
            found = True
        elif op in (">", ">="):
            lower = sp.Max(lower, radius_sq)
            found = True
    if not found:
        return None
    lower = sp.simplify(lower)
    upper = sp.simplify(upper)
    if upper == sp.oo:
        raise NotImplementedError("unbounded radial integrals are not yet supported")
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
        atoms = _atoms(condition)
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
        rest = sp.simplify(poly_y.as_expr() - coeff * y)
        boundary = sp.simplify(-rest / coeff)
        if coeff.is_negative:
            flipped = {"<": ">", "<=": ">=", ">": "<", ">=": "<="}[op]
        elif coeff.is_positive:
            flipped = op
        else:
            return None
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
        # The first implementation avoids piecewise max/min in the slice height.
        if sp.simplify(lower - lo) != 0 or sp.simplify(upper - hi) != 0:
            x_conditions.extend([lower >= lo, upper <= hi])
    x_condition = sp.And(*x_conditions) if x_conditions else sp.true
    intervals = _one_dimensional_intervals(x_condition, x, None)
    return lower, upper, intervals


def _integrate_vertical_slice(
    integrand: sp.Expr,
    condition: sp.Expr,
    x: sp.Symbol,
    y: sp.Symbol,
    bounds: Mapping[sp.Symbol, tuple[sp.Expr, sp.Expr]],
) -> sp.Expr | None:
    data = _vertical_slice_data(condition, x, y, bounds)
    if data is None:
        return None
    lower, upper, intervals = data
    if not intervals:
        return sp.Integer(0)
    inner = sp.integrate(integrand, (y, lower, upper))
    if isinstance(inner, sp.Integral):
        raise NotImplementedError("SymPy could not evaluate the inner vertical-slice integral")
    total = sp.Integer(0)
    for lo, hi in intervals:
        if lo == -sp.oo or hi == sp.oo:
            raise NotImplementedError("unbounded vertical-slice integrals are not yet supported")
        val = sp.integrate(inner, (x, lo, hi))
        if isinstance(val, sp.Integral):
            raise NotImplementedError("SymPy could not evaluate the outer vertical-slice integral")
        total += val
    return sp.simplify(total)


def _merge_interval_bound(
    lower: sp.Expr,
    upper: sp.Expr,
    relation_expr: sp.Expr,
    op: str,
    variable: sp.Symbol,
) -> tuple[sp.Expr, sp.Expr] | None:
    """Update a one-variable interval bound from a linear relation."""

    try:
        poly = sp.Poly(relation_expr, variable)
    except _RECOVERABLE_ERRORS:
        return None
    if poly.degree() > 1 or relation_expr.free_symbols - {variable}:
        return None
    coeff = sp.simplify(poly.coeff_monomial(variable))
    if coeff == 0:
        return (lower, upper) if _truth_at(make_zero_relation(relation_expr, op), {}) else None
    rest = sp.simplify(poly.as_expr() - coeff * variable)
    boundary = sp.simplify(-rest / coeff)
    normalized_op = op
    if coeff.is_negative:
        normalized_op = {"<": ">", "<=": ">=", ">": "<", ">=": "<="}[op]
    elif not coeff.is_positive:
        return None
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
) -> tuple[tuple[sp.Symbol, sp.Expr, sp.Expr], ...] | None:
    """Recognize axis-aligned boxes from independent linear bounds."""

    try:
        atoms = _atoms(condition)
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
        var = involved[0]
        merged = _merge_interval_bound(lower[var], upper[var], expr, op, var)
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
        if bool(sp.simplify(hi - lo) < 0):
            return ()
        limits.append((var, lo, hi))
    return tuple(limits)


def _reduce_box(
    integrand: sp.Expr,
    condition: sp.Expr,
    variables: Sequence[sp.Symbol],
    bounds: Mapping[sp.Symbol, tuple[sp.Expr, sp.Expr]],
) -> tuple[RegionIntegralPiece, ...] | None:
    limits = _box_limits_from_condition(condition, variables, bounds)
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
        atoms = _atoms(condition)
    except NotImplementedError:
        return None
    needed = {"x_nonnegative": False, "y_nonnegative": False, "sum_le_one": False}
    for atom in atoms:
        if not getattr(atom, "is_Relational", False) or isinstance(atom, (Equality, Unequality)):
            return None
        expr, op = _relation_parts(atom)
        if _canonical_linear_coefficients(expr, (x, y)) is None:
            return None
        if sp.simplify(expr - x) == 0 and op in (">", ">="):
            needed["x_nonnegative"] = True
        elif sp.simplify(expr + x) == 0 and op in ("<", "<="):
            needed["x_nonnegative"] = True
        elif sp.simplify(expr - y) == 0 and op in (">", ">="):
            needed["y_nonnegative"] = True
        elif sp.simplify(expr + y) == 0 and op in ("<", "<="):
            needed["y_nonnegative"] = True
        elif sp.simplify(expr - (x + y - 1)) == 0 and op in ("<", "<="):
            needed["sum_le_one"] = True
        elif sp.simplify(expr + (x + y - 1)) == 0 and op in (">", ">="):
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
        atoms = _atoms(condition)
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


def _integral_piece_value(piece: RegionIntegralPiece) -> sp.Expr:
    value = sp.integrate(piece.integrand, *piece.limits)
    if isinstance(value, sp.Integral) or value.has(sp.Integral):
        raise NotImplementedError("SymPy could not evaluate one of the reduced iterated integrals")
    return sp.simplify(value)


def _numeric_piece_value(piece: RegionIntegralPiece, precision: int) -> sp.Expr:
    value = piece.as_integral().evalf(precision)
    if isinstance(value, sp.Integral) or value.has(sp.Integral):
        raise NotImplementedError(
            "SymPy could not numerically evaluate one of the reduced iterated integrals"
        )
    return value


def _evaluate_reduced_integral(
    reduced: ReducedRegionIntegral,
    *,
    method: str,
    precision: int,
) -> tuple[sp.Expr, bool, str]:
    """Evaluate reduced pieces according to the requested evaluation method.

    Returns ``(value, exact, evaluation_method)``. ``method="symbolic"``
    never falls back to numerical quadrature. ``method="auto"`` tries the
    symbolic path first and falls back to numerical evaluation only if symbolic
    integration leaves an unevaluated Integral.
    """

    if method not in {"symbolic", "numeric", "auto"}:
        raise ValueError('method must be one of "symbolic", "numeric", or "auto"')

    if method in {"symbolic", "auto"}:
        try:
            symbolic = sp.simplify(
                sum((_integral_piece_value(piece) for piece in reduced.pieces), sp.Integer(0))
            )
            return symbolic, True, "symbolic"
        except NotImplementedError:
            if method == "symbolic":
                raise

    numeric = sum(
        (_numeric_piece_value(piece, precision) for piece in reduced.pieces), sp.Float(0, precision)
    )
    return sp.N(numeric, precision), False, "numeric"


def _reduce_1d(
    integrand: sp.Expr,
    condition: sp.Expr,
    variable: sp.Symbol,
    bound: tuple[sp.Expr, sp.Expr] | None,
) -> tuple[RegionIntegralPiece, ...]:
    pieces: list[RegionIntegralPiece] = []
    for lo, hi in _one_dimensional_intervals(condition, variable, bound):
        pieces.append(
            RegionIntegralPiece(
                integrand=integrand,
                limits=((variable, lo, hi),),
                method="one_dimensional_cell_interval",
            )
        )
    return tuple(pieces)


def _reduce_radial_vertical_pieces(
    integrand: sp.Expr,
    condition: sp.Expr,
    x: sp.Symbol,
    y: sp.Symbol,
) -> tuple[RegionIntegralPiece, ...] | None:
    radii = _radial_radii_squared(condition, x, y)
    if radii is None:
        return None
    lower_sq, upper_sq = radii
    if lower_sq == upper_sq:
        return ()
    if upper_sq == sp.oo:
        raise NotImplementedError(
            "unbounded radial regions cannot yet be reduced to finite iterated pieces"
        )

    pieces: list[RegionIntegralPiece] = []

    def add_disk_piece(radius_sq: sp.Expr, signed_integrand: sp.Expr, label: str) -> None:
        if sp.simplify(radius_sq) == 0:
            return
        radius = sp.sqrt(radius_sq)
        height = sp.sqrt(radius_sq - x**2)
        pieces.append(
            RegionIntegralPiece(
                integrand=signed_integrand,
                limits=((y, -height, height), (x, -radius, radius)),
                method=label,
                diagnostics={"radius_squared": radius_sq},
            )
        )

    add_disk_piece(upper_sq, integrand, "radial_region_outer_disk_vertical_slice")
    if sp.simplify(lower_sq) != 0:
        add_disk_piece(lower_sq, -integrand, "radial_region_inner_disk_subtraction")
    return tuple(pieces)


def _reduce_vertical_slice(
    integrand: sp.Expr,
    condition: sp.Expr,
    x: sp.Symbol,
    y: sp.Symbol,
    bounds: Mapping[sp.Symbol, tuple[sp.Expr, sp.Expr]],
) -> tuple[RegionIntegralPiece, ...] | None:
    data = _vertical_slice_data(condition, x, y, bounds)
    if data is None:
        return None
    lower, upper, intervals = data
    return tuple(
        RegionIntegralPiece(
            integrand=integrand,
            limits=((y, lower, upper), (x, lo, hi)),
            method="vertical_slice_iterated_integral",
            diagnostics={"lower": lower, "upper": upper},
        )
        for lo, hi in intervals
    )


def _reduce_cylindrical_vertical_bounds_2d(
    integrand: sp.Expr,
    condition: sp.Expr,
    x: sp.Symbol,
    y: sp.Symbol,
) -> tuple[RegionIntegralPiece, ...] | None:
    """Reduce supported CAD-like 2D cylindrical formulas to vertical pieces."""

    try:
        cells = decompose_cylindrical_formula_to_vertical_bounds_2d(condition, (x, y))
    except NotImplementedError:
        return None
    pieces: list[RegionIntegralPiece] = []
    for cell in cells:
        lo, hi = cell.x_interval
        if lo == hi:
            continue
        if lo == -sp.oo or hi == sp.oo:
            return None
        for lower, upper in cell.y_bounds:
            if lower == upper:
                continue
            if lower == -sp.oo or upper == sp.oo:
                return None
            pieces.append(
                RegionIntegralPiece(
                    integrand=integrand,
                    limits=((y, lower, upper), (x, lo, hi)),
                    method="cylindrical_formula_vertical_bounds_2d",
                    diagnostics={
                        "x_interval": (lo, hi),
                        "y_bounds": (lower, upper),
                        "source_formula": cell.source_formula,
                    },
                )
            )
    return tuple(pieces)


def _reduce_cad_extracted_vertical_bounds_2d(
    integrand: sp.Expr,
    condition: sp.Expr,
    x: sp.Symbol,
    y: sp.Symbol,
) -> tuple[RegionIntegralPiece, ...] | None:
    """Reduce arbitrary complete-CAD 2D cells to vertical integral pieces."""

    try:
        from .cad.cells import extract_vertical_bounds_from_cad_2d

        cells = extract_vertical_bounds_from_cad_2d(condition, (x, y), full_dimensional_only=True)
    except _RECOVERABLE_ERRORS:
        return None
    pieces: list[RegionIntegralPiece] = []
    for cell in cells:
        lo, hi = cell.x_interval
        if lo == hi or lo == -sp.oo or hi == sp.oo:
            continue
        for lower, upper in cell.y_bounds:
            if lower == upper or lower == -sp.oo or upper == sp.oo:
                continue
            pieces.append(
                RegionIntegralPiece(
                    integrand=integrand,
                    limits=((y, lower, upper), (x, lo, hi)),
                    method="complete_cad_vertical_bounds_2d",
                    diagnostics={
                        "x_interval": (lo, hi),
                        "y_bounds": (lower, upper),
                        "source_formula": cell.source_formula,
                    },
                )
            )
    return tuple(pieces) if pieces else None


def _normalize_measure_dimension(
    measure_dimension: object, ambient_dimension: int, condition: sp.Expr
) -> int | str:
    """Normalize measure-dimension options for the first intrinsic layer."""

    if measure_dimension in (None, "ambient"):
        return ambient_dimension
    if measure_dimension == "top":
        return _infer_region_dimension(condition, ambient_dimension)
    if measure_dimension == "intrinsic":
        return _infer_region_dimension(condition, ambient_dimension)
    if isinstance(measure_dimension, int):
        if measure_dimension < 0 or measure_dimension > ambient_dimension:
            raise ValueError("measure_dimension must be between 0 and the ambient dimension")
        return measure_dimension
    raise ValueError('measure_dimension must be None, "ambient", "intrinsic", "top", or an integer')


def _infer_region_dimension(condition: sp.Expr, ambient_dimension: int) -> int:
    """Infer a useful dimension for common semialgebraic regions.

    This intentionally small heuristic recognizes full-dimensional inequality
    regions, univariate finite point sets, two-dimensional finite point sets,
    and plane algebraic curves. It is not a replacement for CAD-based cell
    dimension computation; unsupported cases default to the ambient dimension.
    """

    try:
        atoms = _atoms(condition)
    except _RECOVERABLE_ERRORS:
        return ambient_dimension
    nontrivial_equalities: list[sp.Expr] = []
    for atom in atoms:
        if isinstance(atom, Equality):
            expr, _ = _relation_parts(atom)
            if sp.simplify(expr) != 0:
                nontrivial_equalities.append(expr)
    if not nontrivial_equalities:
        return ambient_dimension
    if ambient_dimension == 1:
        return 0
    if ambient_dimension == 2:
        if len(nontrivial_equalities) >= 2:
            return 0
        return 1
    return max(0, ambient_dimension - len(nontrivial_equalities))


def _zero_dimensional_points(
    condition: sp.Expr, variables: Sequence[sp.Symbol]
) -> tuple[dict[sp.Symbol, sp.Expr], ...]:
    """Return finite real points for common zero-dimensional formulas."""

    try:
        atoms = _atoms(condition)
    except _RECOVERABLE_ERRORS as exc:
        raise NotImplementedError(
            "zero-dimensional integration currently supports conjunctions of relations"
        ) from exc

    equalities: list[sp.Expr] = []
    for atom in atoms:
        if isinstance(atom, Equality):
            expr, _ = _relation_parts(atom)
            if sp.simplify(expr) != 0:
                equalities.append(expr)
    if not equalities:
        raise NotImplementedError(
            "zero-dimensional integration requires at least one nontrivial equality"
        )

    if len(variables) == 1:
        var = variables[0]
        roots: set[str] = set()
        points: list[dict[sp.Symbol, sp.Expr]] = []
        for eq in equalities:
            if not eq.free_symbols <= {var}:
                raise NotImplementedError(
                    "zero-dimensional univariate formulas must only use the integration variable"
                )
            for root in _finite_real_roots(eq, var):
                subs = {var: root}
                if _truth_at(condition, subs):
                    key = sp.sstr(sp.simplify(root))
                    if key not in roots:
                        roots.add(key)
                        points.append(subs)
        return tuple(points)

    if len(variables) == 2:
        x, y = variables
        if len(equalities) == 1:
            raise NotImplementedError(
                "a single plane equation is generally one-dimensional, not zero-dimensional"
            )
        try:
            raw = sp.solve(equalities, tuple(variables), dict=True)
        except _RECOVERABLE_ERRORS as exc:
            raise NotImplementedError(
                "could not solve the zero-dimensional equality system"
            ) from exc
        points = []
        seen: set[str] = set()
        for sol in raw:
            if not all(var in sol for var in variables):
                continue
            subs = {var: sp.simplify(sol[var]) for var in variables}
            if any(bool(sp.im(sp.N(val, 50)) != 0) for val in subs.values()):
                continue
            if _truth_at(condition, subs):
                key = tuple(sp.sstr(subs[var]) for var in variables)
                if str(key) not in seen:
                    seen.add(str(key))
                    points.append(subs)
        return tuple(points)

    raise NotImplementedError(
        "zero-dimensional integration currently supports one or two ambient variables"
    )


def _integrate_zero_dimensional(
    integrand: sp.Expr, condition: sp.Expr, variables: Sequence[sp.Symbol]
) -> sp.Expr:
    total = sp.Integer(0)
    for point in _zero_dimensional_points(condition, variables):
        total += sp.simplify(integrand.subs(point))
    return sp.simplify(total)


def _circle_radius_squared(condition: sp.Expr, x: sp.Symbol, y: sp.Symbol) -> sp.Expr | None:
    """Recognize x**2 + y**2 == r**2 centered at the origin."""

    try:
        atoms = _atoms(condition)
    except _RECOVERABLE_ERRORS:
        return None
    radius_sq: sp.Expr | None = None
    other_atoms: list[sp.Expr] = []
    for atom in atoms:
        if isinstance(atom, Equality):
            expr, _ = _relation_parts(atom)
            try:
                poly = sp.Poly(expr, x, y)
            except _RECOVERABLE_ERRORS:
                return None
            coeff_x2 = poly.coeff_monomial(x**2)
            coeff_y2 = poly.coeff_monomial(y**2)
            if (
                coeff_x2 != 0
                and sp.simplify(coeff_x2 - coeff_y2) == 0
                and all(monom in {(2, 0), (0, 2), (0, 0)} for monom in poly.monoms())
            ):
                candidate = sp.simplify(-poly.coeff_monomial(1) / coeff_x2)
                radius_sq = candidate if radius_sq is None else radius_sq
            else:
                other_atoms.append(atom)
        else:
            other_atoms.append(atom)
    if radius_sq is None or other_atoms:
        return None
    if bool(radius_sq < 0):
        return sp.Integer(0)
    return sp.simplify(radius_sq)


def _integrate_circle_intrinsic(
    integrand: sp.Expr, condition: sp.Expr, x: sp.Symbol, y: sp.Symbol
) -> sp.Expr | None:
    radius_sq = _circle_radius_squared(condition, x, y)
    if radius_sq is None:
        return None
    if radius_sq == 0:
        return sp.simplify(integrand.subs({x: 0, y: 0}))
    radius = sp.sqrt(radius_sq)
    theta = sp.Symbol("theta", real=True)
    expr = sp.simplify(
        radius * integrand.subs({x: radius * sp.cos(theta), y: radius * sp.sin(theta)})
    )
    value = sp.integrate(expr, (theta, 0, 2 * sp.pi))
    if isinstance(value, sp.Integral) or value.has(sp.Integral):
        raise NotImplementedError("SymPy could not evaluate the circle intrinsic integral")
    return sp.simplify(value)


def _graph_curve_data(
    condition: sp.Expr, x: sp.Symbol, y: sp.Symbol
) -> tuple[sp.Expr, tuple[tuple[sp.Expr, sp.Expr], ...]] | None:
    """Recognize a plane graph y = g(x) together with x-conditions."""

    try:
        atoms = _atoms(condition)
    except _RECOVERABLE_ERRORS:
        return None
    graph: sp.Expr | None = None
    x_conditions: list[sp.Expr] = []
    for atom in atoms:
        if isinstance(atom, Equality):
            expr, _ = _relation_parts(atom)
            if y not in expr.free_symbols:
                x_conditions.append(sp.Eq(expr, 0))
                continue
            try:
                poly_y = sp.Poly(expr, y)
            except _RECOVERABLE_ERRORS:
                return None
            if poly_y.degree() != 1:
                return None
            coeff = sp.simplify(poly_y.coeff_monomial(y))
            rest = sp.simplify(poly_y.as_expr() - coeff * y)
            candidate = sp.simplify(-rest / coeff)
            if candidate.free_symbols - {x}:
                return None
            graph = candidate if graph is None else graph
            if sp.simplify(graph - candidate) != 0:
                return None
        elif getattr(atom, "is_Relational", False):
            expr, op = _relation_parts(atom)
            if y in expr.free_symbols:
                # A first implementation avoids inequalities along graph curves.
                return None
            x_conditions.append(make_zero_relation(expr, op))
        else:
            return None
    if graph is None:
        return None
    x_condition = sp.And(*x_conditions) if x_conditions else sp.true
    intervals = _one_dimensional_intervals(x_condition, x, None)
    return graph, intervals


def _integrate_graph_curve_intrinsic(
    integrand: sp.Expr, condition: sp.Expr, x: sp.Symbol, y: sp.Symbol
) -> sp.Expr | None:
    data = _graph_curve_data(condition, x, y)
    if data is None:
        return None
    graph, intervals = data
    if not intervals:
        return sp.Integer(0)
    speed = sp.sqrt(1 + sp.diff(graph, x) ** 2)
    expr = sp.simplify(integrand.subs(y, graph) * speed)
    total = sp.Integer(0)
    for lo, hi in intervals:
        if lo == -sp.oo or hi == sp.oo:
            raise NotImplementedError("unbounded intrinsic graph integrals are not yet supported")
        value = sp.integrate(expr, (x, lo, hi))
        if isinstance(value, sp.Integral) or value.has(sp.Integral):
            raise NotImplementedError("SymPy could not evaluate the intrinsic graph integral")
        total += value
    return sp.simplify(total)


def _integrate_intrinsic_dimension(
    integrand: sp.Expr,
    condition: sp.Expr,
    variables: Sequence[sp.Symbol],
    measure_dimension: int,
    *,
    bounds: Mapping[sp.Symbol, tuple[sp.Expr, sp.Expr]],
) -> tuple[sp.Expr, str]:
    """Evaluate an intrinsic-dimensional integral with certified geometry."""

    ambient = len(variables)
    if measure_dimension == ambient:
        reduced = reduce_region_integral(integrand, condition, variables, bounds=bounds)
        assert isinstance(reduced, ReducedRegionIntegral)
        value, _, _ = _evaluate_reduced_integral(reduced, method="symbolic", precision=50)
        return value, reduced.method
    if measure_dimension == 0:
        return _integrate_zero_dimensional(
            integrand, condition, variables
        ), "zero_dimensional_counting_measure"
    if ambient == 2 and measure_dimension == 1:
        x, y = variables
        circle = _integrate_circle_intrinsic(integrand, condition, x, y)
        if circle is not None:
            return circle, "circle_intrinsic_length_measure"
        graph = _integrate_graph_curve_intrinsic(integrand, condition, x, y)
        if graph is not None:
            return graph, "graph_curve_intrinsic_length_measure"
    # Section cells use the induced Hausdorff metric, so they cannot be treated
    # as zero-width limits of ambient Lebesgue integrals.
    try:
        from .cad.cells import extract_cylindrical_solution
        from .cad.integration import intrinsic_solution_integrals

        solution = extract_cylindrical_solution(condition, variables, selected_only=True)
        pieces = intrinsic_solution_integrals(
            solution, integrand, dimension=measure_dimension, evaluate=False, require_verified=True
        )
        if pieces:
            values = []
            for piece in pieces:
                value = piece.doit()
                if isinstance(value, sp.Integral) or getattr(value, "has", lambda *_: False)(
                    sp.Integral
                ):
                    raise NotImplementedError(
                        "SymPy could not evaluate an intrinsic CAD-cell integral"
                    )
                values.append(value)
            return sp.simplify(sum(values, sp.Integer(0))), "cylindrical_solution_intrinsic_measure"
    except NotImplementedError:
        raise
    except _RECOVERABLE_ERRORS:
        pass
    raise NotImplementedError(
        "intrinsic integration supports finite point sets and verified cylindrical CAD graph cells; "
        "unsupported singular/non-graph strata fail conservatively"
    )


def _reduce_cylindrical_solution_cells_nd(
    integrand: sp.Expr,
    condition: sp.Expr,
    variables: Sequence[sp.Symbol],
) -> tuple[RegionIntegralPiece, ...] | None:
    """Reduce arbitrary-dimensional full CAD cells to nested iterated integrals.

    This exploits the public cylindrical-solution representation. It is exact
    for full-dimensional cells whose coordinate bounds can be expressed by the
    current CAD extractor. Lower-dimensional cells are ignored for ambient
    Lebesgue integration.
    """

    try:
        from .cad.cells import extract_cylindrical_solution, extract_explicit_cylindrical_solution
        from .cad.integration import full_dimensional_cell_integral

        cyl = extract_explicit_cylindrical_solution(condition, variables)
        if cyl is None:
            cyl = extract_cylindrical_solution(condition, variables, selected_only=True)
    except _RECOVERABLE_ERRORS:
        return None
    pieces: list[RegionIntegralPiece] = []
    n = len(tuple(variables))
    for cell in getattr(cyl, "cells", ()):
        if getattr(cell, "dimension", None) != n:
            continue
        try:
            adapted = full_dimensional_cell_integral(cell, integrand, require_verified=True)
        except _RECOVERABLE_ERRORS:
            continue
        limits = adapted.limits
        if len(limits) != n:
            continue
        pieces.append(
            RegionIntegralPiece(
                integrand=integrand,
                limits=limits,
                method="cylindrical_solution_cell_iterated_integral",
                diagnostics={
                    "cell_index": getattr(cell, "index", None),
                    "cell_dimension": getattr(cell, "dimension", None),
                    "typed_bounds_verified": adapted.certified_bounds,
                },
            )
        )
    return tuple(pieces) if pieces else None


def reduce_region_integral(
    integrand: object,
    condition: object,
    variables: Sequence[sp.Symbol | str],
    *,
    bounds: Sequence[tuple[sp.Symbol | str, object, object]]
    | Mapping[sp.Symbol | str, tuple[object, object]]
    | None = None,
    return_integrals: bool = False,
) -> ReducedRegionIntegral | tuple[sp.Integral, ...]:
    """Reduce a supported region integral to explicit iterated integrals.

    This is the structural layer used by ``integrate_over_region``. It does
    not call ``sympy.integrate`` unless callers later ask to evaluate the
    returned pieces. The reducer uses standard-region shortcuts, exact
    one-dimensional cells, radial/ellipse/simplex transformations, explicit
    vertical slices, and certified CAD-extracted cylindrical bounds where
    available. Unsupported or uncertified reductions are declined.
    """

    formula = _normalize_formula(condition)
    expr = sp.sympify(integrand)
    vars_ = _normalize_variables(variables, formula, expr)
    bound_map = _normalize_bounds(bounds, vars_)

    method = ""
    if len(vars_) == 1:
        pieces = _reduce_1d(expr, formula, vars_[0], bound_map.get(vars_[0]))
        method = "one_dimensional_cell_integration"
    elif len(vars_) == 2:
        x, y = vars_
        pieces = _reduce_box(expr, formula, vars_, bound_map)
        if pieces is not None:
            method = "axis_aligned_box_iterated_integral"
        if pieces is None and not bound_map:
            pieces = _reduce_standard_simplex_2d(expr, formula, x, y)
            if pieces is not None:
                method = "unit_simplex_iterated_integral"
        if pieces is None and not bound_map:
            pieces = _reduce_radial_vertical_pieces(expr, formula, x, y)
            if pieces is not None:
                method = "radial_region_as_signed_vertical_slices"
        if pieces is None and not bound_map:
            pieces = _reduce_axis_aligned_ellipse(expr, formula, x, y)
            if pieces is not None:
                method = "axis_aligned_ellipse_affine_unit_disk"
        if pieces is None:
            pieces = _reduce_vertical_slice(expr, formula, x, y, bound_map)
            if pieces is not None:
                method = "vertical_slice_iterated_integral"
        if pieces is None and not bound_map:
            pieces = _reduce_cylindrical_vertical_bounds_2d(expr, formula, x, y)
            if pieces is not None:
                has_nonlinear_stack = any(
                    sp.Poly(atom.lhs - atom.rhs, y).degree() > 1
                    for atom in getattr(formula, "args", (formula,))
                    if getattr(atom, "is_Relational", False)
                )
                method = (
                    "complete_cad_vertical_bounds_2d"
                    if has_nonlinear_stack
                    else "cylindrical_formula_vertical_bounds_2d"
                )
        if pieces is None and not bound_map:
            pieces = _reduce_cad_extracted_vertical_bounds_2d(expr, formula, x, y)
            if pieces is not None:
                method = "complete_cad_vertical_bounds_2d"
        if pieces is None:
            raise NotImplementedError(
                "reduce_region_integral currently supports exact 1D sets, axis-aligned boxes, the 2D unit simplex, "
                "2D axis-aligned ellipses, 2D origin-centered radial regions, simple 2D vertical-slice regions, "
                "and supported 2D cylindrical formulas"
            )
    else:
        pieces = _reduce_box(expr, formula, vars_, bound_map)
        if pieces is not None:
            method = "axis_aligned_box_iterated_integral"
        if pieces is None and not bound_map:
            pieces = _reduce_cylindrical_solution_cells_nd(expr, formula, vars_)
            if pieces is not None:
                method = "cylindrical_solution_cell_iterated_integral_nd"
        if pieces is None:
            raise NotImplementedError(
                "reduce_region_integral currently supports higher-dimensional axis-aligned boxes and "
                "full-dimensional cells exposed by the cylindrical CAD solution extractor"
            )
    reduced = ReducedRegionIntegral(
        integrand=expr,
        condition=formula,
        variables=vars_,
        pieces=tuple(pieces),
        method=method,
        diagnostics={"bounds": bound_map},
    )
    return reduced.as_integrals() if return_integrals else reduced


def integrate_over_region(
    integrand: object,
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
) -> sp.Expr | RegionIntegralResult:
    """Integrate ``integrand`` over a supported semialgebraic region.

    The region is first reduced to explicit iterated-integral pieces using
    ``reduce_region_integral``. The ``method`` option controls evaluation:

    ``method="symbolic"``
        Default. Evaluate exactly with SymPy and raise ``NotImplementedError``
        if any reduced piece remains unevaluated.

    ``method="numeric"``
        Numerically evaluate the reduced iterated integrals with ``evalf``.
        The returned value is approximate and ``RegionIntegralResult.exact`` is
        false when ``return_result=True``.

    ``method="auto"``
        Try symbolic evaluation first, then fall back to numerical evaluation.

    By default, integrals are with respect to ambient Lebesgue measure in
    ``variables``. ``measure_dimension="intrinsic"`` or an integer dimension
    enables exact intrinsic integration on certified regular CAD strata, along
    with specialized finite-point and standard-curve cases. Singular strata are
    kept explicit and are never silently treated as regular manifolds.
    """

    expr = sp.sympify(integrand)
    if isinstance(condition, StandardRegion):
        vars_ = _normalize_variables(variables, expr)
        if bounds is not None:
            raise NotImplementedError(
                "extra bounds are not yet supported for explicit StandardRegion objects"
            )
        if measure_dimension not in (None, "ambient", "intrinsic", "top", condition.dimension()):
            raise NotImplementedError(
                "explicit StandardRegion integration uses the region natural measure dimension"
            )
        value = integrate_over_standard_region(
            expr, condition, vars_, method=method, precision=precision
        )
        result = RegionIntegralResult(
            value=value,
            integrand=expr,
            condition=sp.Symbol(type(condition).__name__),
            variables=vars_,
            method=f"standard_region:{type(condition).__name__}",
            diagnostics={"requested_method": method, "precision": precision, "region": condition},
            exact=(method != "numeric"),
            evaluated=True,
            error_estimate=None,
        )
        return result if return_result else result.value

    formula = _normalize_formula(condition)
    vars_ = _normalize_variables(variables, formula, expr)
    bound_map = _normalize_bounds(bounds, vars_)
    dim = _normalize_measure_dimension(measure_dimension, len(vars_), formula)

    if dim != len(vars_):
        if method == "numeric":
            raise NotImplementedError("numeric intrinsic integration is not implemented yet")
        value, intrinsic_method = _integrate_intrinsic_dimension(
            expr, formula, vars_, dim, bounds=bound_map
        )
        result = RegionIntegralResult(
            value=value,
            integrand=expr,
            condition=formula,
            variables=vars_,
            method=intrinsic_method,
            diagnostics={
                "bounds": bound_map,
                "requested_method": method,
                "evaluation_method": "symbolic",
                "precision": precision,
                "measure_dimension": dim,
                "measure_dimension_request": measure_dimension,
            },
            exact=True,
            evaluated=True,
            error_estimate=None,
        )
        return result if return_result else result.value

    if len(vars_) == 2 and not bound_map and method in {"symbolic", "auto"}:
        x, y = vars_
        fast_cases = (
            (_integrate_radial_polynomial(expr, formula, x, y), "radial_polynomial_moments"),
            (
                _integrate_axis_aligned_ellipse_polynomial(expr, formula, x, y),
                "axis_aligned_ellipse_polynomial_moments",
            ),
        )
        for fast_value, fast_method in fast_cases:
            if fast_value is not None:
                reduced = reduce_region_integral(expr, formula, vars_, bounds=bound_map)
                assert isinstance(reduced, ReducedRegionIntegral)
                result = RegionIntegralResult(
                    value=fast_value,
                    integrand=expr,
                    condition=formula,
                    variables=vars_,
                    method=reduced.method,
                    diagnostics={
                        **(reduced.diagnostics or {}),
                        "pieces": reduced.pieces,
                        "requested_method": method,
                        "evaluation_method": "symbolic",
                        "accelerated_by": fast_method,
                        "precision": precision,
                        "measure_dimension": dim,
                        "measure_dimension_request": measure_dimension,
                    },
                    exact=True,
                    evaluated=True,
                    error_estimate=None,
                )
                return result if return_result else result.value

    reduced = reduce_region_integral(expr, formula, vars_, bounds=bound_map)
    assert isinstance(reduced, ReducedRegionIntegral)
    value, exact, evaluation_method = _evaluate_reduced_integral(
        reduced, method=method, precision=precision
    )
    result = RegionIntegralResult(
        value=value,
        integrand=reduced.integrand,
        condition=reduced.condition,
        variables=reduced.variables,
        method=reduced.method,
        diagnostics={
            **(reduced.diagnostics or {}),
            "pieces": reduced.pieces,
            "requested_method": method,
            "evaluation_method": evaluation_method,
            "precision": precision,
            "measure_dimension": dim,
            "measure_dimension_request": measure_dimension,
        },
        exact=exact,
        evaluated=True,
        error_estimate=None,
    )
    return result if return_result else result.value


__all__ = [
    "RegionIntegralResult",
    "RegionIntegralPiece",
    "ReducedRegionIntegral",
    "reduce_region_integral",
    "integrate_over_region",
]
