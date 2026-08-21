from __future__ import annotations

from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass

import sympy as sp
from sympy.logic.boolalg import And as SymAnd
from sympy.logic.boolalg import BooleanFalse, BooleanTrue
from sympy.logic.boolalg import Not as SymNot
from sympy.logic.boolalg import Or as SymOr

from .interval_decomposition import (
    one_dimensional_intervals as _shared_intervals,
)
from .relations import make_zero_relation
from .relations import split_relation as _relation_parts


@dataclass(frozen=True)
class MeasureResult:
    """Exact measure result for a semialgebraic set.

    The ``value`` field is the Lebesgue measure in the ambient variables used
    by the call. ``method`` records the reconstruction strategy that produced
    the result.
    """

    value: sp.Expr
    variables: tuple[sp.Symbol, ...]
    method: str
    diagnostics: Mapping[str, object] | None = None


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
    if isinstance(condition, SymOr) or isinstance(condition, SymNot):
        raise NotImplementedError(
            "2D measure currently supports conjunctions; 1D supports Boolean combinations"
        )
    if getattr(condition, "is_Relational", False):
        return (condition,)
    raise TypeError(f"unsupported formula expression: {condition!r}")


def _one_dimensional_intervals(
    condition: sp.Expr,
    variable: sp.Symbol,
    bound: tuple[sp.Expr, sp.Expr] | None,
) -> tuple[tuple[sp.Expr, sp.Expr], ...]:
    """Return exact true intervals for a one-dimensional measure problem."""

    return _shared_intervals(
        condition,
        variable,
        bound,
        extra_symbol_error=(
            "1D measure received a formula containing symbols outside the integration variable"
        ),
    )


def _interval_measure(intervals: Iterable[tuple[sp.Expr, sp.Expr]]) -> sp.Expr:
    total = sp.Integer(0)
    for lo, hi in intervals:
        if lo == -sp.oo or hi == sp.oo:
            return sp.oo
        total += sp.simplify(hi - lo)
    return sp.simplify(total)


def _measure_1d(
    condition: sp.Expr, variable: sp.Symbol, bound: tuple[sp.Expr, sp.Expr] | None
) -> sp.Expr:
    return _interval_measure(_one_dimensional_intervals(condition, variable, bound))


def _radial_measure(condition: sp.Expr, x: sp.Symbol, y: sp.Symbol) -> sp.Expr | None:
    try:
        atoms = _atoms(condition)
    except NotImplementedError:
        return None
    lower = sp.Integer(0)
    upper = sp.oo
    found = False
    for atom in atoms:
        if atom is sp.false:
            return sp.Integer(0)
        if not getattr(atom, "is_Relational", False):
            return None
        expr, op = _relation_parts(atom)
        try:
            poly = sp.Poly(expr, x, y)
        except Exception:
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
            return sp.Integer(0)
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
        return sp.oo
    if bool(upper <= lower):
        return sp.Integer(0)
    if bool(upper <= 0):
        return sp.Integer(0)
    if bool(lower < 0):
        lower = sp.Integer(0)
    return sp.simplify(sp.pi * (upper - lower))


def _vertical_slice_measure(
    condition: sp.Expr,
    x: sp.Symbol,
    y: sp.Symbol,
    bounds: Mapping[sp.Symbol, tuple[sp.Expr, sp.Expr]],
) -> sp.Expr | None:
    """Compute planar measure by reducing a supported region to vertical slice bounds."""
    try:
        atoms = _atoms(condition)
    except NotImplementedError:
        return None
    lower_bounds: list[sp.Expr] = []
    upper_bounds: list[sp.Expr] = []
    x_conditions: list[sp.Expr] = []
    for atom in atoms:
        if atom is sp.false:
            return sp.Integer(0)
        expr, op = _relation_parts(atom)
        if op == "==":
            # A nontrivial equality cuts dimension in R^2, so its planar measure is zero.
            if sp.simplify(expr) != 0:
                return sp.Integer(0)
            continue
        if op == "!=":
            # Removing a polynomial hypersurface does not change planar measure.
            continue
        if y not in expr.free_symbols:
            x_conditions.append(make_zero_relation(expr, op))
            continue
        try:
            poly_y = sp.Poly(expr, y)
        except Exception:
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
            # Coefficient depending on x would require sign-splitting by CAD.
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
    x_condition = sp.And(*x_conditions) if x_conditions else sp.true
    intervals = _one_dimensional_intervals(x_condition, x, None)
    if not intervals:
        return sp.Integer(0)
    height = sp.simplify(upper - lower)
    total = sp.Integer(0)
    for lo, hi in intervals:
        if lo == -sp.oo or hi == sp.oo:
            return sp.oo
        total += sp.integrate(height, (x, lo, hi))
    return sp.simplify(total)


def semialgebraic_measure(
    condition: object,
    variables: Sequence[sp.Symbol | str],
    *,
    bounds: Sequence[tuple[sp.Symbol | str, object, object]]
    | Mapping[sp.Symbol | str, tuple[object, object]]
    | None = None,
    measure_dimension: object = "ambient",
    return_result: bool = False,
) -> sp.Expr | MeasureResult:
    """Return the exact measure of a supported semialgebraic set.

    The measure implementation delegates to the structural region-integral
    reducer with integrand ``1``. This keeps ``semialgebraic_measure`` aligned
    with all standard-shape recognition supported by ``integrate_over_region``:
    one-dimensional cell intervals, axis-aligned boxes, the 2D unit simplex,
    axis-aligned ellipses, origin-centered disks/annuli, simple vertical slices,
    and the first intrinsic-dimensional cases.
    """

    from .region_integrate import integrate_over_region

    result = integrate_over_region(
        1,
        condition,
        variables,
        bounds=bounds,
        measure_dimension=measure_dimension,
        return_result=True,
    )
    method_name = result.method
    if method_name == "one_dimensional_cell_integration":
        method_name = "one_dimensional_cell_sampling"
    measure_result = MeasureResult(
        result.value,
        result.variables,
        method_name,
        {**(result.diagnostics or {}), "delegated_to": "integrate_over_region"},
    )
    return measure_result if return_result else measure_result.value


__all__ = ["MeasureResult", "semialgebraic_measure"]
