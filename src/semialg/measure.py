from __future__ import annotations

from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass

import sympy as sp
from sympy.core.relational import (
    Equality,
    GreaterThan,
    LessThan,
    StrictGreaterThan,
    StrictLessThan,
    Unequality,
)
from sympy.logic.boolalg import And as SymAnd
from sympy.logic.boolalg import BooleanFalse, BooleanTrue
from sympy.logic.boolalg import Not as SymNot
from sympy.logic.boolalg import Or as SymOr

from .formula import to_sympy


@dataclass(frozen=True)
class MeasureResult:
    """Exact measure result for a semialgebraic set.

    The ``value`` field is the Lebesgue measure in the ambient variables used
    by the call. ``method`` records the intentionally limited first-pass
    reconstruction strategy that produced the result.
    """

    value: sp.Expr
    variables: tuple[sp.Symbol, ...]
    method: str
    diagnostics: Mapping[str, object] | None = None


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


def _normalize_bounds(
    bounds: Sequence[tuple[sp.Symbol | str, object, object]]
    | Mapping[sp.Symbol | str, tuple[object, object]]
    | None,
    variables: Sequence[sp.Symbol],
) -> dict[sp.Symbol, tuple[sp.Expr, sp.Expr]]:
    if bounds is None:
        return {}
    out: dict[sp.Symbol, tuple[sp.Expr, sp.Expr]] = {}
    if isinstance(bounds, Mapping):
        items = [(key, value[0], value[1]) for key, value in bounds.items()]
    else:
        items = list(bounds)
    by_name = {var.name: var for var in variables}
    for raw_var, lo, hi in items:
        if isinstance(raw_var, str):
            var = by_name.get(raw_var, sp.Symbol(raw_var, real=True))
        else:
            var = raw_var
        out[var] = (sp.sympify(lo), sp.sympify(hi))
    return out


def _relation_parts(atom: sp.Expr) -> tuple[sp.Expr, str]:
    if isinstance(atom, Equality):
        return sp.expand(atom.lhs - atom.rhs), "=="
    if isinstance(atom, Unequality):
        return sp.expand(atom.lhs - atom.rhs), "!="
    if isinstance(atom, StrictLessThan):
        return sp.expand(atom.lhs - atom.rhs), "<"
    if isinstance(atom, LessThan):
        return sp.expand(atom.lhs - atom.rhs), "<="
    if isinstance(atom, StrictGreaterThan):
        return sp.expand(atom.lhs - atom.rhs), ">"
    if isinstance(atom, GreaterThan):
        return sp.expand(atom.lhs - atom.rhs), ">="
    raise TypeError(f"expected a relational atom, got {atom!r}")


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


def _truth_at(condition: sp.Expr, subs: Mapping[sp.Symbol, sp.Expr]) -> bool:
    value = condition.subs(subs)
    if value is sp.true or isinstance(value, BooleanTrue):
        return True
    if value is sp.false or isinstance(value, BooleanFalse):
        return False
    simplified = sp.simplify(value)
    if simplified is sp.true or isinstance(simplified, BooleanTrue):
        return True
    if simplified is sp.false or isinstance(simplified, BooleanFalse):
        return False
    numeric = simplified.evalf(80)
    if numeric is sp.true:
        return True
    if numeric is sp.false:
        return False
    return bool(numeric)


def _finite_real_roots(poly_expr: sp.Expr, variable: sp.Symbol) -> tuple[sp.Expr, ...]:
    poly = sp.Poly(poly_expr, variable)
    if poly.is_zero:
        return ()
    try:
        roots = sp.real_roots(poly.as_expr())
    except Exception:
        roots = [
            root
            for root in sp.nroots(poly.as_expr(), n=80, maxsteps=200)
            if abs(sp.im(root)) < sp.Rational(1, 10) ** 40
        ]
    distinct: list[sp.Expr] = []
    seen: set[str] = set()
    for root in roots:
        root_expr = sp.re(root) if not isinstance(root, sp.Expr) else root
        key = sp.sstr(root_expr)
        if key not in seen:
            distinct.append(root_expr)
            seen.add(key)
    return tuple(sorted(distinct, key=lambda z: float(sp.N(z, 50))))


def _relational_polynomials(condition: sp.Expr) -> tuple[sp.Expr, ...]:
    if condition is sp.true or isinstance(condition, BooleanTrue):
        return ()
    if condition is sp.false or isinstance(condition, BooleanFalse):
        return ()
    if getattr(condition, "is_Relational", False):
        expr, _ = _relation_parts(condition)
        return (expr,)
    if isinstance(condition, (SymAnd, SymOr)):
        out: list[sp.Expr] = []
        for arg in condition.args:
            out.extend(_relational_polynomials(arg))
        return tuple(out)
    if isinstance(condition, SymNot):
        return _relational_polynomials(condition.args[0])
    raise TypeError(f"unsupported formula expression: {condition!r}")


def _sample_between(left: sp.Expr, right: sp.Expr) -> sp.Expr:
    if left == -sp.oo and right == sp.oo:
        return sp.Integer(0)
    if left == -sp.oo:
        return sp.simplify(right - 1)
    if right == sp.oo:
        return sp.simplify(left + 1)
    return sp.simplify((left + right) / 2)


def _interval_measure(intervals: Iterable[tuple[sp.Expr, sp.Expr]]) -> sp.Expr:
    total = sp.Integer(0)
    for lo, hi in intervals:
        if lo == -sp.oo or hi == sp.oo:
            return sp.oo
        total += sp.simplify(hi - lo)
    return sp.simplify(total)


def _one_dimensional_intervals(
    condition: sp.Expr,
    variable: sp.Symbol,
    bound: tuple[sp.Expr, sp.Expr] | None,
) -> tuple[tuple[sp.Expr, sp.Expr], ...]:
    if condition is sp.false or isinstance(condition, BooleanFalse):
        return ()
    lo, hi = bound if bound is not None else (-sp.oo, sp.oo)
    roots: list[sp.Expr] = []
    if lo is not -sp.oo:
        roots.append(lo)
    if hi is not sp.oo:
        roots.append(hi)
    for poly in _relational_polynomials(condition):
        if variable not in poly.free_symbols and sp.simplify(poly) != 0:
            continue
        if not poly.free_symbols <= {variable}:
            raise ValueError(
                "1D measure received a formula containing symbols outside the integration variable"
            )
        roots.extend(_finite_real_roots(poly, variable))
    ordered: list[sp.Expr] = []
    seen: set[str] = set()
    for root in sorted(
        roots,
        key=lambda z: (
            float(sp.N(z, 50))
            if z not in (-sp.oo, sp.oo)
            else (-float("inf") if z is -sp.oo else float("inf"))
        ),
    ):
        if root == -sp.oo or root == sp.oo:
            continue
        if lo != -sp.oo and float(sp.N(root, 50)) < float(sp.N(lo, 50)):
            continue
        if hi != sp.oo and float(sp.N(root, 50)) > float(sp.N(hi, 50)):
            continue
        key = sp.sstr(root)
        if key not in seen:
            ordered.append(root)
            seen.add(key)
    cuts = [lo, *ordered, hi]
    intervals: list[tuple[sp.Expr, sp.Expr]] = []
    for left, right in zip(cuts, cuts[1:], strict=False):
        if left == right:
            continue
        sample = _sample_between(left, right)
        if _truth_at(condition, {variable: sample}):
            intervals.append((left, right))
    return tuple(intervals)


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


def _as_relation(expr: sp.Expr, op: str) -> sp.Expr:
    if op == "<":
        return expr < 0
    if op == "<=":
        return expr <= 0
    if op == ">":
        return expr > 0
    if op == ">=":
        return expr >= 0
    if op == "==":
        return sp.Eq(expr, 0)
    if op == "!=":
        return sp.Ne(expr, 0)
    raise ValueError(op)


def _vertical_slice_measure(
    condition: sp.Expr,
    x: sp.Symbol,
    y: sp.Symbol,
    bounds: Mapping[sp.Symbol, tuple[sp.Expr, sp.Expr]],
) -> sp.Expr | None:
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
            x_conditions.append(_as_relation(expr, op))
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
