from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass, field

import sympy as sp
from sympy import S
from sympy.core.relational import Relational

from .periodic import compute_periodic_window, detect_real_period, recon_periodic_represent
from .semantics import ResultSemantics

_RECOVERABLE_ERRORS = (
    ArithmeticError,
    TypeError,
    ValueError,
    NotImplementedError,
    RuntimeError,
    sp.PolynomialError,
)


@dataclass(frozen=True)
class CertifiedIntervalRoot:
    left: sp.Expr
    right: sp.Expr
    midpoint: sp.Expr
    residual_midpoint: sp.Expr
    sign_change_certified: bool = False


@dataclass(frozen=True)
class RootIsolationResult:
    variable: sp.Symbol
    equation: sp.Expr
    domain: object
    roots: tuple[sp.Expr, ...] = ()
    representative_roots: tuple[sp.Expr, ...] = ()
    certified_intervals: tuple[CertifiedIntervalRoot, ...] = ()
    periodic_formula: sp.Expr | None = None
    intervals_of_zero: tuple[tuple[sp.Expr, sp.Expr], ...] = ()
    complete: bool = False
    method: str = "univariate_transcendental_isolation"
    result_semantics: ResultSemantics = ResultSemantics.UNKNOWN
    validity_window: tuple[sp.Expr, sp.Expr] | None = None
    metadata: dict = field(default_factory=dict)


@dataclass(frozen=True)
class SampledTruthDecomp:
    variable: sp.Symbol
    support_points: tuple[sp.Expr, ...]
    true_points: tuple[sp.Expr, ...]
    true_intervals: tuple[tuple[sp.Expr, sp.Expr], ...] = ()
    method: str = "sampled_truth_decomposition"
    result_semantics: ResultSemantics = ResultSemantics.WINDOW_APPROXIMATION
    validity_window: tuple[sp.Expr, sp.Expr] | None = None
    metadata: dict = field(default_factory=dict)


def _extract_equation(expr: sp.Expr) -> sp.Expr:
    if isinstance(expr, sp.Equality):
        return sp.simplify(expr.lhs - expr.rhs)
    if isinstance(expr, Relational) or isinstance(expr, sp.logic.boolalg.Boolean):
        raise TypeError(
            "Root isolation requires an equality or scalar residual, not a Boolean relation"
        )
    return sp.simplify(expr)


def _numeric_function(expr: sp.Expr, variable: sp.Symbol):
    try:
        return sp.lambdify(variable, expr, "mpmath")
    except _RECOVERABLE_ERRORS:
        return None


def _real_eval(func, x: float):
    try:
        value = complex(func(x))
    except _RECOVERABLE_ERRORS:
        return None
    if abs(value.imag) > 1e-8:
        return None
    return float(value.real)


def _bisect_sign_change(func, left: float, right: float, *, tol: float = 1e-10, max_iter: int = 80):
    fl = _real_eval(func, left)
    fr = _real_eval(func, right)
    if fl is None or fr is None:
        return None
    if fl == 0.0:
        mid = left
    elif fr == 0.0:
        mid = right
    elif fl * fr > 0:
        return None
    else:
        lo, hi = left, right
        vlo = fl
        for _ in range(max_iter):
            mid = 0.5 * (lo + hi)
            vmid = _real_eval(func, mid)
            if vmid is None:
                break
            if abs(vmid) <= tol or abs(hi - lo) <= tol:
                return CertifiedIntervalRoot(
                    left=sp.nsimplify(lo),
                    right=sp.nsimplify(hi),
                    midpoint=sp.nsimplify(mid),
                    residual_midpoint=sp.nsimplify(abs(vmid)),
                    sign_change_certified=False,
                )
            if vlo * vmid <= 0:
                hi = mid
            else:
                lo, vlo = mid, vmid
        mid = 0.5 * (lo + hi)
    vmid = _real_eval(func, mid)
    if vmid is None:
        return None
    return CertifiedIntervalRoot(
        left=sp.nsimplify(left),
        right=sp.nsimplify(right),
        midpoint=sp.nsimplify(mid),
        residual_midpoint=sp.nsimplify(abs(vmid)),
        sign_change_certified=False,
    )


def _certified_brackets(
    expr: sp.Expr, variable: sp.Symbol, lo: float, hi: float, *, samples: int = 160
):
    func = _numeric_function(expr, variable)
    if func is None:
        return ()
    xs = [lo + (hi - lo) * i / samples for i in range(samples + 1)]
    intervals = []
    prev_x = xs[0]
    prev_v = _real_eval(func, prev_x)
    if prev_v is not None and abs(prev_v) <= 1e-8:
        intervals.append(
            CertifiedIntervalRoot(
                sp.nsimplify(prev_x),
                sp.nsimplify(prev_x),
                sp.nsimplify(prev_x),
                sp.Integer(0),
                False,
            )
        )
    for x in xs[1:]:
        cur_v = _real_eval(func, x)
        if prev_v is not None and cur_v is not None:
            if cur_v == 0.0:
                intervals.append(
                    CertifiedIntervalRoot(
                        sp.nsimplify(x), sp.nsimplify(x), sp.nsimplify(x), sp.Integer(0), False
                    )
                )
            elif prev_v * cur_v < 0:
                bracket = _bisect_sign_change(func, prev_x, x)
                if bracket is not None:
                    intervals.append(bracket)
        prev_x, prev_v = x, cur_v
    # de-duplicate by midpoint
    uniq = []
    seen = set()
    for ci in intervals:
        key = sp.srepr(sp.nsimplify(ci.midpoint))
        if key not in seen:
            seen.add(key)
            uniq.append(ci)
    return tuple(uniq)


def _support_points_intvs(intervals: Sequence[CertifiedIntervalRoot], lo: float, hi: float):
    mids = [float(sp.N(ci.midpoint)) for ci in intervals if ci.left != ci.right]
    boundaries = [lo] + mids + [hi]
    support = []
    for a, b in zip(boundaries[:-1], boundaries[1:], strict=True):
        if b - a > 1e-8:
            support.append(sp.nsimplify((a + b) / 2.0))
    return tuple(support)


def isolate_univar_roots(
    expr: sp.Expr,
    variable: sp.Symbol,
    *,
    domain: object = S.Reals,
) -> RootIsolationResult:
    """Isolate real roots of a supported univariate transcendental expression."""
    equation = _extract_equation(expr)
    try:
        solset = sp.solveset(equation, variable, domain=domain)
    except _RECOVERABLE_ERRORS:
        solset = None

    if isinstance(solset, sp.FiniteSet):
        roots = tuple(sorted(solset, key=sp.default_sort_key))
        return RootIsolationResult(
            variable=variable,
            equation=equation,
            domain=domain,
            roots=roots,
            representative_roots=roots,
            complete=True,
            method="solveset_finite",
            result_semantics=ResultSemantics.EXACT,
        )

    period = detect_real_period(equation, variable) if domain == S.Reals else None
    if isinstance(solset, sp.ImageSet) and period is not None:
        window = compute_periodic_window(equation, variable)
        reps = []
        representatives_complete = False
        if window is not None:
            try:
                rep_set = sp.solveset(
                    equation, variable, domain=sp.Interval(window.lower_bound, window.upper_bound)
                )
                if isinstance(rep_set, sp.FiniteSet):
                    reps = tuple(sorted(rep_set, key=sp.default_sort_key))
                    representatives_complete = True
            except _RECOVERABLE_ERRORS:
                reps = []
        formula = (
            recon_periodic_represent(variable, reps, period)
            if reps
            else sp.ConditionSet(variable, sp.Eq(equation, 0), domain)
        )
        return RootIsolationResult(
            variable=variable,
            equation=equation,
            domain=domain,
            representative_roots=tuple(reps),
            periodic_formula=formula,
            complete=representatives_complete,
            method="periodic_solveset",
            result_semantics=ResultSemantics.EXACT if representatives_complete else "subset",
            metadata={"period": period},
        )

    if domain in (S.Reals, S.Complexes):
        try:
            poly = sp.Poly(equation, variable)
            if poly.is_univariate:
                all_roots = tuple(sorted(poly.all_roots(), key=sp.default_sort_key))
                roots = (
                    tuple(root for root in all_roots if root.is_real is True)
                    if domain == S.Reals
                    else all_roots
                )
                return RootIsolationResult(
                    variable=variable,
                    equation=equation,
                    domain=domain,
                    roots=roots,
                    representative_roots=roots,
                    complete=True,
                    method="poly_all_roots",
                    result_semantics=ResultSemantics.EXACT,
                )
        except (ArithmeticError, TypeError, ValueError, NotImplementedError, sp.PolynomialError):
            pass

    if domain == S.Reals:
        certified = _certified_brackets(equation, variable, -10.0, 10.0)
        if certified:
            roots = tuple(ci.midpoint for ci in certified)
            return RootIsolationResult(
                variable=variable,
                equation=equation,
                domain=domain,
                roots=roots,
                representative_roots=roots,
                certified_intervals=certified,
                complete=False,
                method="numerical_sign_change_isolation",
                result_semantics=ResultSemantics.WINDOW_SUBSET,
                validity_window=(sp.Integer(-10), sp.Integer(10)),
                metadata={"window": (-10.0, 10.0), "interval_count": len(certified)},
            )

    return RootIsolationResult(
        variable=variable,
        equation=equation,
        domain=domain,
        complete=False,
        method="unsolved_univariate",
        result_semantics=ResultSemantics.UNKNOWN,
        metadata={"solveset": solset},
    )


def evaluate_form_points(
    formula: sp.Expr, variable: sp.Symbol, support_points: Sequence[sp.Expr]
) -> SampledTruthDecomp:
    support_points = tuple(support_points)
    truths = []
    for pt in support_points:
        try:
            if bool(sp.simplify(formula.subs(variable, pt))):
                truths.append(pt)
        except _RECOVERABLE_ERRORS:
            continue
    return SampledTruthDecomp(
        variable=variable, support_points=support_points, true_points=tuple(truths)
    )


def decomp_univar_inequality(
    formula: sp.Expr,
    variable: sp.Symbol,
    *,
    domain: object = S.Reals,
    search_window: tuple[float, float] = (-10.0, 10.0),
) -> SampledTruthDecomp:
    """Decompose a supported univariate transcendental inequality into sign intervals."""
    if domain != S.Reals:
        return SampledTruthDecomp(
            variable=variable,
            support_points=(),
            true_points=(),
            true_intervals=(),
            method="unsupported_complex_inequality",
        )

    atoms = list(formula.args) if isinstance(formula, sp.And) else [formula]
    boundary_exprs = []
    for atom in atoms:
        if isinstance(
            atom,
            (sp.Equality, sp.StrictLessThan, sp.LessThan, sp.StrictGreaterThan, sp.GreaterThan),
        ):
            boundary_exprs.append(sp.simplify(atom.lhs - atom.rhs))

    certified_roots = []
    for expr in boundary_exprs:
        for ci in _certified_brackets(expr, variable, search_window[0], search_window[1]):
            certified_roots.append(ci)
    # unique by midpoint
    uniq = []
    seen = set()
    for ci in certified_roots:
        key = sp.srepr(sp.nsimplify(ci.midpoint))
        if key not in seen:
            seen.add(key)
            uniq.append(ci)
    certified_roots = tuple(sorted(uniq, key=lambda ci: float(sp.N(ci.midpoint))))
    support = _support_points_intvs(certified_roots, search_window[0], search_window[1])

    true_points = []
    true_intervals = []
    for pt in support:
        try:
            truth = bool(sp.simplify(formula.subs(variable, pt)))
        except _RECOVERABLE_ERRORS:
            truth = False
        if truth:
            true_points.append(pt)
    mids = [float(sp.N(ci.midpoint)) for ci in certified_roots]
    boundaries = [search_window[0]] + mids + [search_window[1]]
    for a, b in zip(boundaries[:-1], boundaries[1:], strict=True):
        if b - a <= 1e-8:
            continue
        mid = sp.nsimplify((a + b) / 2.0)
        try:
            truth = bool(sp.simplify(formula.subs(variable, mid)))
        except _RECOVERABLE_ERRORS:
            truth = False
        if truth:
            true_intervals.append((sp.nsimplify(a), sp.nsimplify(b)))

    return SampledTruthDecomp(
        variable=variable,
        support_points=tuple(support),
        true_points=tuple(true_points),
        true_intervals=tuple(true_intervals),
        method="numerical_interval_decomposition"
        if certified_roots
        else "empty_numerical_interval_decomposition",
        result_semantics=ResultSemantics.WINDOW_APPROXIMATION,
        validity_window=(sp.nsimplify(search_window[0]), sp.nsimplify(search_window[1])),
        metadata={"certified_roots": certified_roots, "search_window": search_window},
    )


__all__ = [
    "CertifiedIntervalRoot",
    "RootIsolationResult",
    "SampledTruthDecomp",
    "isolate_univar_roots",
    "evaluate_form_points",
    "decomp_univar_inequality",
]
