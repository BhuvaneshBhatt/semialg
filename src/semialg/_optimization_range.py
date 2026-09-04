from __future__ import annotations

from collections.abc import Iterable, Mapping, Sequence

import sympy as sp
from sympy.core.sympify import SympifyError
from sympy.logic.boolalg import And as SymAnd
from sympy.logic.boolalg import Or as SymOr
from sympy.polys.polyerrors import GeneratorsNeeded, PolynomialError

from ._function_algebraization import (
    algebraic_back_substitute,
    exact_algebraize_function_problem,
)
from .context import with_computation_context
from .exact_arithmetic import compare_exact_reals
from .normalization import normalize_constraints, normalize_problem_variables
from .optimization_results import FunctionRangeResult, OptimizationResult
from .symbol_resolution import resolve_symbol

_EXPECTED_ERRORS = (
    TypeError,
    ValueError,
    ArithmeticError,
    NotImplementedError,
    SympifyError,
    PolynomialError,
    GeneratorsNeeded,
)
FormulaLike = sp.Expr | sp.logic.boolalg.Boolean | bool
_finite_compare = compare_exact_reals


def _range_formula(
    symbol: sp.Symbol, lo: sp.Expr, hi: sp.Expr, lo_attained: bool, hi_attained: bool
) -> sp.Expr:
    pieces: list[sp.Expr] = []
    if lo != -sp.oo:
        pieces.append(symbol >= lo if lo_attained else symbol > lo)
    if hi != sp.oo:
        pieces.append(symbol <= hi if hi_attained else symbol < hi)
    if not pieces:
        return sp.true
    return sp.And(*pieces)


def _simplify_range_condition(formula: sp.Expr, value_symbol: sp.Symbol) -> sp.Expr:
    if formula is sp.true or formula is sp.false:
        return formula
    if isinstance(formula, SymAnd):
        return sp.And(*(_simplify_range_condition(arg, value_symbol) for arg in formula.args))
    if isinstance(formula, SymOr):
        return sp.Or(*(_simplify_range_condition(arg, value_symbol) for arg in formula.args))
    if getattr(formula, "is_Relational", False) and formula.free_symbols <= {value_symbol}:
        return formula
    return formula


def _intervals_from_range_set(
    range_set: sp.Set,
) -> tuple[tuple[sp.Expr, sp.Expr, bool, bool], ...] | None:
    """Return intervals as ``(left, right, left_closed, right_closed)``.

    The helper intentionally accepts only one-dimensional set forms that have
    direct range-bound semantics: intervals, finite point sets, and finite
    unions of those. More exotic sets leave metadata unknown without changing
    the primary range formula.
    """

    if range_set is sp.S.EmptySet:
        return ()
    if range_set is sp.S.Reals or range_set is sp.S.UniversalSet:
        return ((-sp.oo, sp.oo, False, False),)
    if isinstance(range_set, sp.Interval):
        return (
            (range_set.start, range_set.end, not range_set.left_open, not range_set.right_open),
        )
    if isinstance(range_set, sp.FiniteSet):
        return tuple(
            (item, item, True, True) for item in sorted(range_set, key=sp.default_sort_key)
        )
    if isinstance(range_set, sp.Union):
        pieces: list[tuple[sp.Expr, sp.Expr, bool, bool]] = []
        for arg in range_set.args:
            intervals = _intervals_from_range_set(arg)
            if intervals is None:
                return None
            pieces.extend(intervals)
        return tuple(sorted(pieces, key=lambda item: sp.default_sort_key(item[0])))
    return None


def _compare_bounds_for_sort(a: sp.Expr, b: sp.Expr) -> int:
    if a == b:
        return 0
    if a == -sp.oo or b == sp.oo:
        return -1
    if a == sp.oo or b == -sp.oo:
        return 1
    return _finite_compare(a, b)


def _range_metadata_from_formula(
    formula: sp.Expr,
    value_symbol: sp.Symbol,
) -> tuple[sp.Expr | None, sp.Expr | None, bool | None, bool | None, bool | None, int | None]:
    """Extract interval metadata from a one-variable range condition.

    The range formula remains authoritative. These fields are summaries for
    interval-like answers and are intentionally ``None`` when extraction is not
    reliable.
    """

    from ._range_special_cases import _domain_set_from_condition

    range_set = _domain_set_from_condition(formula, value_symbol)
    if range_set is None:
        try:
            range_set = formula.as_set()
        except _EXPECTED_ERRORS:
            return None, None, None, None, None, None
    intervals = _intervals_from_range_set(range_set)
    if intervals is None:
        return None, None, None, None, None, None
    interval_count = len(intervals)
    if interval_count == 0:
        return None, None, False, False, True, 0
    left = intervals[0][0]
    right = intervals[0][1]
    left_closed = intervals[0][2]
    right_closed = intervals[0][3]
    for item in intervals[1:]:
        if _compare_bounds_for_sort(item[0], left) < 0:
            left, left_closed = item[0], item[2]
        if _compare_bounds_for_sort(right, item[1]) < 0:
            right, right_closed = item[1], item[3]
    lower_attained = bool(left_closed) if left not in (-sp.oo, sp.oo) else False
    upper_attained = bool(right_closed) if right not in (-sp.oo, sp.oo) else False
    is_interval = interval_count == 1
    return left, right, lower_attained, upper_attained, is_interval, interval_count


def _range_result_from_formula(
    expression: sp.Expr,
    formula: sp.Expr,
    value_symbol: sp.Symbol,
    variables: tuple[sp.Symbol, ...],
    method: str,
    diagnostics: Mapping[str, object],
) -> FunctionRangeResult:
    simplified_formula = _simplify_range_condition(formula, value_symbol)
    lo, hi, lo_attained, hi_attained, is_interval, interval_count = _range_metadata_from_formula(
        simplified_formula,
        value_symbol,
    )
    return FunctionRangeResult(
        expression,
        simplified_formula,
        value_symbol,
        variables,
        lo,
        hi,
        lo_attained,
        hi_attained,
        (),
        (),
        method,
        {
            **dict(diagnostics),
            "metadata_source": "range_formula_as_set" if interval_count is not None else "unknown",
        },
        is_interval,
        interval_count,
    )


def _range_via_optimization_bounds(
    expr: sp.Expr,
    condition: sp.Expr,
    vars_: tuple[sp.Symbol, ...],
    val_sym: sp.Symbol,
    domain: str,
) -> FunctionRangeResult:
    from .optimization import semialgebraic_maximize, semialgebraic_minimize

    minimum = semialgebraic_minimize(expr, condition, vars_, domain=domain, return_result=True)
    maximum = semialgebraic_maximize(expr, condition, vars_, domain=domain, return_result=True)
    if not isinstance(minimum, OptimizationResult) or not isinstance(maximum, OptimizationResult):
        raise TypeError("range optimization returned an unexpected result type")
    formula = sp.simplify(
        _range_formula(val_sym, minimum.value, maximum.value, minimum.attained, maximum.attained)
    )
    _, _, _, _, is_interval, interval_count = _range_metadata_from_formula(formula, val_sym)
    return FunctionRangeResult(
        expr,
        formula,
        val_sym,
        vars_,
        minimum.value,
        maximum.value,
        minimum.attained,
        maximum.attained,
        minimum.points,
        maximum.points,
        "optimization_bounds",
        {"constraints": sp.sstr(condition), "metadata_source": "optimization_bounds"},
        is_interval,
        interval_count,
    )


@with_computation_context
def _compare_bounds(left: sp.Expr, right: sp.Expr) -> int:
    """Compatibility wrapper for the lazily split exact bound comparator."""

    from ._range_special_cases import _compare_bounds as impl

    return impl(left, right)


def function_range(
    expression: sp.Expr,
    constraints: FormulaLike | Iterable[FormulaLike] | None = None,
    variables: Sequence[sp.Symbol | str] | None = None,
    *,
    value_symbol: sp.Symbol | str | None = None,
    domain: str = "reals",
    method: str = "qe",
    return_result: bool = False,
    parameters: Sequence[sp.Symbol | str] | None = None,
    return_stratified: bool = False,
    eliminate_quantifiers: bool = False,
) -> sp.Expr | FunctionRangeResult | object:
    """Return a quantifier-free formula describing a real function range.

    The preferred direct backend uses the semialgebraic image formulation
    ``exists variables. constraints and value_symbol == expression``. It first
    applies guarded exact graph-elimination shortcuts for common univariate
    images, then tries complete CAD/QE, and may use optimization bounds for
    ``method='auto'`` or ``method='bounds'``.
    """

    from ._function_graph_image import (
        _try_complete_cad_image,
        _try_exact_graph_image,
        _try_semialgebraic_graph_image,
        _try_solved_graph_image,
    )
    from ._range_special_cases import _try_direct_special_image

    if domain.lower() not in {"real", "reals", "r", "rr"}:
        raise NotImplementedError("function_range supports only the real domain")
    requested_method = method.lower()
    if requested_method not in {"qe", "cad", "auto", "bounds", "optimization"}:
        raise ValueError("method must be 'qe', 'cad', 'auto', 'bounds', or 'optimization'")
    expr = sp.sympify(expression)
    condition = normalize_constraints(constraints)
    vars_ = normalize_problem_variables(variables, sp.Tuple(condition, expr))
    val_sym = resolve_symbol(value_symbol or "t")
    if eliminate_quantifiers and not return_stratified:
        raise ValueError("eliminate_quantifiers=True requires return_stratified=True")
    if return_stratified:
        from ._optimization_parametric import (
            _normalize_parameters_for_problem,
            _stratified_function_range,
        )

        params = _normalize_parameters_for_problem(parameters or (), expr, condition)
        if not params:
            raise ValueError("return_stratified=True requires at least one parameter")
        vars_ = tuple(var for var in vars_ if var not in set(params))
        return _stratified_function_range(
            expr,
            condition,
            vars_,
            params,
            val_sym,
            eliminate_quantifiers=eliminate_quantifiers,
        )

    original_expr = expr
    original_condition = condition
    original_vars = vars_

    # Algebraic equalities can remove dependent variables before exact
    # transcendental algebraization (for example y=x**2, exp(y)).
    backsub = algebraic_back_substitute(expr, condition, vars_)
    reduced_expr = backsub.expression
    reduced_condition = backsub.constraints
    reduced_vars = backsub.variables

    algebraized = exact_algebraize_function_problem(reduced_expr, reduced_condition, reduced_vars)
    if algebraized is not None:
        image = _try_exact_graph_image(
            algebraized.expression,
            algebraized.constraints,
            algebraized.variables,
            val_sym,
        )
        if image is not None:
            result = _range_result_from_formula(
                original_expr,
                image,
                val_sym,
                original_vars,
                f"algebraized_{algebraized.kind}",
                {
                    "constraints": sp.sstr(original_condition),
                    "requested_method": requested_method,
                    "back_substitutions": tuple(
                        (sp.sstr(v), sp.sstr(value)) for v, value in backsub.substitutions
                    ),
                    "algebraization": algebraized.kind,
                    "algebraization_diagnostics": algebraized.diagnostics,
                },
            )
            return result if return_result else result.formula

    # Back-substitution is also an exact simplification for purely
    # semialgebraic image problems, even when no transcendental rule applies.
    if backsub.substitutions:
        image = _try_exact_graph_image(reduced_expr, reduced_condition, reduced_vars, val_sym)
        if image is not None:
            result = _range_result_from_formula(
                original_expr,
                image,
                val_sym,
                original_vars,
                "qe_image_after_algebraic_back_substitution",
                {
                    "constraints": sp.sstr(original_condition),
                    "requested_method": requested_method,
                    "back_substitutions": tuple(
                        (sp.sstr(v), sp.sstr(value)) for v, value in backsub.substitutions
                    ),
                },
            )
            return result if return_result else result.formula

    if not vars_:
        formula = sp.Eq(val_sym, expr)
        result = FunctionRangeResult(
            expr,
            formula,
            val_sym,
            vars_,
            expr,
            expr,
            True,
            True,
            (),
            (),
            "constant",
            {"metadata_source": "constant"},
            True,
            1,
        )
        return result if return_result else result.formula

    if requested_method in {"bounds", "optimization"}:
        result = _range_via_optimization_bounds(expr, condition, vars_, val_sym, domain)
        return result if return_result else result.formula

    direct_special = _try_direct_special_image(expr, condition, vars_, val_sym)
    if direct_special is not None:
        direct_method = "direct_univariate_image"
        try:
            if len(vars_) == 1 and sp.Poly(sp.expand(expr), vars_[0]).degree() >= 2:
                direct_method = "optimization_bounds_direct_univariate_image"
        except _EXPECTED_ERRORS:
            pass
        result = _range_result_from_formula(
            expr,
            direct_special,
            val_sym,
            vars_,
            direct_method,
            {"constraints": sp.sstr(condition), "requested_method": requested_method},
        )
        return result if return_result else result.formula

    semialgebraic_graph = _try_semialgebraic_graph_image(expr, condition, vars_, val_sym)
    if semialgebraic_graph is not None:
        result = _range_result_from_formula(
            expr,
            semialgebraic_graph,
            val_sym,
            vars_,
            "qe_image_semialgebraic_graph",
            {"constraints": sp.sstr(condition), "requested_method": requested_method},
        )
        return result if return_result else result.formula

    direct = _try_solved_graph_image(expr, condition, vars_, val_sym)
    if direct is not None:
        result = _range_result_from_formula(
            expr,
            direct,
            val_sym,
            vars_,
            "qe_image_solved_graph",
            {"constraints": sp.sstr(condition), "requested_method": requested_method},
        )
        return result if return_result else result.formula

    cad_formula = (
        _try_complete_cad_image(expr, condition, vars_, val_sym)
        if requested_method == "cad"
        else None
    )
    if cad_formula is not None:
        result = _range_result_from_formula(
            expr,
            cad_formula,
            val_sym,
            vars_,
            "qe_image_complete_cad",
            {"constraints": sp.sstr(condition), "requested_method": requested_method},
        )
        return result if return_result else result.formula

    if requested_method == "cad":
        raise NotImplementedError(
            "CAD/QE-image range computation failed for this expression/domain"
        )

    result = _range_via_optimization_bounds(expr, condition, vars_, val_sym, domain)
    result = FunctionRangeResult(
        result.expression,
        result.formula,
        result.value_symbol,
        result.variables,
        result.infimum,
        result.supremum,
        result.minimum_attained,
        result.maximum_attained,
        result.minimizers,
        result.maximizers,
        "optimization_bounds_after_qe_image_fallback",
        {**dict(result.diagnostics), "requested_method": requested_method},
        result.is_interval,
        result.interval_count,
    )
    return result if return_result else result.formula


__all__ = ["function_range"]
