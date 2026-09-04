"""Monotonicity engines used by :mod:`semialg.function_analysis`."""

from __future__ import annotations

from collections.abc import Mapping

import sympy as sp

from .conditional import ParameterStratifiedResult, conditional_result
from .decision import is_satisfiable
from .domain_solve import function_domain
from .function_analysis import (
    FunctionMonotonicityResult,
    _algebraic_domain_formula,
    _classify_monotonicity,
    _domain_graph,
    _function_graph,
    _FunctionAnalysisContext,
)
from .function_graph import UnsupportedFunctionGraph
from .internal_symbols import fresh_real_dummy
from .normalization import normalize_formula
from .parameters import solvability_conditions
from .regions.operations import region_dimension
from .sampling import sample_point


def _monotonicity_violation_formula(
    expression: sp.Expr,
    domain: sp.Expr,
    variable: sp.Symbol,
    *,
    relation: str,
) -> tuple[sp.Expr, tuple[sp.Symbol, ...], Mapping[str, object]]:
    """Build an existential pairwise violation formula for monotonicity."""

    left_variable = fresh_real_dummy(f"semialg_monotonic_left_{variable.name}")
    right_variable = fresh_real_dummy(f"semialg_monotonic_right_{variable.name}")
    left_domain = domain.xreplace({variable: left_variable})
    right_domain = domain.xreplace({variable: right_variable})
    left_domain_formula, left_domain_aux = _domain_graph(left_domain)
    right_domain_formula, right_domain_aux = _domain_graph(right_domain)
    left_expression = expression.xreplace({variable: left_variable})
    right_expression = expression.xreplace({variable: right_variable})
    left_graph, left_value, left_aux = _function_graph(
        left_expression, "semialg_monotonic_left_value"
    )
    right_graph, right_value, right_aux = _function_graph(
        right_expression, "semialg_monotonic_right_value"
    )
    comparisons = {
        "increasing": left_value > right_value,
        "strictly_increasing": left_value >= right_value,
        "decreasing": left_value < right_value,
        "strictly_decreasing": left_value <= right_value,
        "constant": sp.Ne(left_value, right_value),
    }
    if relation not in comparisons:
        raise ValueError(f"unsupported monotonicity relation: {relation!r}")
    formula = sp.And(
        left_variable < right_variable,
        left_domain_formula,
        right_domain_formula,
        left_graph,
        right_graph,
        comparisons[relation],
    )
    quantified = tuple(
        dict.fromkeys(
            (
                left_variable,
                right_variable,
                left_value,
                right_value,
                *left_domain_aux,
                *right_domain_aux,
                *left_aux,
                *right_aux,
            )
        )
    )
    return (
        formula,
        quantified,
        {
            "left_variable": left_variable,
            "right_variable": right_variable,
            "left_value": left_value,
            "right_value": right_value,
        },
    )


def _monotonicity_property(
    expression: sp.Expr,
    domain: sp.Expr,
    variable: sp.Symbol,
    *,
    relation: str,
) -> tuple[bool | None, Mapping[str, object] | None]:
    try:
        violation, quantified, metadata = _monotonicity_violation_formula(
            expression, domain, variable, relation=relation
        )
    except UnsupportedFunctionGraph:
        return None, None
    result = is_satisfiable(violation, quantified, return_result=True)
    if not bool(result):
        return True, None
    witness = result.witness or {}
    return False, {
        "left": witness.get(metadata["left_variable"]),
        "right": witness.get(metadata["right_variable"]),
        "left_value": witness.get(metadata["left_value"]),
        "right_value": witness.get(metadata["right_value"]),
        "relation": relation,
    }


def _monotonicity_parameter_condition(
    expression: sp.Expr,
    domain: sp.Expr,
    variable: sp.Symbol,
    parameters: tuple[sp.Symbol, ...],
    *,
    relation: str,
) -> sp.Expr | None:
    try:
        violation, quantified, _metadata = _monotonicity_violation_formula(
            expression, domain, variable, relation=relation
        )
    except UnsupportedFunctionGraph:
        return None
    violation_condition = solvability_conditions(violation, quantified, parameters)
    return sp.simplify_logic(sp.Not(violation_condition), force=True)


def _function_monotonicity_with_parameters(
    expression: sp.Expr,
    variable: sp.Symbol,
    parameters: tuple[sp.Symbol, ...],
    explicit_domain: sp.Expr,
) -> ParameterStratifiedResult:
    """Partition parameter space by exact monotonicity conditions for a univariate family."""
    natural_domain = function_domain(expression, (variable, *parameters))
    effective_domain = _algebraic_domain_formula(
        normalize_formula(sp.And(explicit_domain, natural_domain)),
        (variable, *parameters),
    )
    conditions = {
        name: _monotonicity_parameter_condition(
            expression, effective_domain, variable, parameters, relation=name
        )
        for name in (
            "constant",
            "increasing",
            "strictly_increasing",
            "decreasing",
            "strictly_decreasing",
        )
    }
    if any(value is None for value in conditions.values()):
        return conditional_result(
            parameters,
            ((sp.true, "unknown"),),
            method="function_monotonicity_unsupported_graph",
        )
    constant = conditions["constant"]
    increasing = conditions["increasing"]
    strict_inc = conditions["strictly_increasing"]
    decreasing = conditions["decreasing"]
    strict_dec = conditions["strictly_decreasing"]
    branches = (
        (constant, "constant"),
        (sp.And(sp.Not(constant), strict_inc), "strictly_increasing"),
        (sp.And(sp.Not(constant), strict_dec), "strictly_decreasing"),
        (
            sp.And(sp.Not(constant), increasing, sp.Not(strict_inc)),
            "increasing",
        ),
        (
            sp.And(sp.Not(constant), decreasing, sp.Not(strict_dec)),
            "decreasing",
        ),
        (sp.And(sp.Not(increasing), sp.Not(decreasing)), "nonmonotonic"),
    )
    return conditional_result(
        parameters,
        branches,
        method="function_monotonicity_pairwise_parameter_qe",
        diagnostics={"conditions": conditions, "domain": effective_domain},
    )


def _graph_function_sign_univariate_domain(
    expression: sp.Expr,
    domain: sp.Expr,
    variable: sp.Symbol,
) -> tuple[str | None, sp.Expr | None]:
    """Classify sign using an exact function graph when direct function_sign cannot."""

    try:
        domain_formula, domain_aux = _domain_graph(domain)
        graph_formula, value, graph_aux = _function_graph(
            expression, "semialg_monotonic_derivative_value"
        )
    except UnsupportedFunctionGraph:
        return None, None
    auxiliaries = tuple(dict.fromkeys((*domain_aux, value, *graph_aux)))
    base = sp.And(domain_formula, graph_formula)
    has_positive = bool(is_satisfiable(sp.And(base, value > 0), (variable, *auxiliaries)))
    has_negative = bool(is_satisfiable(sp.And(base, value < 0), (variable, *auxiliaries)))
    has_zero = bool(is_satisfiable(sp.And(base, sp.Eq(value, 0)), (variable, *auxiliaries)))
    if has_positive and has_negative:
        sign = "mixed"
    elif has_positive:
        sign = "nonnegative" if has_zero else "positive"
    elif has_negative:
        sign = "nonpositive" if has_zero else "negative"
    else:
        sign = "zero"
    return sign, None


def _function_monotonicity_without_parameters(
    expression: sp.Expr,
    variable: sp.Symbol,
    explicit_domain: sp.Expr,
    *,
    analysis: _FunctionAnalysisContext | None = None,
) -> FunctionMonotonicityResult:
    """Classify monotonicity on a parameter-free univariate domain using exact derivative and pairwise-order tests."""
    analysis = analysis or _FunctionAnalysisContext(expression, (variable,), explicit_domain)
    natural_domain = analysis.natural_domain
    effective_domain = analysis.domain
    if not is_satisfiable(effective_domain, (variable,)):
        return FunctionMonotonicityResult(
            "constant",
            expression,
            variable,
            effective_domain,
            natural_domain,
            True,
            True,
            True,
            True,
            True,
            "empty_domain_vacuous",
            details={"empty_domain": True},
        )

    analysis_expression = expression
    if expression.has(sp.Abs, sp.sign):
        try:
            if not bool(is_satisfiable(sp.And(effective_domain, variable < 0), (variable,))):
                analysis_expression = sp.refine(expression, sp.Q.nonnegative(variable))
            elif not bool(is_satisfiable(sp.And(effective_domain, variable > 0), (variable,))):
                analysis_expression = sp.refine(expression, sp.Q.nonpositive(variable))
        except (TypeError, ValueError, NotImplementedError, sp.PolynomialError):
            pass

    derivative = (
        analysis.derivative(variable)
        if analysis_expression == expression
        else sp.diff(analysis_expression, variable)
    )
    derivative_sign = None
    context = analysis.semialgebraic_context
    sign_domain = effective_domain
    try:
        relative = context.strict_feasible(relative=True, return_result=True)
        if relative.feasible:
            sign_domain = relative.strict_formula
    except (TypeError, ValueError, NotImplementedError, sp.PolynomialError):
        pass
    derivative_zero_formula = None
    try:
        derivative_sign = analysis.sign(derivative, domain=sign_domain)
    except (TypeError, ValueError, NotImplementedError, sp.PolynomialError):
        derivative_sign, derivative_zero_formula = _graph_function_sign_univariate_domain(
            derivative, sign_domain, variable
        )

    try:
        domain_certificate = analysis.domain_convexity()
        domain_convex = domain_certificate.outcome
    except (TypeError, ValueError, NotImplementedError, sp.PolynomialError):
        domain_convex = None

    derivative_fast_path = domain_convex is True and not expression.has(sp.Piecewise)
    if derivative_fast_path and derivative_sign in {
        "positive",
        "negative",
        "nonnegative",
        "nonpositive",
        "zero",
        "mixed",
    }:
        details = {
            "derivative": derivative,
            "derivative_sign": derivative_sign,
            "domain_convex": True,
        }
        if derivative_sign == "zero":
            return FunctionMonotonicityResult(
                "constant",
                expression,
                variable,
                effective_domain,
                natural_domain,
                True,
                True,
                False,
                False,
                True,
                "derivative_sign",
                details=details,
            )
        if derivative_sign == "mixed":
            return FunctionMonotonicityResult(
                "nonmonotonic",
                expression,
                variable,
                effective_domain,
                natural_domain,
                False,
                False,
                False,
                False,
                False,
                "derivative_sign",
                details=details,
            )
        if derivative_sign in {"positive", "negative"}:
            increasing = derivative_sign == "positive"
            return FunctionMonotonicityResult(
                "strictly_increasing" if increasing else "strictly_decreasing",
                expression,
                variable,
                effective_domain,
                natural_domain,
                increasing,
                not increasing,
                increasing,
                not increasing,
                False,
                "derivative_sign",
                details=details,
            )

        # For a nonzero univariate rational derivative, its zero set is finite on
        # every interval component; no general region-dimension computation is
        # needed to certify strictness.  Fall back to dimension only for graph-
        # supported algebraic derivatives.
        zero_dimension = None
        strict = False
        try:
            rational = derivative.is_rational_function(variable)
        except (TypeError, ValueError):
            rational = False
        if rational and sp.cancel(derivative) != 0:
            zero_dimension = 0
            strict = True
        else:
            zero_set = derivative_zero_formula or sp.And(effective_domain, sp.Eq(derivative, 0))
            try:
                zero_dimension = region_dimension(zero_set, (variable,))
                strict = zero_dimension < 1
            except (TypeError, ValueError, NotImplementedError, sp.PolynomialError):
                pass
        details["derivative_zero_set_dimension"] = zero_dimension
        if derivative_sign == "nonnegative":
            return FunctionMonotonicityResult(
                "strictly_increasing" if strict else "increasing",
                expression,
                variable,
                effective_domain,
                natural_domain,
                True,
                False,
                strict,
                False,
                False,
                "derivative_sign+rational_zero_structure"
                if rational
                else "derivative_sign+zero_set_dimension",
                details=details,
            )
        return FunctionMonotonicityResult(
            "strictly_decreasing" if strict else "decreasing",
            expression,
            variable,
            effective_domain,
            natural_domain,
            False,
            True,
            False,
            strict,
            False,
            "derivative_sign+rational_zero_structure"
            if rational
            else "derivative_sign+zero_set_dimension",
            details=details,
        )

    # On disconnected domains, use the already-cheap exact monotonicity
    # partition to reject impossible global directions before any pairwise QE.
    # Exact sample pairs are certificates of failure; only a surviving candidate
    # direction is sent to the quantified fallback.
    partition = None
    try:
        partition = analysis.monotonic_partition(variable)
    except (TypeError, ValueError, NotImplementedError, sp.PolynomialError):
        pass
    if partition is not None and partition.pieces:
        classes = [classification for classification, _ in partition.pieces]
        local_inc = all(c in {"constant", "increasing", "strictly_increasing"} for c in classes)
        local_dec = all(c in {"constant", "decreasing", "strictly_decreasing"} for c in classes)
        cross_inc_ok = local_inc
        cross_dec_ok = local_dec
        cross_counterexample = None
        samples = []
        for _classification, region in partition.pieces:
            try:
                point = sample_point(region, (variable,))[variable]
                samples.append((point, sp.simplify(expression.subs(variable, point))))
            except (TypeError, ValueError, NotImplementedError, sp.PolynomialError, KeyError):
                samples = []
                break
        for i in range(len(samples)):
            for j in range(i + 1, len(samples)):
                x1, y1 = samples[i]
                x2, y2 = samples[j]
                try:
                    order = bool(sp.Lt(x1, x2))
                except TypeError:
                    order = False
                if not order:
                    continue
                if cross_inc_ok and sp.simplify(y1 > y2) is sp.true:
                    cross_inc_ok = False
                    cross_counterexample = {
                        "left": x1,
                        "right": x2,
                        "left_value": y1,
                        "right_value": y2,
                    }
                if cross_dec_ok and sp.simplify(y1 < y2) is sp.true:
                    cross_dec_ok = False
                    cross_counterexample = {
                        "left": x1,
                        "right": x2,
                        "left_value": y1,
                        "right_value": y2,
                    }
        if not cross_inc_ok and not cross_dec_ok:
            return FunctionMonotonicityResult(
                "nonmonotonic",
                expression,
                variable,
                effective_domain,
                natural_domain,
                False,
                False,
                False,
                False,
                False,
                "partition_cross_component_counterexample",
                counterexample=cross_counterexample,
                details={
                    "derivative": derivative,
                    "derivative_sign": derivative_sign,
                    "partition": partition,
                },
            )

    # Exact definition fallback.  Query only properties that remain plausible
    # after derivative/partition analysis instead of unconditionally solving all five.
    plausible = []
    if partition is None:
        plausible = [
            "constant",
            "increasing",
            "strictly_increasing",
            "decreasing",
            "strictly_decreasing",
        ]
    else:
        classes = [classification for classification, _ in partition.pieces]
        if all(c == "constant" for c in classes):
            plausible.append("constant")
        if all(c in {"constant", "increasing", "strictly_increasing"} for c in classes):
            plausible.extend(("strictly_increasing", "increasing"))
        if all(c in {"constant", "decreasing", "strictly_decreasing"} for c in classes):
            plausible.extend(("strictly_decreasing", "decreasing"))
        if not plausible:
            return FunctionMonotonicityResult(
                "nonmonotonic",
                expression,
                variable,
                effective_domain,
                natural_domain,
                False,
                False,
                False,
                False,
                False,
                "monotonicity_partition",
                details={
                    "derivative": derivative,
                    "derivative_sign": derivative_sign,
                    "partition": partition,
                },
            )
    outcomes = {
        name: None
        for name in (
            "constant",
            "increasing",
            "strictly_increasing",
            "decreasing",
            "strictly_decreasing",
        )
    }
    counterexamples = {}
    for relation in dict.fromkeys(plausible):
        outcomes[relation], counterexamples[relation] = _monotonicity_property(
            expression, effective_domain, variable, relation=relation
        )
        if outcomes[relation] is True and relation in {
            "strictly_increasing",
            "strictly_decreasing",
            "constant",
        }:
            break
    # Fill logical implications without more QE.
    if outcomes["strictly_increasing"] is True:
        (
            outcomes["increasing"],
            outcomes["decreasing"],
            outcomes["strictly_decreasing"],
            outcomes["constant"],
        ) = True, False, False, False
    elif outcomes["strictly_decreasing"] is True:
        (
            outcomes["decreasing"],
            outcomes["increasing"],
            outcomes["strictly_increasing"],
            outcomes["constant"],
        ) = True, False, False, False
    elif outcomes["constant"] is True:
        outcomes["increasing"], outcomes["decreasing"] = True, True
    classification = _classify_monotonicity(**outcomes)
    if classification == "unknown" and partition is not None:
        classification = "nonmonotonic"
        outcomes["increasing"] = (
            outcomes["increasing"] if outcomes["increasing"] is not None else False
        )
        outcomes["decreasing"] = (
            outcomes["decreasing"] if outcomes["decreasing"] is not None else False
        )
    return FunctionMonotonicityResult(
        classification,
        expression,
        variable,
        effective_domain,
        natural_domain,
        outcomes["increasing"],
        outcomes["decreasing"],
        outcomes["strictly_increasing"],
        outcomes["strictly_decreasing"],
        outcomes["constant"],
        "partition+selective_pairwise_qe",
        counterexample=next((v for v in counterexamples.values() if v), None),
        details={
            "derivative": derivative,
            "derivative_sign": derivative_sign,
            "partition": partition,
            "queried_properties": tuple(counterexamples),
        },
    )


__all__ = []
