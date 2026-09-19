"""Convexity engines used by :mod:`semialg.function_analysis`."""

from __future__ import annotations

from collections.abc import Mapping

import sympy as sp

from .conditional import ParameterStratifiedResult, conditional_result
from .context import SemialgebraicContext, computation_context
from .convexity import ConvexityCertificate, convexity_certificate
from .decision import is_satisfiable
from .domain_solve import function_domain, normalize_domain_sensitive_constraints
from .function_analysis import (
    FunctionConvexityResult,
    _algebraic_domain_formula,
    _classify,
    _copy_expression,
    _copy_variables,
    _domain_graph,
    _function_graph,
    _FunctionAnalysisContext,
    _truth_condition,
)
from .function_graph import UnsupportedFunctionGraph, semialgebraic_function_graph
from .internal_symbols import fresh_real_dummy
from .matrix_analysis import MatrixDefinitenessResult, matrix_definiteness
from .normalization import normalize_formula
from .parameters import solvability_conditions
from .presolve import PresolveResult, presolve_semialgebraic


def _jensen_violation_formula(
    expression: sp.Expr,
    domain: sp.Expr,
    variables: tuple[sp.Symbol, ...],
    *,
    sense: str,
) -> tuple[sp.Expr, tuple[sp.Symbol, ...], Mapping[str, object]]:
    """Return an existential counterexample formula for Jensen convexity/concavity."""

    left_variables = _copy_variables(variables, "semialg_convex_left")
    right_variables = _copy_variables(variables, "semialg_convex_right")
    interpolation = fresh_real_dummy("semialg_convex_t")
    interpolation_point = tuple(
        sp.expand(interpolation * left + (1 - interpolation) * right)
        for left, right in zip(left_variables, right_variables, strict=True)
    )

    left_domain = _copy_expression(domain, variables, left_variables)
    right_domain = _copy_expression(domain, variables, right_variables)
    middle_domain = domain.xreplace(dict(zip(variables, interpolation_point, strict=True)))
    left_domain_formula, left_domain_auxiliaries = _domain_graph(left_domain)
    right_domain_formula, right_domain_auxiliaries = _domain_graph(right_domain)
    middle_domain_formula, middle_domain_aux = _domain_graph(middle_domain)

    left_expression = _copy_expression(expression, variables, left_variables)
    right_expression = _copy_expression(expression, variables, right_variables)
    middle_expression = expression.xreplace(dict(zip(variables, interpolation_point, strict=True)))
    left_graph, left_value, left_auxiliaries = _function_graph(
        left_expression, "semialg_convex_left_value"
    )
    right_graph, right_value, right_auxiliaries = _function_graph(
        right_expression, "semialg_convex_right_value"
    )
    middle_graph, middle_value, middle_auxiliaries = _function_graph(
        middle_expression, "semialg_convex_middle_value"
    )

    chord_value = sp.expand(interpolation * left_value + (1 - interpolation) * right_value)
    violation = middle_value > chord_value if sense == "convex" else middle_value < chord_value
    formula = sp.And(
        interpolation > 0,
        interpolation < 1,
        left_domain_formula,
        right_domain_formula,
        middle_domain_formula,
        left_graph,
        right_graph,
        middle_graph,
        violation,
    )
    quantified_variables = tuple(
        dict.fromkeys(
            (
                *left_variables,
                *right_variables,
                interpolation,
                left_value,
                right_value,
                middle_value,
                *left_domain_auxiliaries,
                *right_domain_auxiliaries,
                *middle_domain_aux,
                *left_auxiliaries,
                *right_auxiliaries,
                *middle_auxiliaries,
            )
        )
    )
    metadata = {
        "left_variables": left_variables,
        "right_variables": right_variables,
        "interpolation": interpolation,
    }
    return formula, quantified_variables, metadata


def _epigraph_property(
    expression: sp.Expr,
    domain: sp.Expr,
    variables: tuple[sp.Symbol, ...],
    *,
    sense: str,
) -> tuple[bool | None, ConvexityCertificate | None]:
    """Decide convexity through an exact epi/hypograph when normalization suffices."""

    value = fresh_real_dummy("semialg_function_value")
    try:
        semialgebraic_function_graph(expression, value)
    except UnsupportedFunctionGraph:
        return None, None
    relation = value >= expression if sense == "convex" else value <= expression
    formula = sp.And(domain, relation)
    try:
        normalized = normalize_domain_sensitive_constraints(formula, (*variables, value)).formula
        certificate = convexity_certificate(normalized, (*variables, value))
    except (
        TypeError,
        ValueError,
        NotImplementedError,
        sp.PolynomialError,
        UnsupportedFunctionGraph,
    ):
        return None, None
    return certificate.outcome, certificate


def _definition_property(
    expression: sp.Expr,
    domain: sp.Expr,
    variables: tuple[sp.Symbol, ...],
    *,
    sense: str,
) -> tuple[bool | None, Mapping[str, object] | None]:
    try:
        violation, quantified_variables, metadata = _jensen_violation_formula(
            expression, domain, variables, sense=sense
        )
    except UnsupportedFunctionGraph:
        return None, None
    satisfiability = is_satisfiable(violation, quantified_variables, return_result=True)
    if not bool(satisfiability):
        return True, None
    witness = satisfiability.witness or {}
    left_variables = metadata["left_variables"]
    right_variables = metadata["right_variables"]
    interpolation = metadata["interpolation"]
    counterexample = {
        "left": {
            original: witness.get(copied)
            for original, copied in zip(variables, left_variables, strict=True)
        },
        "right": {
            original: witness.get(copied)
            for original, copied in zip(variables, right_variables, strict=True)
        },
        "t": witness.get(interpolation),
    }
    return False, counterexample


def _definition_parameter_condition(
    expression: sp.Expr,
    domain: sp.Expr,
    variables: tuple[sp.Symbol, ...],
    parameters: tuple[sp.Symbol, ...],
    *,
    sense: str,
) -> sp.Expr | None:
    try:
        violation, quantified_variables, _metadata = _jensen_violation_formula(
            expression, domain, variables, sense=sense
        )
    except UnsupportedFunctionGraph:
        return None
    violation_condition = solvability_conditions(violation, quantified_variables, parameters)
    return sp.simplify_logic(sp.Not(violation_condition), force=True)


def _domain_convex_parameter_condition(
    domain: sp.Expr,
    variables: tuple[sp.Symbol, ...],
    parameters: tuple[sp.Symbol, ...],
) -> sp.Expr | None:
    left_variables = _copy_variables(variables, "semialg_domain_left")
    right_variables = _copy_variables(variables, "semialg_domain_right")
    interpolation = fresh_real_dummy("semialg_domain_t")
    interpolation_point = tuple(
        sp.expand(interpolation * left + (1 - interpolation) * right)
        for left, right in zip(left_variables, right_variables, strict=True)
    )
    try:
        left_formula, left_aux = _domain_graph(_copy_expression(domain, variables, left_variables))
        right_formula, right_aux = _domain_graph(
            _copy_expression(domain, variables, right_variables)
        )
        middle_formula, middle_aux = _domain_graph(
            domain.xreplace(dict(zip(variables, interpolation_point, strict=True)))
        )
    except UnsupportedFunctionGraph:
        return None
    violation = sp.And(
        interpolation > 0,
        interpolation < 1,
        left_formula,
        right_formula,
        sp.Not(middle_formula),
    )
    quantified_variables = tuple(
        dict.fromkeys(
            (
                *left_variables,
                *right_variables,
                interpolation,
                *left_aux,
                *right_aux,
                *middle_aux,
            )
        )
    )
    violation_condition = solvability_conditions(violation, quantified_variables, parameters)
    return sp.simplify_logic(sp.Not(violation_condition), force=True)


def _presolve_affine_equalities(
    expression: sp.Expr,
    domain: sp.Expr,
    variables: tuple[sp.Symbol, ...],
) -> tuple[sp.Expr, sp.Expr, tuple[sp.Symbol, ...], PresolveResult | None]:
    """Reduce explicit affine equalities before Hessian testing when safe."""

    try:
        presolved = presolve_semialgebraic(domain, variables, eliminate=variables)
    except (TypeError, ValueError, NotImplementedError, sp.PolynomialError):
        return expression, domain, variables, None
    if not presolved.substitutions:
        return expression, domain, variables, presolved
    reduced_expression = expression
    for variable, replacement in presolved.substitutions:
        reduced_expression = sp.simplify(reduced_expression.subs(variable, replacement))
    return reduced_expression, presolved.formula, presolved.variables, presolved


def _hessian_analysis(
    expression: sp.Expr,
    domain: sp.Expr,
    variables: tuple[sp.Symbol, ...],
    context: SemialgebraicContext,
) -> tuple[
    bool | None,
    bool | None,
    MatrixDefinitenessResult | None,
    MatrixDefinitenessResult | None,
    Mapping[str, object],
]:
    """Use exact Hessian signs for polynomial or graph-supported smooth functions."""

    if not variables:
        return True, True, None, None, {"affine_by_reduction": True, "polynomial": True}
    try:
        polynomial = True
        sp.Poly(expression, *variables, domain="EX")
    except (sp.PolynomialError, TypeError, ValueError):
        polynomial = False

    try:
        hessian = sp.ImmutableMatrix(sp.hessian(expression, variables))
    except (TypeError, ValueError, NotImplementedError):
        return None, None, None, None, {"polynomial": polynomial}
    unsupported_heads = (sp.DiracDelta, sp.Derivative)
    if any(entry.has(*unsupported_heads) for entry in hessian):
        return None, None, None, None, {"hessian": hessian, "polynomial": polynomial}

    if not polynomial:
        # Require every nonzero Hessian entry to have an exact semialgebraic
        # graph before delegating its sign to matrix_analysis.
        try:
            for entry in hessian:
                if entry != 0:
                    semialgebraic_function_graph(entry, fresh_real_dummy("semialg_hessian_entry"))
        except UnsupportedFunctionGraph:
            return None, None, None, None, {"hessian": hessian, "polynomial": False}

    hessian_domain = domain
    relative_result = None
    try:
        relative_result = context.strict_feasible(relative=True, return_result=True)
        if relative_result.feasible:
            hessian_domain = relative_result.strict_formula
    except (TypeError, ValueError, NotImplementedError, sp.PolynomialError):
        pass
    hessian_context = (
        context
        if hessian_domain == domain
        else SemialgebraicContext(hessian_domain, variables, computation=context.computation)
    )

    entry_signs = tuple(hessian_context.function_sign(entry) for entry in hessian if entry != 0)
    if not entry_signs or all(sign == "zero" for sign in entry_signs):
        return (
            True,
            True,
            None,
            None,
            {
                "hessian": hessian,
                "hessian_domain": hessian_domain,
                "hessian_entry_signs": entry_signs,
                "relative_strict_feasibility": relative_result,
                "polynomial": polynomial,
            },
        )

    try:
        convexity_result = hessian_context.matrix_definiteness(
            hessian, requested="positive_semidefinite", return_result=True
        )
        concavity_result = None
        concavity_outcome = None
        # On a full-dimensional domain, a nonzero PSD Hessian cannot also be
        # NSD.  Avoid a second matrix/CAD query in the overwhelmingly common
        # convex case.  Zero Hessians were already handled above.
        full_dimensional = hessian_domain in (sp.true, True)
        if not full_dimensional:
            try:
                full_dimensional = bool(hessian_context.strict_feasible(relative=False))
            except (TypeError, ValueError, NotImplementedError, sp.PolynomialError):
                full_dimensional = False
        if convexity_result.outcome is True and full_dimensional:
            concavity_outcome = False
        else:
            concavity_result = hessian_context.matrix_definiteness(
                hessian, requested="negative_semidefinite", return_result=True
            )
            concavity_outcome = concavity_result.outcome
    except (TypeError, ValueError, NotImplementedError, sp.PolynomialError):
        return (
            None,
            None,
            None,
            None,
            {
                "hessian": hessian,
                "hessian_domain": hessian_domain,
                "hessian_entry_signs": entry_signs,
                "relative_strict_feasibility": relative_result,
                "polynomial": polynomial,
            },
        )
    return (
        convexity_result.outcome,
        concavity_outcome,
        convexity_result,
        concavity_result,
        {
            "hessian": hessian,
            "hessian_domain": hessian_domain,
            "hessian_entry_signs": entry_signs,
            "relative_strict_feasibility": relative_result,
            "polynomial": polynomial,
        },
    )


def _univariate_curvature_sign_analysis(
    analysis: _FunctionAnalysisContext,
    expression: sp.Expr,
    variable: sp.Symbol,
) -> tuple[bool | None, bool | None, Mapping[str, object]] | None:
    """Classify 1D curvature from one exact second-derivative sign query."""
    if expression.has(sp.Abs, sp.Piecewise, sp.sign):
        return None
    try:
        second = analysis.derivative(variable, 2)
    except (TypeError, ValueError, NotImplementedError):
        return None
    if second.has(sp.DiracDelta, sp.Derivative):
        return None
    try:
        direct_semialgebraic = bool(second.is_rational_function(variable))
    except (TypeError, ValueError):
        direct_semialgebraic = False
    if not direct_semialgebraic:
        try:
            semialgebraic_function_graph(
                second, fresh_real_dummy("semialg_second_derivative_value")
            )
        except UnsupportedFunctionGraph:
            return None
    # Hessian necessity needs a domain with nonempty ambient interior.  On a
    # singleton/lower-dimensional slice every function is affine relative to
    # that set, so defer to the relative-domain machinery.
    try:
        if not bool(analysis.semialgebraic_context.strict_feasible(relative=False)):
            return None
    except (TypeError, ValueError, NotImplementedError, sp.PolynomialError):
        if analysis.domain not in (sp.true, True):
            return None
    try:
        sign = analysis.sign(second)
    except (TypeError, ValueError, NotImplementedError, sp.PolynomialError):
        return None
    mapping = {
        "positive": (True, False),
        "nonnegative": (True, False),
        "negative": (False, True),
        "nonpositive": (False, True),
        "zero": (True, True),
        "mixed": (False, False),
    }
    if sign not in mapping:
        return None
    convex, concave = mapping[sign]
    return (
        convex,
        concave,
        {
            "hessian": sp.ImmutableMatrix([[second]]),
            "second_derivative": second,
            "second_derivative_sign": sign,
            "univariate_curvature_fast_path": True,
            "polynomial": bool(expression.is_polynomial(variable)),
        },
    )


def _function_convexity_without_parameters(
    expression: sp.Expr,
    variables: tuple[sp.Symbol, ...],
    explicit_domain: sp.Expr,
    *,
    analysis: _FunctionAnalysisContext | None = None,
) -> FunctionConvexityResult:
    """Classify convexity and concavity using domain convexity, curvature tests, and exact counterexample searches."""
    analysis = analysis or _FunctionAnalysisContext(expression, variables, explicit_domain)
    natural_domain = analysis.natural_domain
    effective_domain = analysis.domain
    context = analysis.semialgebraic_context

    with computation_context(context.computation):
        domain_certificate = analysis.domain_convexity()
        if domain_certificate.outcome is False:
            return FunctionConvexityResult(
                "nonconvex_domain",
                expression,
                variables,
                effective_domain,
                natural_domain,
                False,
                False,
                False,
                "domain_convexity",
                domain_certificate=domain_certificate,
                counterexample=domain_certificate.details or domain_certificate.witness,
            )

        if not context.satisfiable():
            return FunctionConvexityResult(
                "affine",
                expression,
                variables,
                effective_domain,
                natural_domain,
                True,
                True,
                True,
                "empty_domain_vacuous",
                domain_certificate=domain_certificate,
            )

        relative_feasibility = None
        try:
            relative_feasibility = context.strict_feasible(relative=True, return_result=True)
        except NotImplementedError:
            pass

        reduced_expression, reduced_domain, reduced_variables, presolved = (
            _presolve_affine_equalities(expression, effective_domain, variables)
        )
        reduced_context = (
            context
            if (
                reduced_expression == expression
                and reduced_domain == effective_domain
                and reduced_variables == variables
            )
            else SemialgebraicContext(
                reduced_domain, reduced_variables, computation=context.computation
            )
        )
        curvature_fast = None
        if (
            len(variables) == 1
            and reduced_expression == expression
            and reduced_domain == effective_domain
            and reduced_variables == variables
        ):
            curvature_fast = _univariate_curvature_sign_analysis(analysis, expression, variables[0])
        if curvature_fast is not None:
            convex, concave, hessian_details = curvature_fast
            convexity_result = concavity_result = None
        else:
            convex, concave, convexity_result, concavity_result, hessian_details = (
                _hessian_analysis(
                    reduced_expression, reduced_domain, reduced_variables, reduced_context
                )
            )

        # Positive Hessian certificates are sufficient on a convex domain.
        # Exact Hessian failures are necessary only for polynomial functions on
        # a full-dimensional domain; lower-dimensional and boundary-singular
        # cases retain the definition-based fallback.
        convex_counterexample = None
        concave_counterexample = None
        method_parts: list[str] = []
        if hessian_details.get("univariate_curvature_fast_path"):
            method_parts.append("univariate_second_derivative_sign")
        elif convex is True or concave is True:
            method_parts.append("hessian_matrix_definiteness")

        ambient_interior = reduced_domain is sp.true or reduced_domain == sp.true
        if not ambient_interior:
            try:
                ambient_interior = bool(reduced_context.strict_feasible(relative=False))
            except (TypeError, ValueError, NotImplementedError, sp.PolynomialError):
                ambient_interior = False
        hessian_is_complete = bool(hessian_details.get("polynomial")) and ambient_interior

        # If one direction is already proved, that is a valid primary
        # classification.  We only spend extra work on the opposite direction
        # when it is already decided by a complete polynomial Hessian test;
        # this avoids turning simple nonsmooth convex functions into two large
        # QE problems merely to prove that they are not affine.
        if convex is not True and not (hessian_is_complete and convex is False):
            convex, epigraph_certificate = _epigraph_property(
                expression, effective_domain, variables, sense="convex"
            )
            if epigraph_certificate is not None:
                method_parts.append("epigraph_set_convexity")
            if convex is None:
                convex, convex_counterexample = _definition_property(
                    expression, effective_domain, variables, sense="convex"
                )
                method_parts.append("jensen_graph_qe")

        if (
            convex is not True
            and concave is not True
            and not (hessian_is_complete and concave is False)
        ):
            concave, hypograph_certificate = _epigraph_property(
                expression, effective_domain, variables, sense="concave"
            )
            if hypograph_certificate is not None:
                method_parts.append("hypograph_set_convexity")
            if concave is None:
                concave, concave_counterexample = _definition_property(
                    expression, effective_domain, variables, sense="concave"
                )
                if "jensen_graph_qe" not in method_parts:
                    method_parts.append("jensen_graph_qe")

    classification = _classify(convex, concave)
    return FunctionConvexityResult(
        classification,
        expression,
        variables,
        effective_domain,
        natural_domain,
        domain_certificate.outcome,
        convex,
        concave,
        "+".join(method_parts) if method_parts else "inconclusive",
        domain_certificate=domain_certificate,
        convexity_certificate=convexity_result,
        concavity_certificate=concavity_result,
        counterexample=convex_counterexample or concave_counterexample,
        details={
            "relative_strict_feasibility": relative_feasibility,
            "presolve": presolved,
            **hessian_details,
        },
    )


def _parameterized_hessian_classification(
    expression: sp.Expr,
    domain: sp.Expr,
    variables: tuple[sp.Symbol, ...],
    parameters: tuple[sp.Symbol, ...],
) -> ParameterStratifiedResult | None:
    """Use exact Hessian stratification when it is necessary and sufficient."""

    try:
        sp.Poly(expression, *variables, domain="EX")
    except (sp.PolynomialError, TypeError, ValueError):
        return None
    if domain.free_symbols & set(parameters):
        return None
    try:
        domain_certificate = convexity_certificate(domain, variables)
    except (TypeError, ValueError, NotImplementedError, sp.PolynomialError):
        return None
    if domain_certificate.outcome is not True:
        return None

    # Necessity of the Hessian criterion is used only when the convex domain
    # has nonempty ambient interior.  The full real space is the common cheap
    # case; affine polyhedra use the strict-feasibility primitive.
    full_dimensional = domain is sp.true or domain == sp.true
    if not full_dimensional:
        try:
            context = SemialgebraicContext(domain, variables)
            full_dimensional = bool(context.strict_feasible(relative=False))
        except (TypeError, ValueError, NotImplementedError, sp.PolynomialError):
            full_dimensional = False
    if not full_dimensional:
        return None

    hessian = sp.ImmutableMatrix(sp.hessian(expression, variables))
    # Common univariate parametric case: a scalar Hessian independent of x is
    # already its own complete parameter stratification.  Avoid four separate
    # matrix-definiteness condition computations.
    if len(variables) == 1 and hessian.shape == (1, 1):
        curvature = sp.simplify(hessian[0, 0])
        if not (curvature.free_symbols & set(variables)):
            branches = (
                (curvature > 0, "strongly_convex"),
                (sp.Eq(curvature, 0), "affine"),
                (curvature < 0, "strongly_concave"),
            )
            return conditional_result(
                parameters,
                branches,
                method="function_convexity_hessian_parameter_stratification",
                diagnostics={"hessian": hessian, "curvature": curvature},
            )
    positive_semidefinite = matrix_definiteness(
        hessian,
        variables,
        domain=domain,
        requested="positive_semidefinite",
        parameters=parameters,
    )
    negative_semidefinite = matrix_definiteness(
        hessian,
        variables,
        domain=domain,
        requested="negative_semidefinite",
        parameters=parameters,
    )
    constant_in_variables = not any(entry.free_symbols & set(variables) for entry in hessian)
    positive_definite = negative_definite = None
    if constant_in_variables:
        positive_definite = matrix_definiteness(
            hessian, variables, domain=domain, requested="positive_definite", parameters=parameters
        )
        negative_definite = matrix_definiteness(
            hessian, variables, domain=domain, requested="negative_definite", parameters=parameters
        )
    if not isinstance(positive_semidefinite, ParameterStratifiedResult) or not isinstance(
        negative_semidefinite, ParameterStratifiedResult
    ):
        return None
    convex_condition = _truth_condition(positive_semidefinite)
    concave_condition = _truth_condition(negative_semidefinite)
    strong_convex_condition = (
        _truth_condition(positive_definite)
        if isinstance(positive_definite, ParameterStratifiedResult)
        else sp.false
    )
    strong_concave_condition = (
        _truth_condition(negative_definite)
        if isinstance(negative_definite, ParameterStratifiedResult)
        else sp.false
    )
    no_strong = sp.And(sp.Not(strong_convex_condition), sp.Not(strong_concave_condition))
    branches = (
        (strong_convex_condition, "strongly_convex"),
        (sp.And(sp.Not(strong_convex_condition), strong_concave_condition), "strongly_concave"),
        (sp.And(no_strong, convex_condition, concave_condition), "affine"),
        (sp.And(no_strong, convex_condition, sp.Not(concave_condition)), "convex"),
        (sp.And(no_strong, concave_condition, sp.Not(convex_condition)), "concave"),
        (sp.And(no_strong, sp.Not(convex_condition), sp.Not(concave_condition)), "neither"),
    )
    return conditional_result(
        parameters,
        branches,
        method="function_convexity_hessian_parameter_stratification",
        diagnostics={
            "domain_method": domain_certificate.method,
            "positive_semidefinite": positive_semidefinite,
            "negative_semidefinite": negative_semidefinite,
            "positive_definite": positive_definite,
            "negative_definite": negative_definite,
        },
    )


def _function_convexity_with_parameters(
    expression: sp.Expr,
    variables: tuple[sp.Symbol, ...],
    parameters: tuple[sp.Symbol, ...],
    explicit_domain: sp.Expr,
) -> ParameterStratifiedResult:
    """Partition parameter space by exact convexity and concavity conditions for a function family."""
    natural_domain = function_domain(expression, (*variables, *parameters))
    effective_domain = _algebraic_domain_formula(
        normalize_formula(sp.And(explicit_domain, natural_domain)),
        (*variables, *parameters),
    )

    fast_result = _parameterized_hessian_classification(
        expression, effective_domain, variables, parameters
    )
    if fast_result is not None:
        return fast_result

    domain_convex_condition = _domain_convex_parameter_condition(
        effective_domain, variables, parameters
    )
    convex_condition = _definition_parameter_condition(
        expression, effective_domain, variables, parameters, sense="convex"
    )
    concave_condition = _definition_parameter_condition(
        expression, effective_domain, variables, parameters, sense="concave"
    )
    if domain_convex_condition is None or convex_condition is None or concave_condition is None:
        return conditional_result(
            parameters,
            ((sp.true, "unknown"),),
            method="function_convexity_unsupported_graph",
        )

    valid_domain = domain_convex_condition
    branches = (
        (sp.Not(valid_domain), "nonconvex_domain"),
        (sp.And(valid_domain, convex_condition, concave_condition), "affine"),
        (
            sp.And(valid_domain, convex_condition, sp.Not(concave_condition)),
            "convex",
        ),
        (
            sp.And(valid_domain, concave_condition, sp.Not(convex_condition)),
            "concave",
        ),
        (
            sp.And(valid_domain, sp.Not(convex_condition), sp.Not(concave_condition)),
            "neither",
        ),
    )
    return conditional_result(
        parameters,
        branches,
        method="function_convexity_jensen_parameter_qe",
        diagnostics={
            "domain_convex_condition": domain_convex_condition,
            "convex_condition": convex_condition,
            "concave_condition": concave_condition,
        },
    )


__all__ = []
