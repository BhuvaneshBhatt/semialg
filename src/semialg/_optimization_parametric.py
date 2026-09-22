from __future__ import annotations

from collections.abc import Mapping, Sequence
from typing import Literal

import sympy as sp
from sympy.polys.polyerrors import PolynomialError

from ._optimization_backends import qe_by_complete_cad
from ._range_special_cases import _graph_formula_for_expression
from ._zero_testing import certified_equal
from .formula import parse_formula
from .internal_symbols import fresh_real_dummy
from .normalization import normalize_parameters
from .optimization_results import (
    OptimizationResult,
    ParametricFunctionRangeResult,
    ParametricOptimizationResult,
)

_EXPECTED_ERRORS = (TypeError, ValueError, ArithmeticError, NotImplementedError, PolynomialError)


def _validated_parameter_sample(
    condition: sp.Expr,
    parameters: tuple[sp.Symbol, ...],
    candidate: Mapping[sp.Symbol, sp.Expr] | None = None,
) -> dict[sp.Symbol, sp.Expr]:
    """Return an exact representative satisfying a parameter guard."""
    from .instances.real_fallbacks import satisfies_formula
    from .sampling import sample_point

    def valid(point):
        if point is None or any(parameter not in point for parameter in parameters):
            return False
        try:
            return satisfies_formula(condition, point, strict=True)
        except _EXPECTED_ERRORS:
            return False

    if valid(candidate):
        return {parameter: sp.simplify(candidate[parameter]) for parameter in parameters}
    try:
        sampled = sample_point(
            condition,
            parameters,
            strategy="complete",
            strict=True,
            exact=True,
            default_sampling_radius=8,
        )
    except _EXPECTED_ERRORS:
        sampled = None
    if valid(sampled):
        return {parameter: sp.simplify(sampled[parameter]) for parameter in parameters}
    raise NotImplementedError("could not certify a representative parameter sample")


def _normalize_parameters_for_problem(
    parameters: Sequence[sp.Symbol | str],
    *expressions: sp.Expr,
) -> tuple[sp.Symbol, ...]:
    return normalize_parameters(parameters, *expressions)


def _parameter_guards(
    condition: sp.Expr,
    variables: tuple[sp.Symbol, ...],
    parameters: tuple[sp.Symbol, ...],
) -> tuple[sp.Expr, tuple[tuple[sp.Expr, Mapping[sp.Symbol, sp.Expr]], ...]]:
    from .cad_algorithms.cells import extract_cylindrical_solution
    from .parameters import solvability_conditions

    parameter_domain = sp.simplify(solvability_conditions(condition, variables, parameters))
    if parameter_domain is sp.false or parameter_domain == sp.false:
        return sp.false, ()
    try:
        solution = extract_cylindrical_solution(parameter_domain, parameters, selected_only=True)
    except (NotImplementedError, ValueError, TypeError, ArithmeticError, PolynomialError):
        solution = None
    if solution is None or not solution.cells:
        sample = _validated_parameter_sample(parameter_domain, parameters)
        return parameter_domain, ((parameter_domain, sample),)
    guards = []
    for cell in solution.cells:
        guard = cell.as_formula(closed=False)
        sample = _validated_parameter_sample(guard, parameters, cell.sample_point())
        guards.append((guard, sample))
    return parameter_domain, tuple(guards)


def _parametric_optimum_relation_from_problem(
    objective: sp.Expr,
    condition: sp.Expr,
    variables: tuple[sp.Symbol, ...],
    parameters: tuple[sp.Symbol, ...],
    value_symbol: sp.Symbol,
    *,
    kind: str,
) -> tuple[sp.Expr, tuple[tuple[str, sp.Symbol], ...]]:
    """Return an exact first-order definition of a parametric infimum/supremum."""
    bound_vars = tuple(fresh_real_dummy(f"semialg_bound_{i}") for i in range(len(variables)))
    tight_vars = tuple(fresh_real_dummy(f"semialg_tight_{i}") for i in range(len(variables)))
    threshold = fresh_real_dummy("semialg_threshold")
    bound_subs = dict(zip(variables, bound_vars, strict=True))
    tight_subs = dict(zip(variables, tight_vars, strict=True))
    bound_condition = condition.xreplace(bound_subs)
    tight_condition = condition.xreplace(tight_subs)
    bound_objective = objective.xreplace(bound_subs)
    tight_objective = objective.xreplace(tight_subs)
    if kind == "min":
        bound_clause = sp.Or(
            sp.Not(bound_condition), bound_objective >= value_symbol, evaluate=False
        )
        tight_clause = sp.Or(
            threshold <= value_symbol,
            sp.And(tight_condition, tight_objective < threshold, evaluate=False),
            evaluate=False,
        )
    else:
        bound_clause = sp.Or(
            sp.Not(bound_condition), bound_objective <= value_symbol, evaluate=False
        )
        tight_clause = sp.Or(
            threshold >= value_symbol,
            sp.And(tight_condition, tight_objective > threshold, evaluate=False),
            evaluate=False,
        )
    quantifiers = (
        *(("forall", var) for var in bound_vars),
        ("forall", threshold),
        *(("exists", var) for var in tight_vars),
    )
    return sp.And(bound_clause, tight_clause, evaluate=False), tuple(quantifiers)


def _parametric_range_definition(
    expression: sp.Expr,
    condition: sp.Expr,
    variables: tuple[sp.Symbol, ...],
    value_symbol: sp.Symbol,
) -> tuple[sp.Expr, tuple[tuple[str, sp.Symbol], ...]]:
    quantified_vars = tuple(fresh_real_dummy(f"semialg_range_{i}") for i in range(len(variables)))
    substitutions = dict(zip(variables, quantified_vars, strict=True))
    specialized_expression = expression.xreplace(substitutions)
    specialized_condition = condition.xreplace(substitutions)
    graph_formula, aux_symbols = _graph_formula_for_expression(specialized_expression, value_symbol)
    formula = sp.And(specialized_condition, graph_formula, evaluate=False)
    quantified = (*quantified_vars, *aux_symbols)
    return formula, tuple(("exists", var) for var in quantified)


def _direct_var_interval(
    constraints: sp.Expr,
    variable: sp.Symbol,
) -> tuple[sp.Expr | None, bool, sp.Expr | None, bool] | None:
    """Extract one direct symbolic interval for ``variable`` conservatively."""

    atoms = constraints.args if isinstance(constraints, sp.And) else (constraints,)
    lower: tuple[sp.Expr, bool] | None = None
    upper: tuple[sp.Expr, bool] | None = None
    fixed: sp.Expr | None = None
    for atom in atoms:
        if variable not in atom.free_symbols:
            continue
        if isinstance(atom, sp.Equality):
            if atom.lhs == variable and variable not in atom.rhs.free_symbols:
                fixed = atom.rhs
                continue
            if atom.rhs == variable and variable not in atom.lhs.free_symbols:
                fixed = atom.lhs
                continue
            return None
        strict = isinstance(atom, (sp.StrictLessThan, sp.StrictGreaterThan))
        if isinstance(atom, (sp.LessThan, sp.StrictLessThan)):
            if atom.lhs == variable and variable not in atom.rhs.free_symbols:
                candidate = (atom.rhs, strict)
                if upper is not None and certified_equal(upper[0], candidate[0]) is not True:
                    return None
                upper = candidate
                continue
            if atom.rhs == variable and variable not in atom.lhs.free_symbols:
                candidate = (atom.lhs, strict)
                if lower is not None and certified_equal(lower[0], candidate[0]) is not True:
                    return None
                lower = candidate
                continue
        if isinstance(atom, (sp.GreaterThan, sp.StrictGreaterThan)):
            if atom.lhs == variable and variable not in atom.rhs.free_symbols:
                candidate = (atom.rhs, strict)
                if lower is not None and certified_equal(lower[0], candidate[0]) is not True:
                    return None
                lower = candidate
                continue
            if atom.rhs == variable and variable not in atom.lhs.free_symbols:
                candidate = (atom.lhs, strict)
                if upper is not None and certified_equal(upper[0], candidate[0]) is not True:
                    return None
                upper = candidate
                continue
        return None
    if fixed is not None:
        return fixed, False, fixed, False
    return (
        None if lower is None else lower[0],
        False if lower is None else lower[1],
        None if upper is None else upper[0],
        False if upper is None else upper[1],
    )


def _affine_in_variable(expr: sp.Expr, variable: sp.Symbol) -> tuple[sp.Expr, sp.Expr] | None:
    try:
        poly = sp.Poly(sp.expand(expr), variable, domain="EX")
    except (sp.PolynomialError, TypeError, ValueError):
        return None
    if poly.degree() > 1:
        return None
    if poly.degree() == 0:
        return sp.Integer(0), poly.as_expr()
    coeff, offset = poly.all_coeffs()
    return sp.simplify(coeff), sp.simplify(offset)


def _constant_sign(expr: sp.Expr, parameters: tuple[sp.Symbol, ...]) -> int | None:
    if expr.free_symbols & set(parameters):
        return None
    simplified = sp.simplify(expr)
    if simplified == 0:
        return 0
    if simplified.is_positive is True:
        return 1
    if simplified.is_negative is True:
        return -1
    return None


def _simple_parametric_range(
    expression: sp.Expr,
    constraints: sp.Expr,
    variables: tuple[sp.Symbol, ...],
    parameters: tuple[sp.Symbol, ...],
    value_symbol: sp.Symbol,
) -> sp.Expr | None:
    """Return a cheap quantifier-free range relation for common 1D fibers."""

    if not (expression.free_symbols & set(variables)):
        return sp.Eq(value_symbol, expression)
    if len(variables) != 1:
        return None
    var = variables[0]
    affine = _affine_in_variable(expression, var)
    interval = _direct_var_interval(constraints, var)
    if affine is None or interval is None:
        return None
    coeff, offset = affine
    sign = _constant_sign(coeff, parameters)
    if sign is None:
        return None
    lower, lower_strict, upper, upper_strict = interval
    if sign == 0:
        return sp.Eq(value_symbol, offset)
    if lower is None or upper is None:
        return None
    low_value = sp.simplify(coeff * lower + offset)
    high_value = sp.simplify(coeff * upper + offset)
    if sign < 0:
        low_value, high_value = high_value, low_value
        lower_strict, upper_strict = upper_strict, lower_strict
    low_atom = value_symbol > low_value if lower_strict else value_symbol >= low_value
    high_atom = value_symbol < high_value if upper_strict else value_symbol <= high_value
    return sp.And(low_atom, high_atom)


def _simple_parametric_optimum(
    objective: sp.Expr,
    constraints: sp.Expr,
    variables: tuple[sp.Symbol, ...],
    parameters: tuple[sp.Symbol, ...],
    value_symbol: sp.Symbol,
    *,
    kind: str,
) -> sp.Expr | None:
    """Return a cheap exact optimum-value relation for common 1D fibers."""

    if not (objective.free_symbols & set(variables)):
        return sp.Eq(value_symbol, objective)
    if len(variables) != 1:
        return None
    var = variables[0]
    affine = _affine_in_variable(objective, var)
    interval = _direct_var_interval(constraints, var)
    if affine is None or interval is None:
        return None
    coeff, offset = affine
    sign = _constant_sign(coeff, parameters)
    if sign is None:
        return None
    lower, _lower_strict, upper, _upper_strict = interval
    if sign == 0:
        return sp.Eq(value_symbol, offset)
    endpoint = lower if (kind == "min") == (sign > 0) else upper
    if endpoint is None:
        return None
    return sp.Eq(value_symbol, sp.simplify(coeff * endpoint + offset))


def _eliminate_parametric_relation(
    formula: sp.Expr,
    quantifiers: tuple[tuple[str, sp.Symbol], ...],
    free_variables: tuple[sp.Symbol, ...],
) -> sp.Expr:
    """Eliminate a parametric first-order relation on explicit opt-in.

    Parametric optimization and range APIs keep their exact
    first-order relation by default because a second CAD/QE can be much more
    expensive than constructing the relation itself.  Callers that request
    elimination get the complete-CAD result or the backend exception; there is
    no silent downgrade to an uneliminated relation.
    """

    quantified_variables = tuple(var for _, var in quantifiers)
    variables = tuple(dict.fromkeys((*free_variables, *quantified_variables)))
    qe = qe_by_complete_cad(
        variables,
        quantifiers,
        parse_formula(formula),
        free_variables=free_variables,
        return_result=True,
    )
    return sp.simplify(qe.formula)


def _stratified_optimization(
    objective: sp.Expr,
    constraints: sp.Expr,
    variables: tuple[sp.Symbol, ...],
    parameters: tuple[sp.Symbol, ...],
    *,
    kind: str,
    domain: str,
    certification: Literal["auto", "complete", "candidate"],
    range_cost_limit: int,
    recursion_limit: int,
    max_boolean_branches: int = 32,
    eliminate_quantifiers: bool = False,
):
    """Build guarded exact optimization relations over parameter strata."""

    from .conditional import ConditionalBranch, conditional_result
    from .optimization import _optimize, _polynomial_problem

    if not _polynomial_problem(objective, constraints, (*parameters, *variables)):
        raise NotImplementedError("parameter-stratified optimization requires polynomial data")
    value_symbol = sp.Symbol("_semialg_optimum_value", real=True)
    optimum_relation, optimum_quantifiers = _parametric_optimum_relation_from_problem(
        objective, constraints, variables, parameters, value_symbol, kind=kind
    )
    parameter_domain, guards = _parameter_guards(constraints, variables, parameters)
    branches = []
    for guard, sample in guards:
        specialized_constraints = sp.simplify(constraints.subs(sample))
        specialized_objective = sp.simplify(objective.subs(sample))
        sample_result = None
        try:
            sample_result = _optimize(
                specialized_objective,
                specialized_constraints,
                variables,
                kind=kind,
                domain=domain,
                return_result=True,
                certification=certification,
                range_cost_limit=range_cost_limit,
                recursion_limit=recursion_limit,
                max_boolean_branches=max_boolean_branches,
            )
        except _EXPECTED_ERRORS:
            pass
        guarded_relation = sp.And(guard, optimum_relation, evaluate=False)
        branch_quantifiers = optimum_quantifiers
        quantifier_free = False
        method = "parametric_first_order_optimum_relation"
        if eliminate_quantifiers:
            direct = _simple_parametric_optimum(
                objective, constraints, variables, parameters, value_symbol, kind=kind
            )
            if direct is not None:
                guarded_relation = sp.And(guard, direct, evaluate=False)
                method = "parametric_direct_optimum_relation"
            else:
                guarded_relation = _eliminate_parametric_relation(
                    guarded_relation, optimum_quantifiers, (*parameters, value_symbol)
                )
                method = "parametric_qe_optimum_relation"
            branch_quantifiers = ()
            quantifier_free = True
        value = ParametricOptimizationResult(
            objective,
            constraints,
            variables,
            parameters,
            value_symbol,
            kind,
            guarded_relation,
            branch_quantifiers,
            sample_result if isinstance(sample_result, OptimizationResult) else None,
            method,
            True,
            quantifier_free,
        )
        branches.append(ConditionalBranch(guard, value, certified=True, sample=sample))
    return conditional_result(
        parameters,
        branches,
        coverage_condition=parameter_domain,
        complete=True,
        disjoint=True,
        certified=True,
        method="parametric_qe_optimization",
        diagnostics={"kind": kind, "branch_count": len(branches)},
        normalize=False,
    )


def _stratified_function_range(
    expression: sp.Expr,
    constraints: sp.Expr,
    variables: tuple[sp.Symbol, ...],
    parameters: tuple[sp.Symbol, ...],
    value_symbol: sp.Symbol,
    *,
    eliminate_quantifiers: bool = False,
):
    """Compute a parameter-stratified exact function-range relation over certified parameter cells."""
    from .conditional import ConditionalBranch, conditional_result

    relation, quantifiers = _parametric_range_definition(
        expression, constraints, variables, value_symbol
    )
    parameter_domain, guards = _parameter_guards(constraints, variables, parameters)
    branches = []
    for guard, sample in guards:
        guarded_formula = sp.And(guard, relation, evaluate=False)
        branch_quantifiers = quantifiers
        quantifier_free = False
        method = "parametric_first_order_range_relation"
        if eliminate_quantifiers:
            direct = _simple_parametric_range(
                expression, constraints, variables, parameters, value_symbol
            )
            if direct is not None:
                guarded_formula = sp.And(guard, direct, evaluate=False)
                method = "parametric_direct_range_relation"
            else:
                guarded_formula = _eliminate_parametric_relation(
                    guarded_formula, quantifiers, (*parameters, value_symbol)
                )
                method = "parametric_qe_range_relation"
            branch_quantifiers = ()
            quantifier_free = True
        value = ParametricFunctionRangeResult(
            expression,
            constraints,
            variables,
            parameters,
            value_symbol,
            guarded_formula,
            branch_quantifiers,
            method,
            True,
            quantifier_free,
        )
        branches.append(ConditionalBranch(guard, value, certified=True, sample=sample))
    return conditional_result(
        parameters,
        branches,
        coverage_condition=parameter_domain,
        complete=True,
        disjoint=True,
        certified=True,
        method="parametric_qe_range",
        diagnostics={"branch_count": len(branches)},
        normalize=False,
    )
