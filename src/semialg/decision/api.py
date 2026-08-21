from __future__ import annotations

from collections.abc import Iterable, Mapping, Sequence

import sympy as sp
from sympy.logic.boolalg import Boolean
from sympy.polys.polyerrors import PolynomialError

from ..algebraic.rational_univariate import solve_formula_with_rur
from ..context import with_computation_context
from ..decision_diagnostics import solution_capability_diagnostics
from ..formula import ParsedPrenexFormula, parse_formula
from ..normalization import conjuncts as _conjuncts
from ..normalization import normalize_formula as _normalize_formula
from ..qe import qe_by_complete_cad
from ..solve import reduce_formula
from ._inputs import (
    as_real_symbol as _as_real_symbol,
)
from ._inputs import (
    normalize_decision_variables as _normalize_variables,
)
from ._inputs import (
    prepare_solve_inputs as _prepare_solve_inputs,
)
from ._metadata import (
    collect_solution_metadata as _collect_solution_metadata,
)
from ._metadata import (
    components_formula as _components_formula,
)
from ._metadata import (
    metadata_request_for_output as _metadata_request_for_output,
)
from ._metadata import (
    one_dim_components as _one_dim_components,
)
from ._outputs import (
    add_standard_solver_diagnostics as _add_standard_solver_diagnostics,
)
from ._outputs import (
    select_solution_output as _select_solution_output,
)
from ._witnesses import find_validated_witness as _find_validated_witness
from .sampling_helpers import _collect_structural_samples, _normalize_sample_request
from .solution import (
    EquivalenceResult,
    ImplicationResult,
    IntervalComponent,
    SatisfiabilityResult,
    SemialgebraicSolution,
    TautologyResult,
)

FormulaLike = sp.Expr | Boolean | bool

_RECOVERABLE_ERRORS = (
    ArithmeticError,
    TypeError,
    ValueError,
    NotImplementedError,
    RuntimeError,
    sp.PolynomialError,
)


def _safe_simplify_expr(expr: sp.Expr) -> sp.Expr:
    """Simplify algebraic expressions without sending Boolean formulas to radsimp."""

    if isinstance(expr, Boolean):
        try:
            return sp.simplify_logic(expr, form="dnf")
        except (TypeError, ValueError, NotImplementedError):
            return expr
    return sp.simplify(expr)


def _merge_variables(
    first: Sequence[sp.Symbol],
    second: Sequence[sp.Symbol],
    *formulas: sp.Expr,
) -> tuple[sp.Symbol, ...]:
    """Merge explicit variables and free symbols while preserving order."""

    out: list[sp.Symbol] = []
    seen: set[sp.Symbol] = set()
    for source in (first, second):
        for sym in source:
            if sym not in seen:
                out.append(sym)
                seen.add(sym)
    for formula in formulas:
        for sym in sorted(getattr(formula, "free_symbols", set()), key=lambda item: item.name):
            if sym not in seen:
                out.append(sym)
                seen.add(sym)
    return tuple(out)


def _fast_parameter_conditions(
    expr: sp.Expr,
    variables: tuple[sp.Symbol, ...],
    parameters: tuple[sp.Symbol, ...],
) -> sp.Expr | None:
    """Fast parameter conditions for common one-variable polynomial atoms."""

    if len(variables) != 1 or not isinstance(expr, sp.core.relational.Relational):
        return None
    var = variables[0]
    try:
        poly = sp.Poly(sp.expand(expr.lhs - expr.rhs), var)
    except _RECOVERABLE_ERRORS:
        return None
    degree = poly.degree()
    if isinstance(expr, sp.Equality):
        if degree == 0:
            return sp.Eq(poly.as_expr(), 0)
        if degree == 1:
            return sp.Ne(poly.LC(), 0)
        if degree == 2:
            coeffs = poly.all_coeffs()
            a2, a1, a0 = coeffs
            return _safe_simplify_expr(a1**2 - 4 * a2 * a0 >= 0)
    return None


def _fast_solution_formula(
    expr: sp.Expr, variables: tuple[sp.Symbol, ...]
) -> tuple[sp.Expr, str, bool | None]:
    """Fast conservative formula/satisfiability path for common solves."""

    if expr is sp.true or expr == sp.true:
        return sp.true, "trivial", True
    if expr is sp.false or expr == sp.false:
        return sp.false, "trivial", False
    if len(variables) == 1:
        components = _one_dim_components(expr, variables[0])
        if components is not None:
            reduced = _components_formula(components)
            return reduced, "one_dimensional_components", bool(components)
        try:
            reduced = sp.reduce_inequalities(list(_conjuncts(expr)), variables[0])
            if reduced is not None:
                return (
                    reduced,
                    "sympy_reduce_inequalities",
                    reduced is not sp.false and reduced != sp.false,
                )
        except _RECOVERABLE_ERRORS:
            pass
    if len(variables) == 2:
        try:
            from ..implicit_geometry import decompose_cylindrical_formula_to_vertical_bounds_2d

            cells = tuple(decompose_cylindrical_formula_to_vertical_bounds_2d(expr, variables))
            if cells:
                return expr, "vertical_bounds_2d", True
        except _RECOVERABLE_ERRORS:
            pass
    return expr, "cad", None


def _point_formula(point: Mapping[sp.Symbol, sp.Expr], variables: Sequence[sp.Symbol]) -> sp.Expr:
    """Return an exact conjunction describing one finite solution point."""

    atoms = [sp.Eq(var, sp.sympify(point[var])) for var in variables]
    return sp.And(*atoms) if atoms else sp.true


def _finite_points_formula(
    points: Sequence[Mapping[sp.Symbol, sp.Expr]], variables: Sequence[sp.Symbol]
) -> sp.Expr:
    """Return an exact finite-set formula from point assignments."""

    if not points:
        return sp.false
    formulas = [_point_formula(point, variables) for point in points]
    return sp.Or(*formulas) if len(formulas) > 1 else formulas[0]


def _try_rur_formula(
    formula: sp.Expr,
    variables: Sequence[sp.Symbol],
    *,
    max_solutions: int | None = None,
):
    """Try the exact RUR finite-system backend for a Boolean formula.

    Returns ``None`` when the formula is outside the supported finite equality
    fragment. A returned object with ``status == 'unknown'`` means RUR saw at
    least one unsupported branch, so callers must not use it as an UNSAT proof.
    """

    if not variables:
        return None
    try:
        result = solve_formula_with_rur(
            formula, tuple(variables), real=True, max_solutions=max_solutions
        )
    except _RECOVERABLE_ERRORS:
        return None
    if result is None:
        return None
    if str(result.status).lower() == "unknown":
        return None
    if result.partial and not result.assignments:
        return None
    return result


def _make_quantified_sentence(
    formula: sp.Expr, variables: Sequence[sp.Symbol]
) -> tuple[tuple[str, sp.Symbol], ...]:
    return tuple(("exists", var) for var in variables)


def _truth_from_qe_result(result) -> bool:
    if result.is_sentence:
        return bool(result.truth_value)
    simplified = _safe_simplify_expr(result.formula)
    if simplified is sp.true or simplified == sp.true:
        return True
    if simplified is sp.false or simplified == sp.false:
        return False
    # A non-sentence result means parameters escaped the requested variable set.
    # Treat satisfiability existentially over remaining free symbols.
    remaining = tuple(sorted(simplified.free_symbols, key=lambda item: item.name))
    if not remaining:
        return bool(simplified)
    return bool(
        qe_by_complete_cad(
            remaining, _make_quantified_sentence(simplified, remaining), parse_formula(simplified)
        ).truth_value
    )


@with_computation_context
def is_satisfiable(
    formula: FormulaLike | Iterable[FormulaLike],
    variables: Sequence[sp.Symbol | str] | None = None,
    *,
    domain: str = "reals",
    strategy: str | None = None,
    return_result: bool = False,
) -> bool | SatisfiabilityResult:
    """Return whether a real semialgebraic formula has a satisfying point.

    By default this preserves the historical boolean API. With
    ``return_result=True`` it returns a ``SatisfiabilityResult`` containing the
    normalized formula, variable order, backend method, and a validated witness
    when one is cheaply available.
    """

    if domain.lower() not in {"real", "reals", "r", "rr"}:
        raise NotImplementedError("is_satisfiable currently supports only the real domain")
    expr = _normalize_formula(formula)
    vars_ = _normalize_variables(variables, expr)
    method = "trivial"
    witness: Mapping[sp.Symbol, sp.Expr] | None = None
    if expr is sp.true or expr == sp.true:
        sat = True
        witness = {var: sp.Integer(0) for var in vars_}
    elif expr is sp.false or expr == sp.false:
        sat = False
    else:
        # First try exact finite-system dispatch. This proves SAT/UNSAT
        # for supported zero-dimensional equality branches without constructing
        # a full CAD of the ambient space.
        rur_result = _try_rur_formula(expr, vars_, max_solutions=1)
        if rur_result is not None and not rur_result.partial:
            sat = bool(rur_result.assignments)
            method = "rational_univariate"
            witness = dict(rur_result.assignments[0]) if rur_result.assignments else None
        else:
            # Try a validated witness before paying for full QE. This never proves
            # unsatisfiability, but it gives cheap structured results for common
            # full-dimensional feasible regions.
            witness = _find_validated_witness(expr, vars_, strategy=strategy)
            if witness is not None:
                sat = True
                method = "validated_sample"
            else:
                if len(vars_) == 1:
                    try:
                        reduced = sp.reduce_inequalities(list(_conjuncts(expr)), vars_[0])
                        comps = _one_dim_components(reduced, vars_[0])
                        if comps is not None:
                            sat = bool(comps)
                            method = "sympy_reduce_inequalities"
                            witness = {vars_[0]: comps[0].sample_point()} if comps else None
                        else:
                            raise ValueError("univariate reduction did not yield components")
                    except _RECOVERABLE_ERRORS:
                        result = qe_by_complete_cad(
                            vars_, _make_quantified_sentence(expr, vars_), parse_formula(expr)
                        )
                        sat = _truth_from_qe_result(result)
                        method = getattr(result, "method", "complete_cad_qe")
                else:
                    result = qe_by_complete_cad(
                        vars_, _make_quantified_sentence(expr, vars_), parse_formula(expr)
                    )
                    sat = _truth_from_qe_result(result)
                    method = getattr(result, "method", "complete_cad_qe")
                if sat and witness is None:
                    witness = _find_validated_witness(expr, vars_, strategy=strategy)
    if return_result:
        return SatisfiabilityResult(
            bool(sat),
            expr,
            vars_,
            witness=witness,
            method=method,
            diagnostics={"domain": domain, "strategy": strategy},
        )
    return bool(sat)


def is_tautology(
    formula: FormulaLike | Iterable[FormulaLike],
    variables: Sequence[sp.Symbol | str] | None = None,
    *,
    domain: str = "reals",
    strategy: str | None = None,
    return_result: bool = False,
) -> bool | TautologyResult:
    """Return whether a real semialgebraic formula is true for all variables.

    With ``return_result=True``, a false result includes a validated
    counterexample whenever the sampling layer can provide one.
    """

    expr = _normalize_formula(formula)
    vars_ = _normalize_variables(variables, expr)
    negated = sp.Not(expr)
    sat = is_satisfiable(negated, vars_, domain=domain, strategy=strategy, return_result=True)
    taut = not bool(sat)
    if return_result:
        return TautologyResult(
            taut,
            expr,
            vars_,
            counterexample=sat.witness if not taut else None,
            method=sat.method,
            diagnostics={"satisfiability": sat.diagnostics},
        )
    return taut


def implies(
    assumptions: FormulaLike | Iterable[FormulaLike],
    conclusion: FormulaLike,
    variables: Sequence[sp.Symbol | str] | None = None,
    *,
    domain: str = "reals",
    strategy: str | None = None,
    return_result: bool = False,
) -> bool | ImplicationResult:
    """Return whether ``assumptions`` imply ``conclusion`` over the reals.

    With ``return_result=True``, invalid implications include a validated
    counterexample satisfying the premise and falsifying the conclusion when
    available.
    """

    premise = _normalize_formula(assumptions)
    consequent = _normalize_formula(conclusion)
    universe = _normalize_variables(variables, sp.And(premise, consequent))
    counterexample_formula = sp.And(premise, sp.Not(consequent))
    sat = is_satisfiable(
        counterexample_formula, universe, domain=domain, strategy=strategy, return_result=True
    )
    valid = not bool(sat)
    if return_result:
        return ImplicationResult(
            valid,
            premise,
            consequent,
            universe,
            counterexample=sat.witness if not valid else None,
            method=sat.method,
            diagnostics={
                "counterexample_formula": sp.sstr(counterexample_formula),
                "satisfiability": sat.diagnostics,
            },
        )
    return valid


def equivalent(
    lhs: FormulaLike,
    rhs: FormulaLike,
    variables: Sequence[sp.Symbol | str] | None = None,
    *,
    domain: str = "reals",
    strategy: str | None = None,
    return_result: bool = False,
) -> bool | EquivalenceResult:
    """Return whether two semialgebraic formulas define the same real set.

    With ``return_result=True``, a false result includes a counterexample from
    the symmetric difference and, when it can be determined cheaply, the failed
    implication direction.
    """

    left = _normalize_formula(lhs)
    right = _normalize_formula(rhs)
    universe = _normalize_variables(variables, sp.And(left, right))
    if left == right:
        if return_result:
            return EquivalenceResult(True, left, right, universe, method="syntactic")
        return True
    try:
        if sp.simplify_logic(sp.Xor(left, right)) is sp.false:
            if return_result:
                return EquivalenceResult(True, left, right, universe, method="logic_simplify")
            return True
    except (TypeError, ValueError, NotImplementedError, AttributeError):
        pass
    if len(universe) == 1:
        try:
            same_set = bool(left.as_set() == right.as_set())
            if same_set:
                if return_result:
                    return EquivalenceResult(True, left, right, universe, method="sympy_set")
                return True
        except (TypeError, ValueError, NotImplementedError, AttributeError):
            pass
    left_not_right = sp.And(left, sp.Not(right))
    right_not_left = sp.And(right, sp.Not(left))
    lnr = is_satisfiable(
        left_not_right, universe, domain=domain, strategy=strategy, return_result=True
    )
    rnl = is_satisfiable(
        right_not_left, universe, domain=domain, strategy=strategy, return_result=True
    )
    equiv = not bool(lnr) and not bool(rnl)
    failed_direction = None
    witness = None
    if bool(lnr) and bool(rnl):
        failed_direction = "both"
        witness = lnr.witness or rnl.witness
    elif bool(lnr):
        failed_direction = "lhs_implies_rhs"
        witness = lnr.witness
    elif bool(rnl):
        failed_direction = "rhs_implies_lhs"
        witness = rnl.witness
    if return_result:
        return EquivalenceResult(
            equiv,
            left,
            right,
            universe,
            counterexample=witness,
            failed_direction=failed_direction,
            method="symmetric_difference",
            diagnostics={"lhs_not_rhs": lnr.diagnostics, "rhs_not_lhs": rnl.diagnostics},
        )
    return equiv


def _parameter_solution_data(
    formula: sp.Expr,
    variables: tuple[sp.Symbol, ...],
    parameters: tuple[sp.Symbol, ...],
    *,
    domain: str,
) -> tuple[sp.Expr | None, object | None]:
    """Compute exact parameter solvability data when parameters are present."""

    if not parameters:
        return None, None
    pieces = _conjuncts(formula)
    param_formula = pieces[0] if len(pieces) == 1 else formula
    conditions = _fast_parameter_conditions(param_formula, variables, parameters)
    if conditions is None:
        from ..parameters import solvability_conditions

        conditions = solvability_conditions(param_formula, variables, parameters, domain=domain)
    decomposition = None
    try:
        from ..parameter_stratification import parameterized_cylindrical_decomposition

        decomposition = parameterized_cylindrical_decomposition(
            formula,
            variables,
            parameters,
            domain=domain,
            specialize_fibers=True,
        )
        if conditions is None:
            conditions = decomposition.parameter_condition
    except (TypeError, ValueError, ArithmeticError, NotImplementedError, PolynomialError):
        decomposition = None
    return conditions, decomposition


def _trivial_solution(
    formula: sp.Expr,
    variables: tuple[sp.Symbol, ...],
    parameters: tuple[sp.Symbol, ...],
    *,
    satisfiable: bool,
    sample_count: int,
    parameter_conditions: sp.Expr | None,
    parameter_decomposition: object | None,
    strategy: str | None,
) -> SemialgebraicSolution:
    """Build the structured result for a constant true or false formula."""

    result_formula = sp.true if satisfiable else sp.false
    samples = (
        tuple({var: sp.Integer(0) for var in variables} for _ in range(1 if sample_count else 0))
        if satisfiable
        else ()
    )
    meta = _collect_solution_metadata(result_formula, variables)
    return SemialgebraicSolution(
        result_formula,
        variables,
        samples,
        satisfiable,
        "trivial",
        solution_capability_diagnostics(formula),
        parameters=parameters,
        simplified_constraints=(),
        parameter_conditions=(
            parameter_conditions
            if parameter_conditions is not None
            else (sp.true if satisfiable else sp.false)
        ),
        parameter_decomposition=parameter_decomposition if satisfiable else None,
        dimension=meta["dimension"],
        bounded=meta["bounded"],
        closed=meta["closed"],
        compact=meta["compact"],
        components=meta["components"],
        cells=meta["cells"],
        cylindrical_solution=meta.get("cylindrical_solution"),
        connectivity=meta.get("connectivity"),
    )


@with_computation_context
def solve_semialgebraic(
    constraints: FormulaLike | Iterable[FormulaLike],
    variables: Sequence[sp.Symbol | str] | None = None,
    *,
    parameters: Sequence[sp.Symbol | str] | None = None,
    domain: str = "reals",
    count: int = 1,
    samples: int | str | None = None,
    sample_mode: str | None = None,
    strategy: str | None = None,
    method: str = "auto",
    variable_order: Sequence[sp.Symbol | str] | None = None,
    projection_order: Sequence[sp.Symbol | str] | None = None,
    normalize_domains: bool = True,
    return_formula: bool = False,
    output: str | None = None,
) -> SemialgebraicSolution | sp.Expr | tuple[object, ...] | bool | None:
    """Reduce, sample, and summarize a semialgebraic system over the reals.

    ``output`` may be used as a convenience selector for common views of the
    solution. The default ``None`` preserves structured result-object behavior.
    Component- and cell-aware sampling is available through
    ``samples="per_component"``, ``samples="per_cell"``, or the equivalent
    ``sample_mode`` keyword. Supported selectors are ``"formula"``, ``"reduced_formula"``,
    ``"piecewise"``, ``"samples"``, ``"components"``, ``"cells"``,
    ``"cylindrical"``, and ``"conditions"``. The ``"conditions"`` selector returns
    the parameter-space condition under which the system is solvable.

    The returned ``SemialgebraicSolution`` includes best-effort metadata such as simplified
    constraints, parameter conditions, dimension, boundedness, compactness,
    exact 1D components, and 2D vertical-bound cells when those analyses are
    supported. Unsupported metadata is reported as ``None`` or an empty tuple
    rather than being guessed.
    """

    expr_original, expr, params, vars_, method_key, domain_normalization = _prepare_solve_inputs(
        constraints,
        variables,
        parameters,
        domain=domain,
        method=method,
        variable_order=variable_order,
        projection_order=projection_order,
        normalize_domains=normalize_domains,
    )
    sample_count, resolved_sample_mode = _normalize_sample_request(count, samples, sample_mode)
    if method_key == "interval" and len(vars_) != 1:
        raise NotImplementedError("method='interval' supports exactly one solve variable")

    parameter_conditions, parameter_decomposition = _parameter_solution_data(
        expr,
        vars_,
        params,
        domain=domain,
    )

    if expr is sp.true or expr == sp.true:
        result = _trivial_solution(
            expr,
            vars_,
            params,
            satisfiable=True,
            sample_count=sample_count,
            parameter_conditions=parameter_conditions,
            parameter_decomposition=parameter_decomposition,
            strategy=strategy,
        )
        return result.formula if return_formula else _select_solution_output(result, output)
    if expr is sp.false or expr == sp.false:
        result = _trivial_solution(
            expr,
            vars_,
            params,
            satisfiable=False,
            sample_count=sample_count,
            parameter_conditions=parameter_conditions,
            parameter_decomposition=parameter_decomposition,
            strategy=strategy,
        )
        return result.formula if return_formula else _select_solution_output(result, output)

    condition_keys = {"conditions", "parameter_conditions", "solvability_conditions"}
    if (
        output is not None
        and output.lower().replace("-", "_") in condition_keys
        and not return_formula
    ):
        return parameter_conditions if parameter_conditions is not None else sp.true

    if params:
        satisfiable = parameter_conditions is not sp.false and parameter_conditions != sp.false
        simplified_constraints = _conjuncts(expr)
        result = SemialgebraicSolution(
            expr if satisfiable else sp.false,
            vars_,
            (),
            bool(satisfiable),
            "parameter_conditions",
            _add_standard_solver_diagnostics(
                solution_capability_diagnostics(
                    expr,
                    selected_output=output,
                    selected_sample_mode=resolved_sample_mode,
                    requested_sample_count=sample_count,
                    has_parameter_conditions=parameter_conditions is not None,
                    has_param_decomp=parameter_decomposition is not None,
                ),
                method="parameter_conditions",
                variables=vars_,
                projection_order=projection_order,
                domain_normalization=domain_normalization,
                metadata={},
                parameter_decomposition=parameter_decomposition,
            ),
            parameters=params,
            simplified_constraints=simplified_constraints,
            parameter_conditions=parameter_conditions,
            parameter_decomposition=parameter_decomposition,
            dimension=None,
            bounded=None,
            closed=None,
            compact=None,
            components=(),
            cells=(),
            cylindrical_solution=None,
            connectivity=None,
        )
        return result.formula if return_formula else _select_solution_output(result, output)

    selected_method: str | None = None
    rur_result = None
    if method_key in {"auto", "rur"}:
        rur_result = _try_rur_formula(
            expr, vars_, max_solutions=sample_count if sample_count else None
        )
        if method_key == "rur" and rur_result is None:
            raise NotImplementedError(
                "method='rur' supports finite zero-dimensional equality branches only"
            )
    if rur_result is not None and not rur_result.partial:
        assignments = tuple(dict(point) for point in rur_result.assignments)
        reduced = _finite_points_formula(assignments, vars_) if assignments else sp.false
        satisfiable = bool(assignments)
        selected_method = "rational_univariate"
        fast_method = selected_method
        fast_satisfiable = satisfiable
        solved = None
        explicit_cyl = None
    elif method_key in {"cad", "qe", "cylindrical", "sampling"}:
        reduced, fast_method, fast_satisfiable = expr, method_key, None
        solved = None
    else:
        reduced, fast_method, fast_satisfiable = _fast_solution_formula(expr, vars_)
        solved = None
    explicit_cyl = None
    if selected_method is None and fast_satisfiable is None and len(vars_) > 1:
        try:
            from ..cad.cells import extract_explicit_cylindrical_solution

            explicit_cyl = extract_explicit_cylindrical_solution(expr, vars_)
        except _RECOVERABLE_ERRORS:
            explicit_cyl = None
    if selected_method is None and fast_satisfiable is None and explicit_cyl is not None:
        reduced = expr
        satisfiable = True
        selected_method = "explicit_cylindrical_bounds"
    elif selected_method is None and fast_satisfiable is None:
        parsed = ParsedPrenexFormula(vars_, (), parse_formula(expr), expr)
        solved = reduce_formula(parsed, domain=domain, return_result=True, strategy=strategy)
        reduced = _safe_simplify_expr(solved.result)
        try:
            satisfiable = is_satisfiable(reduced, vars_, domain=domain, strategy=strategy)
        except PolynomialError:
            # Reconstructed CAD formulas may contain algebraic boundary functions.
            # A non-false CAD reconstruction already carries selected cells, so use
            # it as the satisfiability witness when polynomial parsing is not available.
            satisfiable = reduced is not sp.false and reduced != sp.false
        selected_method = getattr(solved, "method", "cad")
    elif selected_method is None:
        satisfiable = bool(fast_satisfiable)
        selected_method = fast_method

    selected_method = selected_method or fast_method
    simplified_constraints: tuple[sp.Expr, ...] = _conjuncts(reduced)
    # Record a simplified constraint tuple without making full redundancy removal
    # part of the critical solve path. The dedicated ``simplify_system`` API
    # remains available for heavier semantic cleanup.

    final_formula = sp.false if not satisfiable else reduced
    meta = _collect_solution_metadata(
        final_formula,
        vars_,
        request=_metadata_request_for_output(output, resolved_sample_mode),
    )
    if satisfiable:
        samples_out = _collect_structural_samples(
            final_formula,
            expr,
            vars_,
            meta,
            count=sample_count,
            mode=resolved_sample_mode,
            strategy=strategy,
        )
    else:
        samples_out = ()
    diagnostics = dict(getattr(solved, "metadata", {}) or {}) if solved is not None else {}
    diagnostics.update(solution_capability_diagnostics(expr))
    diagnostics["selected_output"] = output
    diagnostics["selected_sample_mode"] = resolved_sample_mode
    diagnostics["requested_sample_count"] = sample_count
    diagnostics["structural_sample_count"] = len(samples_out)
    diagnostics["simplified_constraint_count"] = len(simplified_constraints)
    diagnostics["has_parameter_conditions"] = parameter_conditions is not None
    diagnostics["has_parameter_decomposition"] = parameter_decomposition is not None
    diagnostics["used_rur"] = selected_method == "rational_univariate"
    if selected_method == "rational_univariate" and rur_result is not None:
        diagnostics["rur_solved_branches"] = rur_result.solved_branches
        diagnostics["rur_skipped_branches"] = rur_result.skipped_branches
        diagnostics["rur_notes"] = tuple(rur_result.notes)
    diagnostics = _add_standard_solver_diagnostics(
        diagnostics,
        method=method_key,
        variables=vars_,
        projection_order=projection_order,
        domain_normalization=domain_normalization,
        metadata=meta,
        parameter_decomposition=parameter_decomposition,
        solved=solved,
    )
    result = SemialgebraicSolution(
        final_formula,
        vars_,
        samples_out,
        satisfiable,
        selected_method,
        diagnostics,
        parameters=params,
        simplified_constraints=simplified_constraints,
        parameter_conditions=parameter_conditions,
        parameter_decomposition=parameter_decomposition,
        dimension=meta["dimension"],
        bounded=meta["bounded"],
        closed=meta["closed"],
        compact=meta["compact"],
        components=meta["components"],
        cells=meta["cells"],
        cylindrical_solution=meta.get("cylindrical_solution"),
        connectivity=meta.get("connectivity"),
    )
    return result.formula if return_formula else _select_solution_output(result, output)


def canonicalize_one_dimensional_formula(
    formula: FormulaLike | Iterable[FormulaLike], variable: sp.Symbol | str
) -> sp.Expr:
    """Return a canonical interval-union formula for supported 1D systems."""

    expr = _normalize_formula(formula)
    var = _as_real_symbol(variable)
    components = _one_dim_components(expr, var)
    if components is None:
        try:
            reduced = sp.reduce_inequalities(list(_conjuncts(expr)), var)
            components = _one_dim_components(reduced, var)
        except _RECOVERABLE_ERRORS:
            components = None
    if components is None:
        return _safe_simplify_expr(expr)
    return _components_formula(components)


__all__ = [
    "IntervalComponent",
    "SemialgebraicSolution",
    "EquivalenceResult",
    "ImplicationResult",
    "TautologyResult",
    "SatisfiabilityResult",
    "equivalent",
    "implies",
    "is_satisfiable",
    "is_tautology",
    "canonicalize_one_dimensional_formula",
    "solve_semialgebraic",
]
