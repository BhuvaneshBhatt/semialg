from __future__ import annotations

from collections.abc import Iterable, Sequence

import sympy as sp
from sympy.logic.boolalg import Boolean
from sympy.polys.polyerrors import CoercionFailed, PolynomialError

from ..context import with_computation_context
from ..decision_diagnostics import solution_capability_diagnostics
from ..formula import ParsedPrenexFormula, parse_formula
from ..formulas.boolean import is_false_expr, is_true_expr
from ..normalization import conjuncts
from ..solve import reduce_formula
from ..solve.planner import (
    affine_presolve,
    exact_linear_feasibility,
    profile_semialgebraic_system,
    reconstruct_affine_solution_formula,
)
from ._inputs import (
    prepare_solve_inputs as _prepare_solve_inputs,
)
from ._metadata import (
    collect_solution_metadata as _collect_solution_metadata,
)
from ._metadata import (
    metadata_request_for_output as _metadata_request_for_output,
)
from ._outputs import (
    add_standard_solver_diagnostics as _add_standard_solver_diagnostics,
)
from ._outputs import (
    select_solution_output as _select_solution_output,
)
from .sampling_helpers import _collect_structural_samples, _normalize_sample_request
from .solution import (
    SemialgebraicSolution,
)

FormulaLike = sp.Expr | Boolean | bool

_RECOVERABLE_ERRORS = (
    ArithmeticError,
    TypeError,
    ValueError,
    NotImplementedError,
    sp.PolynomialError,
    CoercionFailed,
)


from ._predicates import is_satisfiable  # noqa: E402
from .api import (  # noqa: E402
    _fast_solution_formula,
    _finite_points_formula,
    _parameter_solution_data,
    _planner_step,
    _run_auto_solve_plan,
    _safe_simplify_expr,
    _trivial_solution,
    _try_rur_formula,
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

    _expr_original, expr, params, vars_, method_key, domain_normalization = _prepare_solve_inputs(
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

    if is_true_expr(expr):
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
    if is_false_expr(expr):
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
        simplified_constraints = conjuncts(expr)
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
    solved = None
    planner_diagnostics: dict[str, object] = {}

    if method_key == "auto":
        planned = _run_auto_solve_plan(
            expr,
            vars_,
            domain=domain,
            strategy=strategy,
            sample_count=sample_count,
        )
        reduced = planned.reduced
        satisfiable = planned.satisfiable
        selected_method = planned.method
        solved = planned.solved
        rur_result = planned.rur_result
        planner_diagnostics = planned.diagnostics
        fast_method = selected_method
        fast_satisfiable = satisfiable
    else:
        # Explicit method requests retain their strict semantics and bypass the
        # automatic planner except for the exact method-specific fast path.
        if method_key == "linear":
            presolved = affine_presolve(expr, vars_)
            working = presolved.formula
            working_vars = presolved.variables
            linear_profile = profile_semialgebraic_system(working, working_vars)
            if not linear_profile.linear or not linear_profile.conjunctive:
                raise NotImplementedError(
                    "method='linear' requires a conjunctive affine real system"
                )
            linear = exact_linear_feasibility(working, working_vars)
            if linear is None:
                raise NotImplementedError(
                    "method='linear' could not eliminate this affine system exactly"
                )
            satisfiable, certificate = linear
            reduced = (
                reconstruct_affine_solution_formula(presolved, working) if satisfiable else sp.false
            )
            selected_method = "linear_fourier_motzkin"
            fast_method = selected_method
            fast_satisfiable = bool(satisfiable)
            planner_diagnostics = {
                "planner_steps": (
                    _planner_step(
                        "linear_fourier_motzkin",
                        True,
                        "explicit linear method request",
                        eliminated_condition=certificate,
                    ),
                )
            }
        else:
            if method_key == "rur":
                rur_result = _try_rur_formula(
                    expr, vars_, max_solutions=sample_count if sample_count else None
                )
                if rur_result is None:
                    raise NotImplementedError(
                        "method='rur' supports finite zero-dimensional equality branches only"
                    )
                assignments = tuple(dict(point) for point in rur_result.assignments)
                reduced = _finite_points_formula(assignments, vars_) if assignments else sp.false
                satisfiable = bool(assignments)
                selected_method = "rational_univariate"
                fast_method = selected_method
                fast_satisfiable = satisfiable
            elif method_key == "interval":
                reduced, fast_method, fast_satisfiable = _fast_solution_formula(expr, vars_)
                if fast_satisfiable is None:
                    raise NotImplementedError(
                        "method='interval' could not reduce this univariate system"
                    )
                satisfiable = bool(fast_satisfiable)
                selected_method = fast_method
            else:
                # cad/qe/cylindrical/sampling intentionally force the general
                # reducer rather than accepting a structurally cheaper method.
                parsed = ParsedPrenexFormula(vars_, (), parse_formula(expr), expr)
                solved = reduce_formula(
                    parsed, domain=domain, return_result=True, strategy=strategy
                )
                reduced = _safe_simplify_expr(solved.result)
                try:
                    satisfiable = is_satisfiable(reduced, vars_, domain=domain, strategy=strategy)
                except PolynomialError:
                    satisfiable = reduced is not sp.false and reduced != sp.false
                selected_method = getattr(solved, "method", method_key)
                fast_method = selected_method
                fast_satisfiable = bool(satisfiable)

    selected_method = selected_method or fast_method
    simplified_constraints: tuple[sp.Expr, ...] = conjuncts(reduced)
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
    diagnostics.update(planner_diagnostics)
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
