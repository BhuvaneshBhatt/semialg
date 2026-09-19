from __future__ import annotations

from collections.abc import Mapping, Sequence

import sympy as sp

from .errors import ExactEvaluationFailure
from .region_integral_results import ReducedRegionIntegral


def _antiderivative_endpoint(
    antiderivative: sp.Expr,
    variable: sp.Symbol,
    endpoint: sp.Expr | None,
    *,
    side: str,
) -> sp.Expr:
    if endpoint is None:
        direction = -sp.oo if side == "left" else sp.oo
        value = sp.limit(antiderivative, variable, direction)
    else:
        value = antiderivative.subs(variable, endpoint)
    if value.has(sp.Integral, sp.Limit):
        raise NotImplementedError("parametric endpoint evaluation remained unevaluated")
    return sp.simplify(value)


def _cad_univariate_param_integral(
    integrand: sp.Expr,
    condition: sp.Expr,
    variable: sp.Symbol,
    parameters: tuple[sp.Symbol, ...],
):
    """Integrate a 1D parametric fiber using a delineable full CAD.

    Building CAD in ``parameters + (variable,)`` stratifies parameter space by
    the projection discriminants/resultants needed to keep every fiber root
    delineable and ordered.  Sector bounds are therefore exact algebraic root
    functions on each returned parameter cell.
    """

    from .algebraic.samples import sample_to_expr
    from .cad_algorithms.decomposition import decomp_collins_complete
    from .conditional import ConditionalBranch, conditional_result
    from .formula import formula_polynomials, parse_formula
    from .qe.complete import evaluate_formula_on_cell
    from .reconstruct.cylindrical import _sector_bound_expr, path_condition

    matrix = parse_formula(condition)
    polynomials = formula_polynomials(matrix)
    if not polynomials:
        raise NotImplementedError("parametric CAD integration needs polynomial boundaries")
    all_variables = (*parameters, variable)
    cad = decomp_collins_complete(polynomials, all_variables)
    param_level = len(parameters)
    fiber_level = param_level + 1
    antiderivative = sp.integrate(integrand, variable)
    if antiderivative.has(sp.Integral):
        raise NotImplementedError("no exact antiderivative is available for the parametric fiber")

    branches: list[ConditionalBranch] = []
    base_cells = cad.cells_by_level.get(param_level, ()) if param_level else (None,)
    for index, base in enumerate(base_cells):
        parent_index = None if base is None else base.index
        children = tuple(
            cell
            for cell in cad.cells_by_level.get(fiber_level, ())
            if cell.parent_index == parent_index
        )
        pieces: list[sp.Expr] = []
        for cell in children:
            if cell.kind != "sector":
                continue
            if not evaluate_formula_on_cell(matrix, cell, all_variables):
                continue
            left = _sector_bound_expr(
                cell,
                variable,
                cad.cells_by_level,
                side="left",
                base_variables=parameters,
            )
            right = _sector_bound_expr(
                cell,
                variable,
                cad.cells_by_level,
                side="right",
                base_variables=parameters,
            )
            lower_value = _antiderivative_endpoint(antiderivative, variable, left, side="left")
            upper_value = _antiderivative_endpoint(antiderivative, variable, right, side="right")
            pieces.append(sp.simplify(upper_value - lower_value))
        value = sp.simplify(sum(pieces, sp.Integer(0)))
        if base is None:
            guard = sp.true
            sample = {}
        else:
            guard = path_condition(base, parameters, cad.cells_by_level, closed=False)
            sample = {
                parameter: sample_to_expr(value_sample)
                for parameter, value_sample in zip(parameters, base.sample, strict=True)
            }
        branches.append(
            ConditionalBranch(
                guard,
                value,
                sample=sample,
                certified=True,
                metadata={
                    "integration_method": "parametric_algebraic_root_cad",
                    "parameter_cell": index,
                    "fiber_sector_count": len(pieces),
                },
            )
        )
    return conditional_result(
        parameters,
        branches,
        coverage_condition=sp.true,
        complete=True,
        disjoint=True,
        certified=True,
        method="parametric_algebraic_root_cad_integration",
        diagnostics={
            "cad_cell_count": len(cad.cells),
            "parameter_cell_count": len(base_cells),
        },
        normalize=True,
    )


def _stratified_region_integral(
    integrand: sp.Expr,
    condition: sp.Expr,
    variables: tuple[sp.Symbol, ...],
    parameters: tuple[sp.Symbol, ...],
    *,
    bounds: Sequence[tuple[sp.Symbol | str, object, object]]
    | Mapping[sp.Symbol | str, tuple[object, object]]
    | None,
    method: str,
    precision: int,
    measure_dimension: object,
):
    """Return an exact piecewise integral over semialgebraic parameter strata.

    This path is distinct from ``integrate_over_parametric_region``: here the
    *region formula* depends on symbolic parameters.  Geometry is stratified in
    parameter space, while the integral itself is reduced with symbolic
    parameter-dependent limits.  The current exact reduction intentionally
    requires a reduction that is valid symbolically (not merely at a sampled
    fiber), which keeps branch values certified.
    """

    from .conditional import ConditionalBranch, conditional_result
    from .parameter_stratification import parameterized_cylindrical_decomposition

    if method == "numeric":
        raise NotImplementedError(
            "numeric parameter-dependent integration requires specializing the parameters first"
        )
    if measure_dimension not in (None, "ambient"):
        raise NotImplementedError(
            "parameter-dependent intrinsic integration is not implemented yet"
        )

    from .region_integrate import _evaluate_reduced_integral, reduce_region_integral

    decomposition = parameterized_cylindrical_decomposition(
        condition, variables, parameters, specialize_fibers=False
    )
    try:
        reduced = reduce_region_integral(
            integrand,
            condition,
            variables,
            bounds=bounds,
            parameters=parameters,
        )
        if not isinstance(reduced, ReducedRegionIntegral):
            raise TypeError("region reduction returned an unexpected result type")
        value, exact, evaluation_method = _evaluate_reduced_integral(
            reduced, method="symbolic", precision=precision
        )
        if not exact:
            raise ExactEvaluationFailure("symbolic parametric integration became inexact")
    except (NotImplementedError, ValueError, TypeError, ArithmeticError, sp.PolynomialError):
        if len(variables) == 1 and not bounds:
            return _cad_univariate_param_integral(integrand, condition, variables[0], parameters)
        raise

    branches: list[ConditionalBranch] = []
    if decomposition.parameter_condition is not sp.false:
        if decomposition.strata:
            for stratum in decomposition.strata:
                branches.append(
                    ConditionalBranch(
                        stratum.condition,
                        sp.simplify(value),
                        sample=stratum.sample,
                        certified=True,
                        metadata={
                            "integration_method": reduced.method,
                            "evaluation_method": evaluation_method,
                            "parameter_stratum": stratum.index,
                        },
                    )
                )
        else:
            branches.append(
                ConditionalBranch(
                    decomposition.parameter_condition,
                    sp.simplify(value),
                    certified=True,
                    metadata={
                        "integration_method": reduced.method,
                        "evaluation_method": evaluation_method,
                    },
                )
            )

    # The integral of the indicator-weighted region is zero where the fiber is
    # empty.  Adding this branch makes the result cover the complete requested
    # parameter space rather than only the feasible projection.
    if decomposition.parameter_condition is not sp.true:
        branches.append(
            ConditionalBranch(
                sp.Not(decomposition.parameter_condition),
                sp.Integer(0),
                certified=True,
                metadata={"empty_fiber": True},
            )
        )

    return conditional_result(
        parameters,
        branches,
        coverage_condition=sp.true,
        complete=True,
        disjoint=bool(decomposition.parameter_space_solution is not None),
        certified=True,
        method="parametric_semialgebraic_region_integration",
        diagnostics={
            "integration_method": reduced.method,
            "evaluation_method": evaluation_method,
            "parameter_condition": decomposition.parameter_condition,
            "stratum_count": len(branches),
        },
        normalize=True,
    )
