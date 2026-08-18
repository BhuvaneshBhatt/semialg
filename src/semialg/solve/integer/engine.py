from __future__ import annotations

from collections.abc import Callable, Sequence
from dataclasses import dataclass

import sympy as sp

from .output_normalization import CanonIntSolveResult, canon_int_result


@dataclass(frozen=True)
class IntegerSolveRequest:
    expr: sp.Expr
    variables: tuple[sp.Symbol, ...]
    max_solutions: int | None = None
    search_bound: int | None = None


@dataclass(frozen=True)
class IntSolverRegistration:
    name: str
    runner: Callable[[IntegerSolveRequest], CanonIntSolveResult | None]
    priority: int = 100


def _adapt_result(name: str, raw_result) -> CanonIntSolveResult:
    if raw_result is None:
        return None
    return canon_int_result(
        getattr(raw_result, "variables", ()),
        formula=getattr(raw_result, "formula", sp.false),
        solutions=getattr(raw_result, "solutions", ()),
        method=getattr(raw_result, "method", name),
        complete=bool(getattr(raw_result, "complete", False)),
        provenance=[name],
        metadata=getattr(raw_result, "metadata", {}),
    )


def _run_sum_of_two_squares(req: IntegerSolveRequest):
    from .diophantine import solve_sum_of_two_squares

    return _adapt_result("sum_of_two_squares", solve_sum_of_two_squares(req.expr, req.variables))


def run_special_fams(req: IntegerSolveRequest):
    from .special_families import solve_int_fams

    return solve_int_fams(req.expr, req.variables)


def _run_thue_family(req: IntegerSolveRequest):
    from .thue import solve_binary_bounded

    return solve_binary_bounded(req.expr, req.variables, search_bound=req.search_bound or 200)


def _run_linear_divisibility(req: IntegerSolveRequest):
    from .diophantine import reduce_int_divis

    return _adapt_result("linear_divisibility", reduce_int_divis(req.expr, req.variables))


def run_rec_lin_divis(req: IntegerSolveRequest):
    from .linear_recursion import rec_reduce_int_lin_sys

    return rec_reduce_int_lin_sys(req.expr, req.variables)


def _run_factorization(req: IntegerSolveRequest):
    from .diophantine import solve_int_sys_via_factor

    return _adapt_result(
        "factorization_branching", solve_int_sys_via_factor(req.expr, req.variables)
    )


def run_factor_recursion(req: IntegerSolveRequest):
    from .factorization import solve_int_recursion

    return solve_int_recursion(req.expr, req.variables)


def _run_groebner(req: IntegerSolveRequest):
    from .diophantine import solve_int_recursion2

    return _adapt_result("groebner_recursion", solve_int_recursion2(req.expr, req.variables))


def _run_groebner_recursion(req: IntegerSolveRequest):
    from .groebner_recursion import rec_reduce_sys

    return rec_reduce_sys(req.expr, req.variables)


def _run_modular_pruning(req: IntegerSolveRequest):
    from .diophantine import solve_int_pruning

    return _adapt_result("modular_pruning", solve_int_pruning(req.expr, req.variables))


def default_int_registry() -> list[IntSolverRegistration]:
    return sorted(
        [
            IntSolverRegistration("sum_of_two_squares", _run_sum_of_two_squares, 10),
            IntSolverRegistration("specialized_families", run_special_fams, 15),
            IntSolverRegistration("thue_family_bounded", _run_thue_family, 20),
            IntSolverRegistration("recursive_linear_divisibility", run_rec_lin_divis, 25),
            IntSolverRegistration("linear_divisibility", _run_linear_divisibility, 30),
            IntSolverRegistration("factorization_recursion", run_factor_recursion, 35),
            IntSolverRegistration("factorization_branching", _run_factorization, 40),
            IntSolverRegistration("groebner_recursion", _run_groebner_recursion, 45),
            IntSolverRegistration("direct_groebner_recursion", _run_groebner, 50),
            IntSolverRegistration("modular_pruning", _run_modular_pruning, 60),
        ],
        key=lambda reg: reg.priority,
    )


def run_int_solver_pipeline(
    expr: sp.Expr,
    variables: Sequence[sp.Symbol],
    *,
    max_solutions: int | None = None,
    search_bound: int | None = None,
) -> CanonIntSolveResult | None:
    req = IntegerSolveRequest(
        expr=expr,
        variables=tuple(variables),
        max_solutions=max_solutions,
        search_bound=search_bound,
    )
    for reg in default_int_registry():
        try:
            result = reg.runner(req)
        except Exception:
            result = None
        if result is not None:
            return result
    return None


__all__ = [
    "IntegerSolveRequest",
    "IntSolverRegistration",
    "default_int_registry",
    "run_int_solver_pipeline",
]
