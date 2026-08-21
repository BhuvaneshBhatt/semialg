from __future__ import annotations

from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass

import sympy as sp
from sympy.logic.boolalg import Boolean

from .cad.cells import CylindricalSolution, extract_cylindrical_solution
from .conditional import ConditionalBranch, ParameterStratifiedResult, conditional_result
from .normalization import normalize_formula as _normalize_formula
from .normalization import normalize_variables as _normalize_variables
from .parameters import solvability_conditions

FormulaLike = sp.Expr | Boolean | bool


@dataclass(frozen=True)
class ParameterStratum:
    """One cylindrical stratum in parameter space with a representative fiber.

    ``parameter_cell`` is a nested CAD cell over the parameter variables.
    ``condition`` is the formula describing that parameter cell. ``sample`` is
    a deterministic representative parameter assignment. ``specialized_formula``
    is the original system after substituting that sample. ``solution`` is a
    cylindrical solution of the sampled fiber when extraction succeeds.
    """

    index: int
    parameters: tuple[sp.Symbol, ...]
    variables: tuple[sp.Symbol, ...]
    parameter_cell: object
    condition: sp.Expr
    sample: Mapping[sp.Symbol, sp.Expr]
    specialized_formula: sp.Expr
    solution: CylindricalSolution | None = None

    @property
    def dimension(self) -> int | None:
        return getattr(self.parameter_cell, "dimension", None)

    @property
    def solution_dimension(self) -> int | None:
        return None if self.solution is None else self.solution.dimension

    def as_pair(self) -> tuple[sp.Expr, CylindricalSolution | None]:
        return (self.condition, self.solution)


@dataclass(frozen=True)
class ParameterizedCylindricalDecomposition:
    """Piecewise cylindrical solution stratified by parameter-space CAD cells."""

    formula: sp.Expr
    variables: tuple[sp.Symbol, ...]
    parameters: tuple[sp.Symbol, ...]
    parameter_condition: sp.Expr
    parameter_space_solution: CylindricalSolution | None
    strata: tuple[ParameterStratum, ...]
    method: str = "parameter_space_cylindrical_decomposition"

    @property
    def empty(self) -> bool:
        return self.parameter_condition is sp.false or self.parameter_condition == sp.false

    @property
    def nonempty(self) -> bool:
        return not self.empty

    @property
    def stratum_count(self) -> int:
        return len(self.strata)

    def conditions(self) -> tuple[sp.Expr, ...]:
        return tuple(stratum.condition for stratum in self.strata)

    def as_stratified_result(self) -> ParameterStratifiedResult:
        """Expose the parameter cells themselves as certified guarded values.

        The representative fiber ``solution`` inside each ``ParameterStratum``
        remains sample data; this method does not claim that sampled fiber is
        symbolically constant throughout the stratum.
        """

        branches = [ConditionalBranch(stratum.condition, stratum) for stratum in self.strata]
        return conditional_result(
            self.parameters,
            branches,
            coverage_condition=self.parameter_condition,
            complete=True,
            disjoint=bool(self.parameter_space_solution is not None),
            certified=bool(self.parameter_space_solution is not None),
            method=f"{self.method}+conditional",
            diagnostics={"stratum_count": len(self.strata)},
            normalize=False,
        )


def parameterized_cylindrical_decomposition(
    constraints: FormulaLike | Iterable[FormulaLike],
    variables: Sequence[sp.Symbol | str],
    parameters: Sequence[sp.Symbol | str],
    *,
    domain: str = "reals",
    specialize_fibers: bool = True,
) -> ParameterizedCylindricalDecomposition:
    """Return a CAD-style parameter stratification for a semialgebraic system.

    This is a conservative parameterized solving path. It computes
    the parameter-space feasibility condition, decomposes that condition into
    cylindrical parameter cells, and optionally attaches a representative
    cylindrical solution for the fiber over each parameter sample.
    """

    expr = _normalize_formula(constraints)
    params = _normalize_variables(parameters, expr, append_context_symbols=False)
    vars_ = _normalize_variables(
        variables,
        expr,
        append_context_symbols=False,
        exclude=params,
    )
    if not params:
        raise ValueError("parameterized_cylindrical_decomposition requires at least one parameter")
    param_condition = solvability_conditions(expr, vars_, params, domain=domain)
    if param_condition is sp.false or param_condition == sp.false:
        return ParameterizedCylindricalDecomposition(expr, vars_, params, sp.false, None, ())

    parameter_solution: CylindricalSolution | None = None
    strata: list[ParameterStratum] = []
    try:
        parameter_solution = extract_cylindrical_solution(
            param_condition, params, selected_only=True
        )
    except (NotImplementedError, ValueError, TypeError, ArithmeticError, sp.PolynomialError):
        parameter_solution = None

    if parameter_solution is None or not parameter_solution.cells:
        # Keep the full parameter condition as one stratum when finer exact
        # decomposition is unavailable.
        sample = {param: sp.Integer(0) for param in params}
        specialized = sp.simplify(expr.subs(sample))
        fiber_solution = None
        if specialize_fibers:
            try:
                fiber_solution = extract_cylindrical_solution(
                    specialized, vars_, selected_only=True
                )
            except (
                NotImplementedError,
                ValueError,
                TypeError,
                ArithmeticError,
                sp.PolynomialError,
            ):
                fiber_solution = None
        strata.append(
            ParameterStratum(
                index=0,
                parameters=params,
                variables=vars_,
                parameter_cell=None,
                condition=param_condition,
                sample=sample,
                specialized_formula=specialized,
                solution=fiber_solution,
            )
        )
    else:
        for i, cell in enumerate(parameter_solution.cells):
            sample = cell.sample_point()
            specialized = sp.simplify(expr.subs(sample))
            fiber_solution = None
            if specialize_fibers:
                try:
                    fiber_solution = extract_cylindrical_solution(
                        specialized, vars_, selected_only=True
                    )
                except Exception:
                    fiber_solution = None
            strata.append(
                ParameterStratum(
                    index=i,
                    parameters=params,
                    variables=vars_,
                    parameter_cell=cell,
                    condition=cell.as_formula(closed=False),
                    sample=sample,
                    specialized_formula=specialized,
                    solution=fiber_solution,
                )
            )

    return ParameterizedCylindricalDecomposition(
        expr, vars_, params, param_condition, parameter_solution, tuple(strata)
    )


__all__ = [
    "ParameterStratum",
    "ParameterizedCylindricalDecomposition",
    "parameterized_cylindrical_decomposition",
]
