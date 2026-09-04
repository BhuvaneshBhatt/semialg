from __future__ import annotations

from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass

import sympy as sp
from sympy.logic.boolalg import Boolean

from ._errors import EXACT_OPERATION_ERRORS as _RECOVERABLE_ERRORS
from .cad_algorithms.cells import CylindricalSolution, extract_cylindrical_solution
from .conditional import ConditionalBranch, ParameterStratifiedResult, conditional_result
from .normalization import normalize_formula, normalize_variables
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

    expr = normalize_formula(constraints)
    params = normalize_variables(parameters, expr, append_context_symbols=False)
    vars_ = normalize_variables(
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
    except _RECOVERABLE_ERRORS:
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
            except _RECOVERABLE_ERRORS:
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
                except _RECOVERABLE_ERRORS:
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
    "ParameterExceptionalPolynomial",
    "ParameterExceptionalAnalysis",
    "exceptional_parameter_analysis",
]


@dataclass(frozen=True)
class ParameterExceptionalPolynomial:
    polynomial: sp.Expr
    source: str
    variable: sp.Symbol | None = None
    parents: tuple[sp.Expr, ...] = ()

    def branches(self) -> tuple[sp.Expr, sp.Expr, sp.Expr]:
        p = sp.expand(self.polynomial)
        return (p < 0, sp.Eq(p, 0), p > 0)


@dataclass(frozen=True)
class ParameterExceptionalAnalysis:
    variables: tuple[sp.Symbol, ...]
    parameters: tuple[sp.Symbol, ...]
    causes: tuple[ParameterExceptionalPolynomial, ...]

    @property
    def exceptional_condition(self) -> sp.Expr:
        if not self.causes:
            return sp.false
        return sp.Or(*(sp.Eq(c.polynomial, 0) for c in self.causes))

    @property
    def branching_polynomials(self) -> tuple[sp.Expr, ...]:
        return tuple(c.polynomial for c in self.causes)


def exceptional_parameter_analysis(
    formula: FormulaLike | Iterable[FormulaLike],
    variables: Sequence[sp.Symbol | str],
    parameters: Sequence[sp.Symbol | str],
) -> ParameterExceptionalAnalysis:
    """Discover parameter polynomials where algebraic problem type can change.

    Causes include leading coefficients (degree drops), coefficient sign
    boundaries, discriminants (multiplicity/root-count changes), and pairwise
    resultants (root collisions/common factors).  Only parameter-only
    polynomials are retained, so every reported branch is a valid parameter
    stratum boundary.
    """

    expr = normalize_formula(formula)
    params = normalize_variables(parameters, expr, append_context_symbols=False)
    vars_ = normalize_variables(variables, expr, append_context_symbols=False, exclude=params)

    def relational_atoms(node: sp.Expr) -> tuple[sp.Expr, ...]:
        if getattr(node, "is_Relational", False):
            return (node,)
        if isinstance(node, (sp.And, sp.Or)):
            return tuple(atom for arg in node.args for atom in relational_atoms(arg))
        if isinstance(node, sp.Not):
            return relational_atoms(node.args[0])
        return tuple()

    polys = [sp.expand(atom.lhs - atom.rhs) for atom in relational_atoms(expr)]
    causes: list[ParameterExceptionalPolynomial] = []
    seen: set[tuple[str, str]] = set()

    def add(poly_expr: sp.Expr, source: str, variable=None, parents=()):
        candidate = sp.expand(poly_expr)
        if candidate == 0 or candidate.free_symbols - set(params):
            return
        try:
            pp = sp.Poly(candidate, *params) if params else sp.Poly(candidate)
            if pp.total_degree() == 0:
                return
            _, primitive = pp.primitive()
            factors = primitive.factor_list()[1]
            factor_polys = [factor for factor, _multiplicity in factors] or [primitive]
        except (sp.PolynomialError, TypeError, ValueError, NotImplementedError):
            return
        for factor_poly in factor_polys:
            factor_expr = sp.expand(factor_poly.as_expr())
            if factor_expr.could_extract_minus_sign():
                factor_expr = -factor_expr
            key = (source, sp.srepr(factor_expr))
            if key not in seen:
                seen.add(key)
                causes.append(
                    ParameterExceptionalPolynomial(factor_expr, source, variable, tuple(parents))
                )

    for expr_poly in polys:
        for var in vars_:
            try:
                p = sp.Poly(expr_poly, var)
            except (sp.PolynomialError, TypeError, ValueError):
                continue
            if p.degree() <= 0:
                continue
            add(p.LC(), "degree_drop", var, (expr_poly,))
            for coeff in p.all_coeffs():
                add(coeff, "coefficient_sign", var, (expr_poly,))
            if p.degree() >= 2:
                try:
                    add(sp.discriminant(p.as_expr(), var), "discriminant", var, (expr_poly,))
                except (sp.PolynomialError, TypeError, ValueError):
                    pass

    for i, left in enumerate(polys):
        for right in polys[i + 1 :]:
            for var in vars_:
                if var not in left.free_symbols or var not in right.free_symbols:
                    continue
                try:
                    add(sp.resultant(left, right, var), "resultant", var, (left, right))
                except (sp.PolynomialError, TypeError, ValueError):
                    pass
    causes.sort(
        key=lambda c: (
            c.source,
            "" if c.variable is None else c.variable.name,
            sp.srepr(c.polynomial),
        )
    )
    return ParameterExceptionalAnalysis(vars_, params, tuple(causes))
