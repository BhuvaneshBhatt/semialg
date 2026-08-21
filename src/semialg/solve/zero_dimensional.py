from __future__ import annotations

from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass

import sympy as sp

from ..algebraic.rational_univariate import (
    FilteredRationalUnivariateSolutions,
    RationalUnivariateError,
    solve_and_filter_zero_dimensional_system_with_rur,
)
from ..dimension_validation import zip_equal


@dataclass(frozen=True)
class ZeroDimensionalSolveResult:
    """Exact solutions of a zero-dimensional polynomial system.

    The result is intentionally small and backend-neutral. ``points`` stores
    coordinate tuples in the requested variable order; ``assignments`` exposes
    the same solutions as dictionaries for callers that need witnesses.
    ``backend`` records the engine that produced the points. ``representation``
    is populated for the RUR backend so advanced callers can reuse the exact
    rational-univariate representation without recomputing it.
    """

    variables: tuple[sp.Symbol, ...]
    points: tuple[tuple[sp.Expr, ...], ...]
    backend: str
    status: str = "satisfied"
    representation: object | None = None
    notes: tuple[str, ...] = ()

    @property
    def assignments(self) -> tuple[Mapping[sp.Symbol, sp.Expr], ...]:
        return tuple(
            dict(zip_equal(self.variables, point, context="RUR solution point"))
            for point in self.points
        )

    @property
    def satisfiable(self) -> bool:
        return bool(self.points)


def _as_polynomial_equation(expr: sp.Expr | bool, variables: Sequence[sp.Symbol]) -> sp.Expr | None:
    if expr in (True, sp.true, sp.S.true):
        return None
    if expr in (False, sp.false, sp.S.false):
        return sp.Integer(1)
    if isinstance(expr, sp.Equality):
        residual = sp.expand(expr.lhs - expr.rhs)
    else:
        residual = sp.expand(sp.sympify(expr))
    if residual == 0:
        return None
    sp.Poly(residual, *variables, domain=sp.QQ)
    return residual


def _normalize_equations(
    equations: Iterable[sp.Expr | bool],
    variables: Sequence[sp.Symbol],
) -> tuple[sp.Expr, ...]:
    normalized: list[sp.Expr] = []
    for equation in equations:
        poly = _as_polynomial_equation(equation, variables)
        if poly is not None:
            normalized.append(poly)
    return tuple(normalized)


def is_zero_dimensional(
    equations: Iterable[sp.Expr | bool],
    variables: Sequence[sp.Symbol],
) -> bool:
    """Return whether rational polynomial equations define a finite complex set.

    This function is deliberately conservative: non-polynomial inputs and
    non-rational coefficients return ``False`` rather than raising. Internally
    it uses SymPy's Groebner-basis zero-dimensionality test, which checks for a
    pure-power leading monomial in every variable direction.
    """

    variable_tuple = tuple(variables)
    if not variable_tuple:
        return False
    try:
        polys = _normalize_equations(equations, variable_tuple)
        if not polys:
            return False
        basis = sp.groebner(polys, *variable_tuple, order="grevlex", domain=sp.QQ)
    except (sp.PolynomialError, ValueError, TypeError):
        return False
    return bool(basis.is_zero_dimensional)


def _cad_or_regular_chains_placeholder(
    equations: tuple[sp.Expr, ...],
    inequalities: sp.Expr | bool | Iterable[sp.Expr | bool] | None,
    variables: tuple[sp.Symbol, ...],
    *,
    real: bool,
) -> ZeroDimensionalSolveResult:
    """Fallback hook for future positive-dimensional solvers.

    The semialg package already has CAD/VS entry points for Boolean formulas,
    but they do not currently expose a point-enumerating API compatible with a
    zero-dimensional exact solver. This hook makes the backend dispatch explicit
    while preventing accidental fake point enumeration for curves/surfaces.
    """

    raise RationalUnivariateError(
        "the equality system is not zero-dimensional; use CAD/quantifier "
        "elimination or a future regular-chains backend for positive-dimensional sets"
    )


def solve_zero_dimensional_system(
    equations: Iterable[sp.Expr | bool],
    inequalities: sp.Expr | bool | Iterable[sp.Expr | bool] | None = None,
    vars: Sequence[sp.Symbol] | None = None,
    *,
    variables: Sequence[sp.Symbol] | None = None,
    backend: str = "rur",
    real: bool = True,
    parameter: sp.Symbol | None = None,
    max_separating_attempts: int = 64,
) -> ZeroDimensionalSolveResult:
    """Solve a finite rational polynomial system exactly.

    Parameters
    ----------
    equations:
        Polynomial equations or ``Eq`` objects. Nonzero expressions are read as
        equations equal to zero.
    inequalities:
        Optional relational atoms or a Boolean combination of atoms used to
        filter candidate equality solutions exactly.
    vars, variables:
        Variable order. ``vars`` matches the public API requested for this
        backend; ``variables`` is accepted as a clearer alias.
    backend:
        ``"rur"`` uses the rational-univariate backend. ``"auto"`` uses RUR
        when the equations are zero-dimensional and otherwise dispatches to the
        positive-dimensional fallback hook.
    real:
        When true, return only real solutions. Complex algebraic solutions are
        available for unconstrained equality systems by setting ``real=False``.
    max_separating_attempts:
        Number of deterministic linear forms to try when searching for a
        separating element of the quotient algebra.
    """

    variable_tuple = tuple(variables if variables is not None else (vars or ()))
    if not variable_tuple:
        raise RationalUnivariateError("a nonempty variable list is required")
    equation_tuple = _normalize_equations(equations, variable_tuple)
    if not equation_tuple:
        raise RationalUnivariateError("at least one nontrivial equation is required")

    backend_name = backend.replace("-", "_").lower()
    constraints = sp.true if inequalities is None else inequalities

    if backend_name in {"auto", "rur", "rational_univariate"}:
        if is_zero_dimensional(equation_tuple, variable_tuple):
            filtered: FilteredRationalUnivariateSolutions = (
                solve_and_filter_zero_dimensional_system_with_rur(
                    equation_tuple,
                    variable_tuple,
                    constraints,
                    real=real,
                    parameter=parameter,
                    max_separating_attempts=max_separating_attempts,
                )
            )
            return ZeroDimensionalSolveResult(
                variables=variable_tuple,
                points=filtered.points,
                backend="rational_univariate",
                status="satisfied" if filtered.points else "unsat",
                representation=filtered.representation,
                notes=(
                    f"quotient_dimension={filtered.representation.dimension}",
                    f"geometric_solution_count={filtered.representation.solution_count}",
                    f"separating_linear_form={sp.sstr(filtered.representation.separating_linear_form)}",
                ),
            )
        if backend_name == "rur":
            raise RationalUnivariateError("RUR backend requires a zero-dimensional equality system")
        return _cad_or_regular_chains_placeholder(
            equation_tuple, constraints, variable_tuple, real=real
        )

    if backend_name in {"cad", "regular_chains", "regular_chains_or_cad"}:
        return _cad_or_regular_chains_placeholder(
            equation_tuple, constraints, variable_tuple, real=real
        )

    raise RationalUnivariateError(f"unknown zero-dimensional solver backend: {backend}")
