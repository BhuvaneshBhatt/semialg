from __future__ import annotations

from collections.abc import Iterable, Mapping, Sequence

import sympy as sp

from ...formulas.boolean import RELATION_TYPES, relation_residual
from .construction import compute_rational_univariate_representation
from .representation import FilteredRationalUnivariateSolutions, RationalUnivariateError
from .solve import solve_rur_representation


def _numeric_sign_with_refinement(expr: sp.Expr) -> int | None:
    """Resolve the sign of an exact algebraic expression by guarded refinement."""

    for precision in (50, 80, 120, 180, 260, 360):
        try:
            value = sp.N(expr, precision)
        except (TypeError, ValueError, ArithmeticError):
            continue
        if value.has(sp.I):
            try:
                complex_value = complex(value)
            except (TypeError, ValueError, ArithmeticError):
                return None
            tolerance = 10 ** (-(precision // 4))
            if abs(complex_value.imag) > tolerance:
                return None
            real_value = complex_value.real
        else:
            real_value = sp.re(value)
        tolerance = sp.Float(10, precision) ** (-(precision // 4))
        try:
            if real_value > tolerance:
                return 1
            if real_value < -tolerance:
                return -1
        except TypeError:
            continue
    return None


def sign_of_algebraic_expression(expr: sp.Expr) -> int:
    """Return the sign of an exact real algebraic expression.

    SymPy exact sign reasoning is used first. Conservative numerical
    refinement is used only when the expression is algebraic but not simplified
    enough for an immediate exact decision.
    """

    simplified = sp.simplify(sp.cancel(expr))
    if simplified == 0 or simplified.is_zero:
        return 0
    if simplified.is_positive:
        return 1
    if simplified.is_negative:
        return -1
    sign = sp.sign(simplified)
    if sign == 1:
        return 1
    if sign == -1:
        return -1
    if sign == 0:
        return 0
    numeric = _numeric_sign_with_refinement(simplified)
    if numeric is not None:
        return numeric
    raise RationalUnivariateError(f"could not determine algebraic sign of {expr!s}")


def evaluate_relation_at_point(relation: sp.Expr, assignment: Mapping[sp.Symbol, sp.Expr]) -> bool:
    """Evaluate a polynomial relational atom at an exact algebraic point."""

    if isinstance(relation, bool):
        return relation
    if relation in (sp.true, sp.S.true):
        return True
    if relation in (sp.false, sp.S.false):
        return False
    if not isinstance(relation, RELATION_TYPES):
        substituted = sp.simplify(relation.subs(assignment))
        if substituted in (sp.true, sp.S.true, True):
            return True
        if substituted in (sp.false, sp.S.false, False):
            return False
        raise RationalUnivariateError(f"unsupported constraint atom: {relation!s}")

    residual = relation_residual(relation).subs(assignment)
    residual_sign = sign_of_algebraic_expression(residual)
    if isinstance(relation, sp.Equality):
        return residual_sign == 0
    if isinstance(relation, sp.Unequality):
        return residual_sign != 0
    if isinstance(relation, sp.StrictLessThan):
        return residual_sign < 0
    if isinstance(relation, sp.LessThan):
        return residual_sign <= 0
    if isinstance(relation, sp.StrictGreaterThan):
        return residual_sign > 0
    if isinstance(relation, sp.GreaterThan):
        return residual_sign >= 0
    raise RationalUnivariateError(f"unsupported relation: {relation!s}")


def evaluate_boolean_formula_at_point(
    formula: sp.Expr | bool, assignment: Mapping[sp.Symbol, sp.Expr]
) -> bool:
    """Evaluate Boolean combinations of polynomial relations at a point."""

    if isinstance(formula, bool):
        return formula
    if formula in (sp.true, sp.S.true):
        return True
    if formula in (sp.false, sp.S.false):
        return False
    if isinstance(formula, sp.And):
        return all(evaluate_boolean_formula_at_point(arg, assignment) for arg in formula.args)
    if isinstance(formula, sp.Or):
        return any(evaluate_boolean_formula_at_point(arg, assignment) for arg in formula.args)
    if isinstance(formula, sp.Not):
        return not evaluate_boolean_formula_at_point(formula.args[0], assignment)
    return evaluate_relation_at_point(formula, assignment)


def filter_rur_solutions_by_constraints(
    points: Iterable[Sequence[sp.Expr]],
    variables: Sequence[sp.Symbol],
    constraints: sp.Expr | bool | Iterable[sp.Expr | bool] = sp.true,
) -> tuple[tuple[sp.Expr, ...], ...]:
    """Filter candidate RUR points by exact polynomial constraints.

    ``constraints`` may be a single SymPy Boolean formula or an iterable of
    relational atoms. The function preserves point order and keeps precisely
    the points satisfying the constraints.
    """

    variable_tuple = tuple(variables)
    if isinstance(constraints, (bool, sp.logic.boolalg.Boolean, *RELATION_TYPES)):
        formula = constraints
    else:
        formula = sp.And(*tuple(constraints))
    accepted: list[tuple[sp.Expr, ...]] = []
    for raw_point in points:
        point = tuple(sp.simplify(value) for value in raw_point)
        assignment = dict(zip(variable_tuple, point, strict=True))
        if evaluate_boolean_formula_at_point(formula, assignment):
            accepted.append(point)
    return tuple(accepted)


def solve_rur_semialgebraic_system(
    equalities: Iterable[sp.Expr],
    variables: Sequence[sp.Symbol],
    constraints: sp.Expr | bool | Iterable[sp.Expr | bool] = sp.true,
    *,
    real: bool = True,
    parameter: sp.Symbol | None = None,
    as_assignments: bool = False,
    max_separating_attempts: int = 64,
) -> tuple[tuple[sp.Expr, ...], ...] | tuple[Mapping[sp.Symbol, sp.Expr], ...]:
    """Solve a finite equality system and filter its solutions by constraints.

    The equality subsystem is solved by RUR. Inequalities, disequalities, and
    Boolean combinations of constraints are then evaluated at each algebraic
    candidate point.
    """

    equality_tuple = tuple(equalities)
    representation = compute_rational_univariate_representation(
        equality_tuple, variables, parameter, max_separating_attempts=max_separating_attempts
    )
    candidate_points = solve_rur_representation(representation, real=real)
    filtered = filter_rur_solutions_by_constraints(candidate_points, variables, constraints)
    if as_assignments:
        variable_tuple = tuple(variables)
        return tuple(dict(zip(variable_tuple, point, strict=True)) for point in filtered)
    return filtered


def solve_and_filter_zero_dimensional_system_with_rur(
    equalities: Iterable[sp.Expr],
    variables: Sequence[sp.Symbol],
    constraints: sp.Expr | bool | Iterable[sp.Expr | bool] = sp.true,
    *,
    real: bool = True,
    parameter: sp.Symbol | None = None,
    max_separating_attempts: int = 64,
) -> FilteredRationalUnivariateSolutions:
    """Return a structured RUR solution object."""

    equality_tuple = tuple(equalities)
    representation = compute_rational_univariate_representation(
        equality_tuple, variables, parameter, max_separating_attempts=max_separating_attempts
    )
    candidate_points = solve_rur_representation(representation, real=real)
    filtered = filter_rur_solutions_by_constraints(candidate_points, variables, constraints)
    return FilteredRationalUnivariateSolutions(
        variables=tuple(variables),
        representation=representation,
        points=filtered,
    )
