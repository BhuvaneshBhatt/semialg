from __future__ import annotations

from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass
from typing import Literal

import sympy as sp
from sympy.core.sympify import SympifyError
from sympy.logic.boolalg import And as SymAnd
from sympy.logic.boolalg import BooleanFalse, BooleanTrue
from sympy.logic.boolalg import Not as SymNot
from sympy.logic.boolalg import Or as SymOr
from sympy.polys.polyerrors import GeneratorsNeeded, PolynomialError

from .context import with_computation_context
from .exact_arithmetic import compare_exact_reals
from .formula import parse_formula
from .formulas.boolean import bounded_dnf_branches
from .interval_decomposition import (
    finite_real_roots as _finite_real_roots,
)
from .interval_decomposition import (
    one_dimensional_intervals as _one_dimensional_intervals,
)
from .interval_decomposition import (
    relational_polynomials as _relational_polynomials,
)
from .normalization import (
    normalize_formula as _shared_normalize_formula,
)
from .normalization import (
    normalize_problem_variables as _shared_normalize_variables,
)
from .optimization_active_sets import (
    jacobian_rank_equations as _jacobian_rank_deficiency_equations,
)
from .optimization_active_sets import (
    kkt_system as _kkt_system,
)
from .optimization_active_sets import (
    pruned_active_subsets as _pruned_active_subsets,
)
from .optimization_geometry import polynomial_locus_dimension
from .optimization_results import (
    FunctionRangeResult,
    OptimizationCertificationPolicy,
    OptimizationResult,
    ParametricFunctionRangeResult,
    ParametricOptimizationResult,
)
from .relations import split_relation as _relation_parts
from .symbol_resolution import resolve_symbol

_EXPECTED_ERRORS = (
    TypeError,
    ValueError,
    ArithmeticError,
    NotImplementedError,
    SympifyError,
    PolynomialError,
    GeneratorsNeeded,
)

FormulaLike = sp.Expr | sp.logic.boolalg.Boolean | bool


def qe_by_complete_cad(*args, **kwargs):
    """Import the CAD backend only when a range computation needs it."""

    from .qe import qe_by_complete_cad as impl

    return impl(*args, **kwargs)


def satisfies_formula(*args, **kwargs):
    """Import instance checking only for optimization candidate validation."""

    from .instances.real_fallbacks import satisfies_formula as impl

    return impl(*args, **kwargs)


def _normalize_variables(
    variables: Sequence[sp.Symbol | str] | None,
    expr: sp.Expr,
) -> tuple[sp.Symbol, ...]:
    """Normalize optimization variables and append remaining problem symbols."""

    return _shared_normalize_variables(variables, expr)


def _normalize_formula(
    formula: FormulaLike | Iterable[FormulaLike] | None,
) -> sp.Expr:
    """Normalize optional optimization constraints to one SymPy formula."""

    if formula is None:
        return sp.true
    if isinstance(formula, Iterable) and not isinstance(
        formula,
        (sp.Basic, sp.logic.boolalg.Boolean, str),
    ):
        return sp.And(*(_shared_normalize_formula(item) for item in formula))
    return _shared_normalize_formula(formula)


@dataclass(frozen=True)
class _Candidate:
    value: sp.Expr
    point: Mapping[sp.Symbol, sp.Expr] | None
    attained: bool


def _atoms(condition: sp.Expr) -> tuple[sp.Expr, ...]:
    if condition is sp.true or isinstance(condition, BooleanTrue):
        return ()
    if condition is sp.false or isinstance(condition, BooleanFalse):
        return (sp.false,)
    if isinstance(condition, SymAnd):
        out: list[sp.Expr] = []
        for arg in condition.args:
            out.extend(_atoms(arg))
        return tuple(out)
    if isinstance(condition, (SymOr, SymNot)):
        raise NotImplementedError(
            "optimization candidate enumeration currently supports conjunctions"
        )
    if getattr(condition, "is_Relational", False):
        return (condition,)
    raise TypeError(f"unsupported formula expression: {condition!r}")


def _is_feasible(condition: sp.Expr, point: Mapping[sp.Symbol, sp.Expr]) -> bool:
    try:
        return bool(satisfies_formula(condition, point, strict=False))
    except _EXPECTED_ERRORS:
        value = sp.simplify(condition.subs(point))
        if value is sp.true or isinstance(value, BooleanTrue):
            return True
        if value is sp.false or isinstance(value, BooleanFalse):
            return False
        return bool(value)


def _finite_compare(a: sp.Expr, b: sp.Expr) -> int:
    return compare_exact_reals(a, b)


def _best_candidates(
    candidates: Iterable[_Candidate], *, kind: str
) -> tuple[sp.Expr, tuple[Mapping[sp.Symbol, sp.Expr], ...], bool]:
    cand = list(candidates)
    if not cand:
        raise ValueError("no feasible candidate points were found")
    best = cand[0]
    for item in cand[1:]:
        cmp = _finite_compare(item.value, best.value)
        if (kind == "min" and cmp < 0) or (kind == "max" and cmp > 0):
            best = item
    tied = [item for item in cand if _finite_compare(item.value, best.value) == 0]
    attained = any(item.attained for item in tied)
    points: list[Mapping[sp.Symbol, sp.Expr]] = []
    seen: set[tuple[tuple[str, str], ...]] = set()
    for item in tied:
        if item.attained and item.point is not None:
            key = tuple(
                sorted((sp.sstr(k), sp.sstr(sp.simplify(v))) for k, v in item.point.items())
            )
            if key not in seen:
                points.append(dict(item.point))
                seen.add(key)
    return sp.simplify(best.value), tuple(points), attained


def _limits_for_interval(
    objective: sp.Expr, variable: sp.Symbol, lo: sp.Expr, hi: sp.Expr
) -> tuple[_Candidate, ...]:
    out: list[_Candidate] = []
    if lo != -sp.oo:
        try:
            out.append(
                _Candidate(sp.simplify(sp.limit(objective, variable, lo, dir="+")), None, False)
            )
        except _EXPECTED_ERRORS:
            pass
    else:
        try:
            out.append(_Candidate(sp.simplify(sp.limit(objective, variable, -sp.oo)), None, False))
        except _EXPECTED_ERRORS:
            pass
    if hi != sp.oo:
        try:
            out.append(
                _Candidate(sp.simplify(sp.limit(objective, variable, hi, dir="-")), None, False)
            )
        except _EXPECTED_ERRORS:
            pass
    else:
        try:
            out.append(_Candidate(sp.simplify(sp.limit(objective, variable, sp.oo)), None, False))
        except _EXPECTED_ERRORS:
            pass
    return tuple(out)


def _univariate_candidates(
    objective: sp.Expr, condition: sp.Expr, variable: sp.Symbol
) -> tuple[_Candidate, ...]:
    """Enumerate exact stationary, endpoint, and limiting candidates in 1D."""

    candidates: list[_Candidate] = []
    intervals = _one_dimensional_intervals(condition, variable, None)
    cuts: set[str] = set()
    cut_values: list[sp.Expr] = []
    for poly in _relational_polynomials(condition):
        if poly.free_symbols <= {variable}:
            for root in _finite_real_roots(poly, variable):
                key = sp.sstr(root)
                if key not in cuts:
                    cuts.add(key)
                    cut_values.append(root)
    # Zero-dimensional feasible pieces and closed interval endpoints.
    for point_value in cut_values:
        point = {variable: sp.simplify(point_value)}
        if _is_feasible(condition, point):
            candidates.append(_Candidate(sp.simplify(objective.subs(point)), point, True))
    derivative = sp.diff(objective, variable)
    derivative_roots: tuple[sp.Expr, ...]
    if derivative == 0:
        derivative_roots = ()
    else:
        try:
            derivative_roots = _finite_real_roots(derivative, variable)
        except _EXPECTED_ERRORS:
            derivative_roots = ()
    for lo, hi in intervals:
        candidates.extend(_limits_for_interval(objective, variable, lo, hi))
        for root in derivative_roots:
            if (lo == -sp.oo or _finite_compare(lo, root) < 0) and (
                hi == sp.oo or _finite_compare(root, hi) < 0
            ):
                point = {variable: sp.simplify(root)}
                if _is_feasible(condition, point):
                    candidates.append(_Candidate(sp.simplify(objective.subs(point)), point, True))
        # Include a simple interior point for constant objectives.
        if derivative == 0:
            if lo == -sp.oo and hi == sp.oo:
                sample = sp.Integer(0)
            elif lo == -sp.oo:
                sample = hi - 1
            elif hi == sp.oo:
                sample = lo + 1
            else:
                sample = sp.simplify((lo + hi) / 2)
            point = {variable: sample}
            if _is_feasible(condition, point):
                candidates.append(_Candidate(sp.simplify(objective.subs(point)), point, True))
    if condition is sp.true or isinstance(condition, BooleanTrue):
        # No relational polynomials means the whole real line.
        if not intervals:
            candidates.extend(_limits_for_interval(objective, variable, -sp.oo, sp.oo))
            for root in derivative_roots:
                point = {variable: sp.simplify(root)}
                candidates.append(_Candidate(sp.simplify(objective.subs(point)), point, True))
    return tuple(candidates)


def _solutions_to_points(
    solutions: object, variables: Sequence[sp.Symbol]
) -> tuple[dict[sp.Symbol, sp.Expr], ...]:
    points: list[dict[sp.Symbol, sp.Expr]] = []
    if isinstance(solutions, dict):
        solutions = [solutions]
    if isinstance(solutions, (list, tuple, set)):
        for sol in solutions:
            if isinstance(sol, dict):
                if all(var in sol for var in variables):
                    points.append({var: sp.simplify(sol[var]) for var in variables})
            elif isinstance(sol, (list, tuple)) and len(sol) >= len(variables):
                points.append({var: sp.simplify(sol[i]) for i, var in enumerate(variables)})
    return tuple(points)


def _real_point(point: Mapping[sp.Symbol, sp.Expr]) -> bool:
    for value in point.values():
        value = sp.simplify(value)
        if value.free_symbols:
            return False
        if value.has(sp.I) and sp.simplify(sp.im(value)) != 0:
            return False
        if value.is_real is False:
            return False
    return True


def _constraint_data(
    condition: sp.Expr, variables: Sequence[sp.Symbol]
) -> tuple[tuple[sp.Expr, ...], tuple[sp.Expr, ...], tuple[sp.Expr, ...]]:
    """Return equality, inequality-boundary, and disequality residuals.

    Non-strict inequalities contribute their zero sets as possible attained
    active boundaries. Strict inequalities are deliberately not enumerated as
    KKT active sets because their boundary points are infeasible; unattained
    extrema on such boundaries are handled by exact range certification.
    """

    equalities: list[sp.Expr] = []
    inequalities: list[sp.Expr] = []
    disequalities: list[sp.Expr] = []
    variable_set = set(variables)
    seen: set[tuple[str, str]] = set()
    for atom in _atoms(condition):
        if atom is sp.false:
            return (), (), ()
        residual, op = _relation_parts(atom)
        if not residual.free_symbols <= variable_set:
            continue
        if op == "!=":
            target = disequalities
        elif op == "==":
            target = equalities
        elif op in {"<", ">"}:
            # An attained feasible point can never lie on a strict boundary.
            continue
        else:
            target = inequalities
        key = (op if op in {"==", "!="} else "ineq", sp.sstr(sp.factor(residual)))
        if key not in seen:
            target.append(sp.expand(residual))
            seen.add(key)
    return tuple(equalities), tuple(inequalities), tuple(disequalities)


def _polynomial_problem(
    objective: sp.Expr, condition: sp.Expr, variables: Sequence[sp.Symbol]
) -> bool:
    if not variables:
        return not (sp.sympify(objective).free_symbols or sp.sympify(condition).free_symbols)
    try:
        sp.Poly(sp.expand(objective), *variables, domain=sp.QQ)
        for atom in _atoms(condition):
            if atom is sp.false:
                continue
            residual, _ = _relation_parts(atom)
            sp.Poly(residual, *variables, domain=sp.QQ)
    except (PolynomialError, ValueError, TypeError):
        return False
    return True


def _solve_exact_equations(
    equations: Sequence[sp.Expr],
    solve_variables: Sequence[sp.Symbol],
    original_variables: Sequence[sp.Symbol],
    condition: sp.Expr,
) -> tuple[dict[sp.Symbol, sp.Expr], ...]:
    """Solve a zero-dimensional polynomial system exactly, preferring RUR."""

    equations = tuple(sp.expand(eq) for eq in equations if sp.expand(eq) != 0)
    if not equations or not solve_variables:
        return ()
    points: list[dict[sp.Symbol, sp.Expr]] = []
    try:
        from .solve.zero_dimensional import is_zero_dimensional, solve_zero_dimensional_system

        if is_zero_dimensional(equations, solve_variables):
            result = solve_zero_dimensional_system(
                equations, vars=tuple(solve_variables), backend="rur", real=True
            )
            for assignment in result.assignments:
                point = {
                    var: sp.simplify(assignment[var])
                    for var in original_variables
                    if var in assignment
                }
                if (
                    len(point) == len(original_variables)
                    and _real_point(point)
                    and _is_feasible(condition, point)
                ):
                    points.append(point)
            return tuple(points)
    except _EXPECTED_ERRORS:
        pass
    try:
        raw = sp.solve(equations, tuple(solve_variables), dict=True)
    except _EXPECTED_ERRORS:
        return ()
    for point in _solutions_to_points(raw, original_variables):
        if _real_point(point) and _is_feasible(condition, point):
            points.append(point)
    return tuple(points)


def _project_kkt_locus(
    equations: Sequence[sp.Expr],
    multipliers: Sequence[sp.Symbol],
    variables: Sequence[sp.Symbol],
) -> tuple[sp.Expr, ...]:
    """Eliminate KKT multipliers, returning equations on original variables."""

    if not multipliers:
        return tuple(sp.expand(eq) for eq in equations if sp.expand(eq) != 0)
    all_vars = (*multipliers, *variables)
    try:
        basis = sp.groebner(tuple(equations), *all_vars, order="lex", domain=sp.QQ)
    except (PolynomialError, ValueError, TypeError):
        return ()
    multiplier_set = set(multipliers)
    projected = [
        sp.expand(poly.as_expr())
        for poly in basis.polys
        if not (poly.as_expr().free_symbols & multiplier_set)
    ]
    return tuple(dict.fromkeys(expr for expr in projected if expr != 0))


def _reduce_linear_equalities(
    objective: sp.Expr,
    condition: sp.Expr,
    variables: tuple[sp.Symbol, ...],
) -> tuple[sp.Expr, sp.Expr, tuple[sp.Symbol, ...], dict[sp.Symbol, sp.Expr]]:
    """Eliminate variables from equations that are globally linear with constant coefficient."""

    obj = sp.expand(objective)
    cond = condition
    remaining = list(variables)
    substitutions: dict[sp.Symbol, sp.Expr] = {}
    changed = True
    while changed:
        changed = False
        try:
            equalities, _, _ = _constraint_data(cond, tuple(remaining))
        except (TypeError, NotImplementedError):
            break
        for eq in equalities:
            for var in tuple(remaining):
                try:
                    poly = sp.Poly(sp.expand(eq), var)
                except (PolynomialError, ValueError, TypeError):
                    continue
                if poly.degree() != 1:
                    continue
                coefficient = sp.expand(poly.coeff_monomial(var))
                rest = sp.expand(poly.coeff_monomial(1))
                # Only divide by a coefficient independent of every remaining
                # optimization variable; otherwise a hidden coefficient-zero
                # parameter stratum would be lost.
                if coefficient == 0 or (coefficient.free_symbols & set(remaining)):
                    continue
                replacement = sp.cancel(-rest / coefficient)
                if var in replacement.free_symbols:
                    continue
                substitutions[var] = sp.simplify(replacement.subs(substitutions))
                obj = sp.expand(obj.subs(var, replacement))
                cond = sp.simplify(cond.subs(var, replacement))
                remaining.remove(var)
                changed = True
                break
            if changed:
                break
    # Compose substitutions so lifted points depend only on retained variables.
    for key in reversed(tuple(substitutions)):
        substitutions[key] = sp.simplify(substitutions[key].subs(substitutions))
    return obj, cond, tuple(remaining), substitutions


def _lift_reduced_points(
    points: Sequence[Mapping[sp.Symbol, sp.Expr]],
    substitutions: Mapping[sp.Symbol, sp.Expr],
    original_variables: Sequence[sp.Symbol],
) -> tuple[Mapping[sp.Symbol, sp.Expr], ...]:
    lifted: list[Mapping[sp.Symbol, sp.Expr]] = []
    for point in points:
        assignment = dict(point)
        pending = dict(substitutions)
        for _ in range(len(pending) + 1):
            progress = False
            for var, expr in list(pending.items()):
                value = sp.simplify(expr.subs(assignment))
                if not (value.free_symbols & set(original_variables)):
                    assignment[var] = value
                    del pending[var]
                    progress = True
            if not pending or not progress:
                break
        if all(var in assignment for var in original_variables):
            lifted.append({var: sp.simplify(assignment[var]) for var in original_variables})
    return tuple(lifted)


def _kkt_candidate_points(
    objective: sp.Expr,
    condition: sp.Expr,
    variables: tuple[sp.Symbol, ...],
    equalities: tuple[sp.Expr, ...],
    inequalities: tuple[sp.Expr, ...],
) -> tuple[dict[sp.Symbol, sp.Expr], ...]:
    """Enumerate exact zero-dimensional KKT/singular active-set candidates."""

    points: list[dict[sp.Symbol, sp.Expr]] = []
    seen: set[tuple[tuple[str, str], ...]] = set()
    subsets = _pruned_active_subsets(equalities, inequalities, variables, condition)
    for subset in subsets:
        active = tuple(equalities) + tuple(subset)
        equations, multipliers = _kkt_system(objective, variables, active)
        solve_vars = (*variables, *multipliers)
        dimension = polynomial_locus_dimension(equations, solve_vars)
        if dimension == 0:
            for point in _solve_exact_equations(equations, solve_vars, variables, condition):
                key = tuple(sorted((sp.sstr(k), sp.sstr(sp.simplify(v))) for k, v in point.items()))
                if key not in seen:
                    points.append(point)
                    seen.add(key)

        if active:
            minors = _jacobian_rank_deficiency_equations(active, variables)
            singular_eqs = (*active, *minors)
            singular_dimension = (
                polynomial_locus_dimension(singular_eqs, variables) if minors else None
            )
            if minors and singular_dimension == 0:
                for point in _solve_exact_equations(singular_eqs, variables, variables, condition):
                    key = tuple(
                        sorted((sp.sstr(k), sp.sstr(sp.simplify(v))) for k, v in point.items())
                    )
                    if key not in seen:
                        points.append(point)
                        seen.add(key)

        if len(active) >= len(variables) and polynomial_locus_dimension(active, variables) == 0:
            for point in _solve_exact_equations(active, variables, variables, condition):
                key = tuple(sorted((sp.sstr(k), sp.sstr(sp.simplify(v))) for k, v in point.items()))
                if key not in seen:
                    points.append(point)
                    seen.add(key)
    return tuple(points)


def _positive_dimensional_kkt_candidates(
    objective: sp.Expr,
    condition: sp.Expr,
    variables: tuple[sp.Symbol, ...],
    equalities: tuple[sp.Expr, ...],
    inequalities: tuple[sp.Expr, ...],
    *,
    kind: str,
    policy: OptimizationCertificationPolicy,
    recursion_depth: int,
    visited_loci: frozenset[str],
) -> tuple[_Candidate, ...]:
    """Optimize recursively over positive-dimensional projected KKT loci."""

    if recursion_depth >= policy.recursion_limit:
        return ()
    out: list[_Candidate] = []
    seen_loci: set[str] = set()
    for subset in _pruned_active_subsets(equalities, inequalities, variables, condition):
        active = tuple(equalities) + tuple(subset)
        equations, multipliers = _kkt_system(objective, variables, active)
        dimension = polynomial_locus_dimension(equations, (*variables, *multipliers))
        if dimension is None or dimension <= 0:
            continue
        projected = _project_kkt_locus(equations, multipliers, variables)
        if not projected:
            continue
        projected_dimension = polynomial_locus_dimension(projected, variables)
        if (
            projected_dimension is None
            or projected_dimension <= 0
            or projected_dimension >= len(variables)
        ):
            continue
        key = "|".join(sorted(sp.srepr(sp.factor(eq)) for eq in projected))
        if key in visited_loci or key in seen_loci:
            continue
        seen_loci.add(key)
        locus_condition = sp.And(condition, *(sp.Eq(eq, 0) for eq in projected), evaluate=False)
        reduced_obj, reduced_cond, reduced_vars, substitutions = _reduce_linear_equalities(
            objective, locus_condition, variables
        )
        try:
            if len(reduced_vars) < len(variables):
                result = _optimize_conjunction(
                    reduced_obj,
                    reduced_cond,
                    reduced_vars,
                    kind=kind,
                    policy=policy,
                    recursion_depth=recursion_depth + 1,
                    visited_loci=visited_loci | frozenset({key}),
                    allow_equality_reduction=True,
                )
                lifted = _lift_reduced_points(result.points, substitutions, variables)
                out.append(_Candidate(result.value, lifted[0] if lifted else None, result.attained))
            else:
                # No safe coordinate elimination is available.  A complete
                # range computation on the lower-dimensional locus is the
                # conservative terminal step for this recursive branch.
                cert = _certify_optimum_by_range(
                    objective, locus_condition, variables, kind=kind, policy=policy
                )
                if cert is not None:
                    value, attained = cert
                    out.append(_Candidate(value, None, attained))
        except (NotImplementedError, ValueError, TypeError, ArithmeticError, PolynomialError):
            continue
    return tuple(out)


def _multivariate_exact_candidates(
    objective: sp.Expr,
    condition: sp.Expr,
    variables: tuple[sp.Symbol, ...],
    *,
    kind: str,
    policy: OptimizationCertificationPolicy,
    recursion_depth: int = 0,
    visited_loci: frozenset[str] = frozenset(),
) -> tuple[_Candidate, ...]:
    equalities, inequalities, _ = _constraint_data(condition, variables)
    points = _kkt_candidate_points(objective, condition, variables, equalities, inequalities)
    candidates = [_Candidate(sp.simplify(objective.subs(point)), point, True) for point in points]
    candidates.extend(
        _positive_dimensional_kkt_candidates(
            objective,
            condition,
            variables,
            equalities,
            inequalities,
            kind=kind,
            policy=policy,
            recursion_depth=recursion_depth,
            visited_loci=visited_loci,
        )
    )
    return tuple(candidates)


def _optimization_candidates(
    objective: sp.Expr,
    constraints: sp.Expr,
    variables: tuple[sp.Symbol, ...],
    *,
    kind: str,
    policy: OptimizationCertificationPolicy,
    recursion_depth: int = 0,
    visited_loci: frozenset[str] = frozenset(),
) -> tuple[_Candidate, ...]:
    if constraints is sp.false or isinstance(constraints, BooleanFalse):
        return ()
    if len(variables) == 1:
        return _univariate_candidates(objective, constraints, variables[0])
    return _multivariate_exact_candidates(
        objective,
        constraints,
        variables,
        kind=kind,
        policy=policy,
        recursion_depth=recursion_depth,
        visited_loci=visited_loci,
    )


def _certify_candidate_by_qe(
    objective: sp.Expr,
    condition: sp.Expr,
    variables: tuple[sp.Symbol, ...],
    value: sp.Expr,
    *,
    kind: str,
) -> bool:
    """Prove that no feasible point has objective strictly better than ``value``."""

    better = objective < value if kind == "min" else objective > value
    sentence = sp.And(condition, better, evaluate=False)
    try:
        result = qe_by_complete_cad(
            variables,
            tuple(("exists", var) for var in variables),
            parse_formula(sentence),
        )
    except (NotImplementedError, ValueError, TypeError, ArithmeticError, PolynomialError):
        return False
    return result.truth_value is False


def _range_certification_cost(
    objective: sp.Expr,
    condition: sp.Expr,
    variables: tuple[sp.Symbol, ...],
) -> int:
    """Conservative symbolic estimate for complete image-CAD cost."""

    if not variables:
        return 0
    degrees: list[int] = []
    polynomial_count = 1
    try:
        degrees.append(max(1, sp.Poly(objective, *variables).total_degree()))
        for atom in _atoms(condition):
            if atom is sp.false:
                continue
            residual, _ = _relation_parts(atom)
            degrees.append(max(1, sp.Poly(residual, *variables).total_degree()))
            polynomial_count += 1
    except _EXPECTED_ERRORS:
        return 10**9
    max_degree = max(degrees, default=1)
    # Image CAD introduces one additional value variable.  This estimate is
    # intentionally monotone and coarse; it is a guardrail, not a complexity proof.
    return int(polynomial_count * (max_degree + 1) ** (len(variables) + 1))


def _certify_optimum_by_range(
    objective: sp.Expr,
    condition: sp.Expr,
    variables: tuple[sp.Symbol, ...],
    *,
    kind: str,
    policy: OptimizationCertificationPolicy,
) -> tuple[sp.Expr, bool] | None:
    """Use complete CAD image computation when allowed by the cost policy."""

    cost = _range_certification_cost(objective, condition, variables)
    if policy.mode == "candidate":
        return None
    if policy.mode == "auto" and cost > policy.range_cost_limit:
        return None
    try:
        result = function_range(
            objective,
            condition,
            variables,
            value_symbol=sp.Symbol("_semialg_opt_value", real=True),
            method="cad",
            return_result=True,
        )
    except (NotImplementedError, ValueError, TypeError, ArithmeticError, PolynomialError):
        return None
    if not isinstance(result, FunctionRangeResult):
        return None
    if kind == "min" and result.infimum is not None:
        return sp.simplify(result.infimum), bool(result.minimum_attained)
    if kind == "max" and result.supremum is not None:
        return sp.simplify(result.supremum), bool(result.maximum_attained)
    return None


def _optimize_conjunction(
    objective: sp.Expr,
    condition: sp.Expr,
    variables: tuple[sp.Symbol, ...],
    *,
    kind: str,
    policy: OptimizationCertificationPolicy,
    recursion_depth: int = 0,
    visited_loci: frozenset[str] = frozenset(),
    allow_equality_reduction: bool = True,
) -> OptimizationResult:
    """Optimize one conjunction with exact KKT and certification machinery.

    Equalities are simplified before active-set enumeration. The routine then
    combines finite exact candidates, positive-dimensional critical loci, and
    CAD-based global certification while preserving attainment information.
    """

    if condition is sp.false or isinstance(condition, BooleanFalse):
        raise ValueError("optimization domain is empty")
    if not _polynomial_problem(objective, condition, variables):
        candidates = _optimization_candidates(
            objective,
            condition,
            variables,
            kind=kind,
            policy=policy,
            recursion_depth=recursion_depth,
            visited_loci=visited_loci,
        )
        value, points, attained = _best_candidates(candidates, kind=kind)
        return OptimizationResult(
            objective,
            variables,
            value,
            points,
            attained,
            kind,
            "critical_point_enumeration",
            {"candidate_count": len(candidates)},
            False,
        )

    if allow_equality_reduction and variables:
        reduced_obj, reduced_cond, reduced_vars, substitutions = _reduce_linear_equalities(
            objective, condition, variables
        )
        if len(reduced_vars) < len(variables):
            reduced_result = _optimize_conjunction(
                reduced_obj,
                reduced_cond,
                reduced_vars,
                kind=kind,
                policy=policy,
                recursion_depth=recursion_depth + 1,
                visited_loci=visited_loci,
                allow_equality_reduction=True,
            )
            seed_points = reduced_result.points
            if not seed_points and reduced_result.attained and not reduced_vars:
                seed_points = ({},)
            lifted_points = _lift_reduced_points(seed_points, substitutions, variables)
            return OptimizationResult(
                objective,
                variables,
                reduced_result.value,
                lifted_points,
                reduced_result.attained,
                kind,
                "equality_reduction+" + reduced_result.method,
                {
                    **dict(reduced_result.diagnostics),
                    "eliminated_variables": tuple(sp.sstr(v) for v in substitutions),
                },
                reduced_result.certified,
            )

    candidates = _optimization_candidates(
        objective,
        condition,
        variables,
        kind=kind,
        policy=policy,
        recursion_depth=recursion_depth,
        visited_loci=visited_loci,
    )
    candidate_value = None
    candidate_points: tuple[Mapping[sp.Symbol, sp.Expr], ...] = ()
    candidate_attained = False
    if candidates:
        candidate_value, candidate_points, candidate_attained = _best_candidates(
            candidates, kind=kind
        )

    if candidate_value is not None and _certify_candidate_by_qe(
        objective, condition, variables, candidate_value, kind=kind
    ):
        return OptimizationResult(
            objective,
            variables,
            candidate_value,
            candidate_points,
            candidate_attained,
            kind,
            "exact_kkt_active_set+cad_decision_certificate",
            {
                "candidate_count": len(candidates),
                "constraints": sp.sstr(condition),
                "global_certificate": "complete_cad_no_better_point",
            },
            True,
        )

    certificate = _certify_optimum_by_range(
        objective, condition, variables, kind=kind, policy=policy
    )
    if certificate is not None:
        value, attained = certificate
        points = (
            candidate_points
            if candidate_value is not None
            and _finite_compare(candidate_value, value) == 0
            and attained
            else ()
        )
        return OptimizationResult(
            objective,
            variables,
            value,
            points,
            attained,
            kind,
            "exact_kkt_active_set+cad_range_certificate",
            {
                "candidate_count": len(candidates),
                "constraints": sp.sstr(condition),
                "global_certificate": "complete_cad_function_range",
            },
            True,
        )
    if candidate_value is None:
        raise NotImplementedError(
            "exact candidate enumeration produced no finite candidate and CAD certification was unavailable"
        )
    return OptimizationResult(
        objective,
        variables,
        candidate_value,
        candidate_points,
        candidate_attained,
        kind,
        "exact_kkt_active_set",
        {
            "candidate_count": len(candidates),
            "constraints": sp.sstr(condition),
            "global_certificate": "candidate_exhaustion_without_cad_range",
        },
        False,
    )


def _normalize_parameters_for_problem(
    parameters: Sequence[sp.Symbol | str],
    *expressions: sp.Expr,
) -> tuple[sp.Symbol, ...]:
    from .symbol_resolution import normalize_variables

    return normalize_variables(
        parameters,
        context=expressions,
        append_context_symbols=False,
    )


def _parameter_guards(
    condition: sp.Expr,
    variables: tuple[sp.Symbol, ...],
    parameters: tuple[sp.Symbol, ...],
) -> tuple[sp.Expr, tuple[tuple[sp.Expr, Mapping[sp.Symbol, sp.Expr]], ...]]:
    from .cad.cells import extract_cylindrical_solution
    from .parameters import solvability_conditions

    parameter_domain = sp.simplify(solvability_conditions(condition, variables, parameters))
    if parameter_domain is sp.false or parameter_domain == sp.false:
        return sp.false, ()
    try:
        solution = extract_cylindrical_solution(parameter_domain, parameters, selected_only=True)
    except (NotImplementedError, ValueError, TypeError, ArithmeticError, PolynomialError):
        solution = None
    if solution is None or not solution.cells:
        sample = {param: sp.Integer(0) for param in parameters}
        return parameter_domain, ((parameter_domain, sample),)
    return parameter_domain, tuple(
        (cell.as_formula(closed=False), cell.sample_point()) for cell in solution.cells
    )


def _parametric_range_relation(
    expression: sp.Expr,
    condition: sp.Expr,
    variables: tuple[sp.Symbol, ...],
    parameters: tuple[sp.Symbol, ...],
    value_symbol: sp.Symbol,
) -> sp.Expr:
    relation, domain_constraint = _relation_for_function_graph(expression, value_symbol)
    image_formula = sp.And(condition, domain_constraint, relation, evaluate=False)
    all_vars = tuple(dict.fromkeys((*parameters, value_symbol, *variables)))
    result = qe_by_complete_cad(
        all_vars,
        tuple(("exists", variable) for variable in variables),
        parse_formula(image_formula),
        free_variables=(*parameters, value_symbol),
    )
    return sp.simplify(result.formula)


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
    bound_vars = tuple(sp.Dummy(f"_semialg_bound_{i}", real=True) for i in range(len(variables)))
    tight_vars = tuple(sp.Dummy(f"_semialg_tight_{i}", real=True) for i in range(len(variables)))
    threshold = sp.Dummy("_semialg_threshold", real=True)
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
    quantified_vars = tuple(
        sp.Dummy(f"_semialg_range_{i}", real=True) for i in range(len(variables))
    )
    substitutions = dict(zip(variables, quantified_vars, strict=True))
    specialized_expression = expression.xreplace(substitutions)
    specialized_condition = condition.xreplace(substitutions)
    graph_formula, aux_symbols = _graph_formula_for_expression(
        specialized_expression, value_symbol, [0]
    )
    formula = sp.And(specialized_condition, graph_formula, evaluate=False)
    quantified = (*quantified_vars, *aux_symbols)
    return formula, tuple(("exists", var) for var in quantified)


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
):
    """Build guarded exact optimization relations over parameter strata."""

    from .conditional import ConditionalBranch, conditional_result

    if not _polynomial_problem(objective, constraints, (*parameters, *variables)):
        raise NotImplementedError(
            "parameter-stratified optimization currently requires polynomial data"
        )
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
            )
        except _EXPECTED_ERRORS:
            pass
        value = ParametricOptimizationResult(
            objective,
            constraints,
            variables,
            parameters,
            value_symbol,
            kind,
            sp.And(guard, optimum_relation, evaluate=False),
            optimum_quantifiers,
            sample_result if isinstance(sample_result, OptimizationResult) else None,
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
):
    from .conditional import ConditionalBranch, conditional_result

    relation, quantifiers = _parametric_range_definition(
        expression, constraints, variables, value_symbol
    )
    parameter_domain, guards = _parameter_guards(constraints, variables, parameters)
    branches = []
    for guard, sample in guards:
        guarded_formula = sp.And(guard, relation, evaluate=False)
        value = ParametricFunctionRangeResult(
            expression,
            constraints,
            variables,
            parameters,
            value_symbol,
            guarded_formula,
            quantifiers,
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


@with_computation_context
def _optimize(
    objective: sp.Expr,
    constraints: FormulaLike | Iterable[FormulaLike] | None,
    variables: Sequence[sp.Symbol | str] | None,
    *,
    kind: str,
    domain: str = "reals",
    return_result: bool = True,
    certification: Literal["auto", "complete", "candidate"] = "auto",
    range_cost_limit: int = 2500,
    recursion_limit: int = 4,
    parameters: Sequence[sp.Symbol | str] | None = None,
    return_stratified: bool = False,
) -> OptimizationResult | sp.Expr | object:
    """Shared exact implementation for minimization and maximization APIs.

    Boolean domains are decomposed into boundedly many conjunctions. Each
    branch is optimized independently and exact branch bounds are compared
    before the combined result and certificate status are assembled.
    """

    if domain.lower() not in {"real", "reals", "r", "rr"}:
        raise NotImplementedError(
            "semialgebraic optimization currently supports only the real domain"
        )
    policy = OptimizationCertificationPolicy(certification, range_cost_limit, recursion_limit)
    obj = sp.sympify(objective)
    condition = _normalize_formula(constraints)
    vars_ = _normalize_variables(variables, sp.Tuple(condition, obj))
    if not vars_:
        value = sp.simplify(obj)
        result = OptimizationResult(obj, vars_, value, (), True, kind, "constant", {}, True)
        return result if return_result else result.value

    expansion = bounded_dnf_branches(condition, max_branches=32)
    if not expansion.complete:
        raise NotImplementedError(
            "optimization Boolean expansion exceeded the bounded branch limit"
        )
    branch_results: list[OptimizationResult] = []
    for branch in expansion.branches:
        branch_condition = sp.And(
            *[piece for piece in branch if piece not in (True, sp.true)], evaluate=False
        )
        if any(piece in (False, sp.false) for piece in branch):
            continue
        try:
            branch_results.append(
                _optimize_conjunction(obj, branch_condition, vars_, kind=kind, policy=policy)
            )
        except ValueError:
            continue
    if not branch_results:
        raise ValueError("optimization domain is empty or unsupported")

    best = branch_results[0]
    for result in branch_results[1:]:
        cmp = _finite_compare(result.value, best.value)
        if (kind == "min" and cmp < 0) or (kind == "max" and cmp > 0):
            best = result
    tied = [result for result in branch_results if _finite_compare(result.value, best.value) == 0]
    points: list[Mapping[sp.Symbol, sp.Expr]] = []
    for result in tied:
        points.extend(result.points)
    attained = any(result.attained for result in tied)
    result = OptimizationResult(
        obj,
        vars_,
        best.value,
        tuple(points),
        attained,
        kind,
        best.method if len(branch_results) == 1 else "exact_branchwise_optimization",
        {
            **dict(best.diagnostics),
            "branch_count": len(branch_results),
            "all_branches_certified": all(item.certified for item in branch_results),
        },
        all(item.certified for item in branch_results),
    )
    return result if return_result else result.value


@with_computation_context
def semialgebraic_minimize(
    objective: sp.Expr,
    constraints: FormulaLike | Iterable[FormulaLike] | None = None,
    variables: Sequence[sp.Symbol | str] | None = None,
    *,
    domain: str = "reals",
    return_result: bool = True,
    certification: Literal["auto", "complete", "candidate"] = "auto",
    range_cost_limit: int = 2500,
    recursion_limit: int = 4,
    parameters: Sequence[sp.Symbol | str] | None = None,
    return_stratified: bool = False,
) -> OptimizationResult | sp.Expr | object:
    """Return an exact minimum/infimum for a polynomial semialgebraic problem.

    Multivariate polynomial problems use exact active-set/KKT enumeration,
    singular active-locus solving, RUR-backed zero-dimensional solving, and a
    complete-CAD image certificate when available.
    """

    if return_stratified:
        obj = sp.sympify(objective)
        condition = _normalize_formula(constraints)
        params = _normalize_parameters_for_problem(parameters or (), obj, condition)
        if not params:
            raise ValueError("return_stratified=True requires at least one parameter")
        vars_ = _normalize_variables(variables, sp.Tuple(condition, obj))
        vars_ = tuple(var for var in vars_ if var not in set(params))
        return _stratified_optimization(
            obj,
            condition,
            vars_,
            params,
            kind="min",
            domain=domain,
            certification=certification,
            range_cost_limit=range_cost_limit,
            recursion_limit=recursion_limit,
        )
    return _optimize(
        objective,
        constraints,
        variables,
        kind="min",
        domain=domain,
        return_result=return_result,
        certification=certification,
        range_cost_limit=range_cost_limit,
        recursion_limit=recursion_limit,
    )


@with_computation_context
def semialgebraic_maximize(
    objective: sp.Expr,
    constraints: FormulaLike | Iterable[FormulaLike] | None = None,
    variables: Sequence[sp.Symbol | str] | None = None,
    *,
    domain: str = "reals",
    return_result: bool = True,
    certification: Literal["auto", "complete", "candidate"] = "auto",
    range_cost_limit: int = 2500,
    recursion_limit: int = 4,
    parameters: Sequence[sp.Symbol | str] | None = None,
    return_stratified: bool = False,
) -> OptimizationResult | sp.Expr | object:
    """Return an exact maximum/supremum for a polynomial semialgebraic problem."""

    if return_stratified:
        obj = sp.sympify(objective)
        condition = _normalize_formula(constraints)
        params = _normalize_parameters_for_problem(parameters or (), obj, condition)
        if not params:
            raise ValueError("return_stratified=True requires at least one parameter")
        vars_ = _normalize_variables(variables, sp.Tuple(condition, obj))
        vars_ = tuple(var for var in vars_ if var not in set(params))
        return _stratified_optimization(
            obj,
            condition,
            vars_,
            params,
            kind="max",
            domain=domain,
            certification=certification,
            range_cost_limit=range_cost_limit,
            recursion_limit=recursion_limit,
        )
    return _optimize(
        objective,
        constraints,
        variables,
        kind="max",
        domain=domain,
        return_result=return_result,
        certification=certification,
        range_cost_limit=range_cost_limit,
        recursion_limit=recursion_limit,
    )


def _range_formula(
    symbol: sp.Symbol, lo: sp.Expr, hi: sp.Expr, lo_attained: bool, hi_attained: bool
) -> sp.Expr:
    pieces: list[sp.Expr] = []
    if lo != -sp.oo:
        pieces.append(symbol >= lo if lo_attained else symbol > lo)
    if hi != sp.oo:
        pieces.append(symbol <= hi if hi_attained else symbol < hi)
    if not pieces:
        return sp.true
    return sp.And(*pieces)


def _relation_for_function_graph(expr: sp.Expr, value_symbol: sp.Symbol) -> tuple[sp.Expr, sp.Expr]:
    """Return ``(relation, domain)`` for ``value_symbol == expr``.

    Polynomial expressions produce ``Eq(value_symbol - expr, 0)``. Rational
    expressions are handled by clearing denominators and adding the denominator
    nonzero condition, keeping the image problem semialgebraic.
    """

    numerator, denominator = sp.fraction(sp.together(expr))
    numerator = sp.expand(numerator)
    denominator = sp.expand(denominator)
    relation = sp.Eq(sp.expand(value_symbol * denominator - numerator), 0)
    domain = sp.true if denominator == 1 else sp.Ne(denominator, 0)
    return relation, domain


def _is_supported_sqrt(expr: sp.Expr) -> bool:
    """Return True for the square-root form supported by graph conversion."""

    return isinstance(expr, sp.Pow) and expr.exp == sp.Rational(1, 2)


def _is_semialgebraic_special(expr: sp.Expr) -> bool:
    """Return True for expression heads handled by semialgebraic graph conversion."""

    return (
        expr.func is sp.Abs
        or expr.func is sp.Max
        or expr.func is sp.Min
        or isinstance(expr, sp.Piecewise)
        or _is_supported_sqrt(expr)
    )


def _has_semialgebraic_special(expr: sp.Expr) -> bool:
    """Whether ``expr`` contains supported non-rational semialgebraic heads."""

    return any(_is_semialgebraic_special(item) for item in sp.preorder_traversal(expr))


def _fresh_graph_symbol(counter: list[int]) -> sp.Symbol:
    symbol = sp.Symbol(f"_semialg_graph_aux_{counter[0]}", real=True)
    counter[0] += 1
    return symbol


def _graph_formula_for_expression(
    expr: sp.Expr,
    target: sp.Symbol,
    counter: list[int],
) -> tuple[sp.Expr, tuple[sp.Symbol, ...]]:
    """Return a semialgebraic graph formula for ``target == expr``.

    The graph conversion supports common expression heads whose graphs are semialgebraic:
    ``Abs``, square roots, ``Min``, ``Max``, and simple ``Piecewise``. General
    arithmetic combinations are handled by replacing special subexpressions by
    auxiliary graph variables and then clearing denominators for the remaining
    rational expression.
    """

    expr = sp.sympify(expr)
    aux_symbols: list[sp.Symbol] = []

    if expr.func is sp.Abs:
        arg_value = _fresh_graph_symbol(counter)
        arg_formula, arg_aux = _graph_formula_for_expression(expr.args[0], arg_value, counter)
        aux_symbols.extend((arg_value, *arg_aux))
        formula = sp.And(
            arg_formula,
            target >= 0,
            sp.Or(sp.Eq(target, arg_value), sp.Eq(target, -arg_value)),
        )
        return formula, tuple(aux_symbols)

    if _is_supported_sqrt(expr):
        base_value = _fresh_graph_symbol(counter)
        base_formula, base_aux = _graph_formula_for_expression(expr.base, base_value, counter)
        aux_symbols.extend((base_value, *base_aux))
        formula = sp.And(
            base_formula,
            base_value >= 0,
            target >= 0,
            sp.Eq(target**2, base_value),
        )
        return formula, tuple(aux_symbols)

    if expr.func in {sp.Max, sp.Min}:
        arg_values: list[sp.Symbol] = []
        arg_formulas: list[sp.Expr] = []
        for arg in expr.args:
            arg_value = _fresh_graph_symbol(counter)
            arg_formula, arg_aux = _graph_formula_for_expression(arg, arg_value, counter)
            arg_values.append(arg_value)
            arg_formulas.append(arg_formula)
            aux_symbols.extend((arg_value, *arg_aux))
        branch_formulas: list[sp.Expr] = []
        for arg_value in arg_values:
            if expr.func is sp.Max:
                order = sp.And(*(arg_value >= other for other in arg_values))
            else:
                order = sp.And(*(arg_value <= other for other in arg_values))
            branch_formulas.append(sp.And(sp.Eq(target, arg_value), order))
        return sp.And(*arg_formulas, sp.Or(*branch_formulas)), tuple(aux_symbols)

    if isinstance(expr, sp.Piecewise):
        branches: list[sp.Expr] = []
        previous_conditions: list[sp.Expr] = []
        for branch_expr, branch_condition in expr.args:
            effective_condition = sp.And(
                branch_condition, *(sp.Not(cond) for cond in previous_conditions)
            )
            branch_formula, branch_aux = _graph_formula_for_expression(branch_expr, target, counter)
            aux_symbols.extend(branch_aux)
            branches.append(sp.And(effective_condition, branch_formula))
            previous_conditions.append(branch_condition)
        return sp.Or(*branches), tuple(aux_symbols)

    if not expr.args:
        relation, domain = _relation_for_function_graph(expr, target)
        return sp.And(domain, relation), ()

    constraints: list[sp.Expr] = []
    replacements: dict[sp.Expr, sp.Symbol] = {}
    all_aux: list[sp.Symbol] = []

    def replace_specials(item: sp.Expr) -> sp.Expr:
        item = sp.sympify(item)
        if _is_semialgebraic_special(item):
            if item in replacements:
                return replacements[item]
            aux = _fresh_graph_symbol(counter)
            formula, nested_aux = _graph_formula_for_expression(item, aux, counter)
            replacements[item] = aux
            constraints.append(formula)
            all_aux.extend((aux, *nested_aux))
            return aux
        if not item.args:
            return item
        new_args = tuple(replace_specials(arg) for arg in item.args)
        if new_args == item.args:
            return item
        return item.func(*new_args)

    transformed = replace_specials(expr)
    relation, domain = _relation_for_function_graph(transformed, target)
    aux_symbols.extend(all_aux)
    return sp.And(*(constraints + [domain, relation])), tuple(dict.fromkeys(aux_symbols))


def _range_formula_from_set(range_set: sp.Set, value_symbol: sp.Symbol) -> sp.Expr | None:
    """Convert a one-dimensional SymPy set into a range condition."""

    try:
        intervals = _intervals_from_range_set(range_set)
    except _EXPECTED_ERRORS:
        intervals = None
    if intervals is not None:
        pieces: list[sp.Expr] = []
        for left, right, left_closed, right_closed in intervals:
            if left == right:
                pieces.append(sp.Eq(value_symbol, left))
                continue
            lower = (
                sp.true
                if left == -sp.oo
                else (value_symbol >= left if left_closed else value_symbol > left)
            )
            upper = (
                sp.true
                if right == sp.oo
                else (value_symbol <= right if right_closed else value_symbol < right)
            )
            pieces.append(sp.And(lower, upper))
        return sp.Or(*pieces) if pieces else sp.false
    try:
        return range_set.as_relational(value_symbol)
    except _EXPECTED_ERRORS:
        return None


def _domain_set_from_condition(condition: sp.Expr, variable: sp.Symbol) -> sp.Set | None:
    """Return the one-variable domain represented by simple interval constraints."""

    if condition is sp.true or isinstance(condition, BooleanTrue):
        return sp.S.Reals
    if condition is sp.false or isinstance(condition, BooleanFalse):
        return sp.S.EmptySet
    atoms = condition.args if isinstance(condition, SymAnd) else (condition,)
    lower = -sp.oo
    upper = sp.oo
    left_open = False
    right_open = False
    for atom in atoms:
        if not getattr(atom, "is_Relational", False):
            return None
        expr, op = _relation_parts(atom)
        try:
            poly = sp.Poly(expr, variable)
        except _EXPECTED_ERRORS:
            return None
        if poly.degree() != 1:
            return None
        coeff = sp.simplify(poly.coeff_monomial(variable))
        const = sp.simplify(poly.eval(0))
        if coeff == 0:
            truth = bool(_relation_from_scalar(const, op))
            if not truth:
                return sp.S.EmptySet
            continue
        bound = sp.simplify(-const / coeff)
        if coeff.is_negative:
            op = {"<": ">", "<=": ">=", ">": "<", ">=": "<="}.get(op, op)
        if op in {">", ">="}:
            if lower == -sp.oo or _finite_compare(lower, bound) < 0:
                lower = bound
                left_open = op == ">"
            elif sp.simplify(lower - bound) == 0:
                left_open = left_open or op == ">"
        elif op in {"<", "<="}:
            if upper == sp.oo or _finite_compare(bound, upper) < 0:
                upper = bound
                right_open = op == "<"
            elif sp.simplify(upper - bound) == 0:
                right_open = right_open or op == "<"
        elif op == "==":
            lower = upper = bound
            left_open = right_open = False
        else:
            return None
    if lower != -sp.oo and upper != sp.oo and _finite_compare(lower, upper) > 0:
        return sp.S.EmptySet
    return sp.Interval(lower, upper, left_open=left_open, right_open=right_open)


def _relation_from_scalar(value: sp.Expr, op: str) -> bool:
    if op == "<":
        return bool(value < 0)
    if op == "<=":
        return bool(value <= 0)
    if op == ">":
        return bool(value > 0)
    if op == ">=":
        return bool(value >= 0)
    if op == "==":
        return bool(sp.simplify(value) == 0)
    if op == "!=":
        return bool(sp.simplify(value) != 0)
    raise ValueError(op)


def _interval_from_bounds(lo: sp.Expr, hi: sp.Expr) -> sp.Interval:
    if _finite_compare(lo, hi) > 0:
        lo, hi = hi, lo
    return sp.Interval(sp.simplify(lo), sp.simplify(hi))


def _poly_range_on_interval(
    expr: sp.Expr, variable: sp.Symbol, domain: sp.Interval
) -> sp.Set | None:
    """Return the range of a univariate polynomial on a closed interval."""

    if not isinstance(domain, sp.Interval):
        return None
    try:
        poly = sp.Poly(sp.expand(expr), variable)
    except _EXPECTED_ERRORS:
        return None
    if domain.start in (-sp.oo, sp.oo) or domain.end in (-sp.oo, sp.oo):
        return None
    candidates = [domain.start, domain.end]
    derivative = sp.diff(poly.as_expr(), variable)
    try:
        for root in sp.solve(sp.Eq(derivative, 0), variable):
            if root.is_real is False:
                continue
            if bool(root >= domain.start) and bool(root <= domain.end):
                candidates.append(root)
    except _EXPECTED_ERRORS:
        pass
    values = [sp.simplify(expr.subs(variable, item)) for item in candidates]
    lo = min(values, key=_bound_sort_key)
    hi = max(values, key=_bound_sort_key)
    return _interval_from_bounds(lo, hi)


def _range_via_abs_affine(
    expr: sp.Expr,
    condition: sp.Expr,
    variables: tuple[sp.Symbol, ...],
    value_symbol: sp.Symbol,
) -> sp.Expr | None:
    """Fast range for affine expressions in ``Abs(variable)``."""

    if len(variables) != 1:
        return None
    variable = variables[0]
    abs_expr = sp.Abs(variable)
    if not expr.has(abs_expr):
        return None
    domain = _domain_set_from_condition(condition, variable)
    if domain is None:
        return None
    u = sp.Symbol("_abs_value", real=True)
    transformed = expr.xreplace({abs_expr: u})
    try:
        poly = sp.Poly(sp.expand(transformed), u)
    except _EXPECTED_ERRORS:
        return None
    if poly.degree() != 1 or any(sym != u for sym in transformed.free_symbols):
        return None
    coeff = sp.simplify(poly.coeff_monomial(u))
    const = sp.simplify(poly.eval(0))
    if coeff == 0:
        return _range_formula_from_set(sp.FiniteSet(const), value_symbol)
    if domain is sp.S.Reals:
        abs_range = sp.Interval(0, sp.oo)
    elif isinstance(domain, sp.Interval):
        values = []
        if domain.start != -sp.oo:
            values.append(sp.Abs(domain.start))
        if domain.end != sp.oo:
            values.append(sp.Abs(domain.end))
        crosses_zero = domain.start <= 0 and domain.end >= 0
        lower = sp.Integer(0) if bool(crosses_zero) else min(values, key=_bound_sort_key)
        upper = (
            sp.oo
            if domain.start == -sp.oo or domain.end == sp.oo
            else max(values, key=_bound_sort_key)
        )
        abs_range = sp.Interval(lower, upper)
    else:
        return None
    left = sp.simplify(coeff * abs_range.start + const) if abs_range.start != -sp.oo else -sp.oo
    right = (
        sp.simplify(coeff * abs_range.end + const)
        if abs_range.end != sp.oo
        else (sp.oo if coeff.is_positive is not False else -sp.oo)
    )
    return _range_formula_from_set(_interval_from_bounds(left, right), value_symbol)


def _range_via_sqrt_quadratic(
    expr: sp.Expr,
    condition: sp.Expr,
    variables: tuple[sp.Symbol, ...],
    value_symbol: sp.Symbol,
) -> sp.Expr | None:
    """Fast range for square roots of nonnegative quadratic radicands."""

    if len(variables) != 1 or not _is_supported_sqrt(expr):
        return None
    variable = variables[0]
    domain = _domain_set_from_condition(condition, variable)
    if domain is None:
        return None
    try:
        base_poly = sp.Poly(sp.expand(expr.base), variable)
    except _EXPECTED_ERRORS:
        return None
    if base_poly.degree() > 2:
        return None
    coeff2 = sp.simplify(base_poly.coeff_monomial(variable**2))
    coeff1 = sp.simplify(base_poly.coeff_monomial(variable))
    coeff0 = sp.simplify(base_poly.coeff_monomial(1))
    if coeff2 == -1 and coeff1 == 0 and coeff0.is_real:
        natural = sp.Interval(-sp.sqrt(coeff0), sp.sqrt(coeff0))
    else:
        return None
    active_domain = natural if domain is sp.S.Reals else domain.intersect(natural)
    if not isinstance(active_domain, sp.Interval):
        return None
    base_range = _poly_range_on_interval(expr.base, variable, active_domain)
    if not isinstance(base_range, sp.Interval):
        return None
    upper = sp.sqrt(sp.simplify(base_range.end))
    return _range_formula_from_set(sp.Interval(0, upper), value_symbol)


def _range_via_sympy_calculus(
    expr: sp.Expr,
    condition: sp.Expr,
    variables: tuple[sp.Symbol, ...],
    value_symbol: sp.Symbol,
) -> sp.Expr | None:
    """Fast image computation for common univariate expression heads."""

    if len(variables) != 1:
        return None
    variable = variables[0]
    domain = _domain_set_from_condition(condition, variable)
    if domain is None:
        return None
    try:
        from sympy.calculus.util import function_range as sympy_range

        range_set = sympy_range(expr, variable, domain)
    except _EXPECTED_ERRORS:
        return None
    return _range_formula_from_set(range_set, value_symbol)


def _bound_sort_key(value: sp.Expr) -> float:
    if value == -sp.oo:
        return -float("inf")
    if value == sp.oo:
        return float("inf")
    return float(sp.N(value, 50))


def _range_via_max_min(
    expr: sp.Expr,
    condition: sp.Expr,
    variables: tuple[sp.Symbol, ...],
    value_symbol: sp.Symbol,
) -> sp.Expr | None:
    """Fast range for two-argument Max/Min with affine arguments."""

    if len(variables) != 1 or expr.func not in {sp.Max, sp.Min} or len(expr.args) != 2:
        return None
    variable = variables[0]
    domain = _domain_set_from_condition(condition, variable)
    if not isinstance(domain, sp.Interval):
        return None
    arg_ranges: list[sp.Set] = []
    for arg in expr.args:
        try:
            poly = sp.Poly(sp.expand(arg), variable)
        except _EXPECTED_ERRORS:
            return None
        if poly.degree() > 1:
            return None
        coeff = sp.simplify(poly.coeff_monomial(variable))
        const = sp.simplify(poly.eval(0))
        if coeff == 0:
            arg_ranges.append(sp.FiniteSet(const))
            continue
        left = -sp.oo if domain.start == -sp.oo else sp.simplify(coeff * domain.start + const)
        right = sp.oo if domain.end == sp.oo else sp.simplify(coeff * domain.end + const)
        lo, hi = (left, right) if _finite_compare(left, right) <= 0 else (right, left)
        arg_ranges.append(sp.Interval(lo, hi))
    if expr.func is sp.Max:
        lower = max((item.inf for item in arg_ranges), key=_bound_sort_key)
        upper = max((item.sup for item in arg_ranges), key=_bound_sort_key)
    else:
        lower = min((item.inf for item in arg_ranges), key=_bound_sort_key)
        upper = min((item.sup for item in arg_ranges), key=_bound_sort_key)
    return _range_formula_from_set(sp.Interval(lower, upper), value_symbol)


def _range_via_piecewise(
    expr: sp.Expr,
    condition: sp.Expr,
    variables: tuple[sp.Symbol, ...],
    value_symbol: sp.Symbol,
) -> sp.Expr | None:
    """Fast range for univariate Piecewise expressions by branch images."""

    if len(variables) != 1 or not isinstance(expr, sp.Piecewise):
        return None
    variable = variables[0]
    pieces: list[sp.Set] = []
    previous: list[sp.Expr] = []
    for branch_expr, branch_cond in expr.args:
        effective = sp.And(condition, branch_cond, *(sp.Not(item) for item in previous))
        previous.append(branch_cond)
        domain = _domain_set_from_condition(effective, variable)
        if domain is None or domain is sp.S.EmptySet:
            continue
        if isinstance(domain, sp.Union):
            domains = domain.args
        else:
            domains = (domain,)
        for branch_domain in domains:
            image = _poly_range_on_interval(branch_expr, variable, branch_domain)
            if image is None:
                return None
            pieces.append(image)
    if not pieces:
        return sp.false
    return _range_formula_from_set(sp.Union(*pieces), value_symbol)


def _try_direct_special_image(
    expr: sp.Expr,
    condition: sp.Expr,
    variables: tuple[sp.Symbol, ...],
    value_symbol: sp.Symbol,
) -> sp.Expr | None:
    """Use fast univariate range formulas before constructing a CAD image problem."""

    for helper in (
        _range_via_abs_affine,
        _range_via_sqrt_quadratic,
        _range_via_max_min,
        _range_via_piecewise,
        _range_via_sympy_calculus,
    ):
        result = helper(expr, condition, variables, value_symbol)
        if result is not None:
            return result
    return None


def _try_semialgebraic_graph_image(
    expr: sp.Expr,
    condition: sp.Expr,
    variables: tuple[sp.Symbol, ...],
    value_symbol: sp.Symbol,
) -> sp.Expr | None:
    """Use CAD/QE on a semialgebraic graph relation for semialgebraic graph ranges."""

    if not _has_semialgebraic_special(expr):
        return None
    try:
        graph_formula, aux_symbols = _graph_formula_for_expression(expr, value_symbol, [0])
        image_formula = sp.And(condition, graph_formula)
        elimination_variables = tuple(dict.fromkeys((*variables, *aux_symbols)))
        parsed = parse_formula(image_formula)
        result = qe_by_complete_cad(
            (value_symbol, *elimination_variables),
            tuple(("exists", variable) for variable in elimination_variables),
            parsed,
            free_variables=(value_symbol,),
        )
    except _EXPECTED_ERRORS:
        return None
    return sp.simplify(result.formula)


def _substitute_formula(expr: sp.Expr, substitution: Mapping[sp.Symbol, sp.Expr]) -> sp.Expr:
    if expr is sp.true or isinstance(expr, BooleanTrue):
        return sp.true
    if expr is sp.false or isinstance(expr, BooleanFalse):
        return sp.false
    return sp.simplify(expr.subs(substitution))


def _try_affine_univariate_image(
    expr: sp.Expr,
    condition: sp.Expr,
    variable: sp.Symbol,
    value_symbol: sp.Symbol,
) -> sp.Expr | None:
    """Fast exact image for an affine one-variable map.

    This covers important disconnected ranges without needing a full CAD run,
    e.g. the image of ``x`` over ``x <= -1 or x >= 1``.
    """

    try:
        poly = sp.Poly(sp.expand(expr), variable)
    except _EXPECTED_ERRORS:
        return None
    if poly.degree() != 1:
        return None
    coefficient = poly.coeff_monomial(variable)
    constant = poly.eval(0)
    if sp.simplify(coefficient) == 0:
        sample = {variable: sp.Integer(0)}
        return (
            sp.Eq(value_symbol, sp.simplify(constant))
            if _is_feasible(condition, sample)
            else sp.false
        )
    inverse = sp.simplify((value_symbol - constant) / coefficient)
    return _substitute_formula(condition, {variable: inverse})


def _try_solved_graph_image(
    expr: sp.Expr,
    condition: sp.Expr,
    variables: tuple[sp.Symbol, ...],
    value_symbol: sp.Symbol,
) -> sp.Expr | None:
    """Try exact image computation by solving the graph equation.

    This conservative helper supports univariate polynomial/rational maps when
    SymPy can solve the graph equation explicitly enough to substitute branches
    back into the domain condition.
    """

    if len(variables) != 1:
        return None
    variable = variables[0]
    affine = _try_affine_univariate_image(expr, condition, variable, value_symbol)
    if affine is not None:
        return sp.simplify(affine)
    numerator, denominator = sp.fraction(sp.together(expr))
    equation = sp.expand(value_symbol * denominator - numerator)
    try:
        solutions = sp.solve(equation, variable)
    except _EXPECTED_ERRORS:
        return None
    if not solutions:
        return None
    branches: list[sp.Expr] = []
    for solution in solutions:
        if variable in solution.free_symbols or solution.has(sp.I):
            return None
        if any(
            isinstance(pow_expr, sp.Pow) and pow_expr.exp.is_integer is False
            for pow_expr in solution.atoms(sp.Pow)
        ):
            return None
        branch_condition = _substitute_formula(condition, {variable: solution})
        branch_domain = _substitute_formula(sp.Ne(denominator, 0), {variable: solution})
        branches.append(sp.And(branch_condition, branch_domain))
    return sp.simplify(sp.Or(*branches))


def _try_complete_cad_image(
    expr: sp.Expr,
    condition: sp.Expr,
    variables: tuple[sp.Symbol, ...],
    value_symbol: sp.Symbol,
) -> sp.Expr | None:
    """Run the CAD/QE image formulation for polynomial/rational expressions."""

    relation, domain_constraint = _relation_for_function_graph(expr, value_symbol)
    image_formula = sp.And(condition, domain_constraint, relation)
    try:
        parsed = parse_formula(image_formula)
        result = qe_by_complete_cad(
            (value_symbol, *variables),
            tuple(("exists", variable) for variable in variables),
            parsed,
            free_variables=(value_symbol,),
        )
    except _EXPECTED_ERRORS:
        return None
    return sp.simplify(result.formula)


def _simplify_range_condition(formula: sp.Expr, value_symbol: sp.Symbol) -> sp.Expr:
    if formula is sp.true or formula is sp.false:
        return formula
    if isinstance(formula, SymAnd):
        return sp.And(*(_simplify_range_condition(arg, value_symbol) for arg in formula.args))
    if isinstance(formula, SymOr):
        return sp.Or(*(_simplify_range_condition(arg, value_symbol) for arg in formula.args))
    if getattr(formula, "is_Relational", False) and formula.free_symbols <= {value_symbol}:
        return formula
    return formula


def _intervals_from_range_set(
    range_set: sp.Set,
) -> tuple[tuple[sp.Expr, sp.Expr, bool, bool], ...] | None:
    """Return intervals as ``(left, right, left_closed, right_closed)``.

    The helper intentionally accepts only one-dimensional set forms that have
    direct range-bound semantics: intervals, finite point sets, and finite
    unions of those. More exotic sets leave metadata unknown without changing
    the primary range formula.
    """

    if range_set is sp.S.EmptySet:
        return ()
    if range_set is sp.S.Reals or range_set is sp.S.UniversalSet:
        return ((-sp.oo, sp.oo, False, False),)
    if isinstance(range_set, sp.Interval):
        return (
            (range_set.start, range_set.end, not range_set.left_open, not range_set.right_open),
        )
    if isinstance(range_set, sp.FiniteSet):
        return tuple(
            (item, item, True, True) for item in sorted(range_set, key=sp.default_sort_key)
        )
    if isinstance(range_set, sp.Union):
        pieces: list[tuple[sp.Expr, sp.Expr, bool, bool]] = []
        for arg in range_set.args:
            intervals = _intervals_from_range_set(arg)
            if intervals is None:
                return None
            pieces.extend(intervals)
        return tuple(sorted(pieces, key=lambda item: sp.default_sort_key(item[0])))
    return None


def _compare_bounds_for_sort(a: sp.Expr, b: sp.Expr) -> int:
    if a == b:
        return 0
    if a == -sp.oo or b == sp.oo:
        return -1
    if a == sp.oo or b == -sp.oo:
        return 1
    return _finite_compare(a, b)


def _range_metadata_from_formula(
    formula: sp.Expr,
    value_symbol: sp.Symbol,
) -> tuple[sp.Expr | None, sp.Expr | None, bool | None, bool | None, bool | None, int | None]:
    """Extract interval metadata from a one-variable range condition.

    The range formula remains authoritative. These fields are summaries for
    interval-like answers and are intentionally ``None`` when extraction is not
    reliable.
    """

    range_set = _domain_set_from_condition(formula, value_symbol)
    if range_set is None:
        try:
            range_set = formula.as_set()
        except _EXPECTED_ERRORS:
            return None, None, None, None, None, None
    intervals = _intervals_from_range_set(range_set)
    if intervals is None:
        return None, None, None, None, None, None
    interval_count = len(intervals)
    if interval_count == 0:
        return None, None, False, False, True, 0
    left = intervals[0][0]
    right = intervals[0][1]
    left_closed = intervals[0][2]
    right_closed = intervals[0][3]
    for item in intervals[1:]:
        if _compare_bounds_for_sort(item[0], left) < 0:
            left, left_closed = item[0], item[2]
        if _compare_bounds_for_sort(right, item[1]) < 0:
            right, right_closed = item[1], item[3]
    lower_attained = bool(left_closed) if left not in (-sp.oo, sp.oo) else False
    upper_attained = bool(right_closed) if right not in (-sp.oo, sp.oo) else False
    is_interval = interval_count == 1
    return left, right, lower_attained, upper_attained, is_interval, interval_count


def _range_result_from_formula(
    expression: sp.Expr,
    formula: sp.Expr,
    value_symbol: sp.Symbol,
    variables: tuple[sp.Symbol, ...],
    method: str,
    diagnostics: Mapping[str, object],
) -> FunctionRangeResult:
    simplified_formula = _simplify_range_condition(formula, value_symbol)
    lo, hi, lo_attained, hi_attained, is_interval, interval_count = _range_metadata_from_formula(
        simplified_formula,
        value_symbol,
    )
    return FunctionRangeResult(
        expression,
        simplified_formula,
        value_symbol,
        variables,
        lo,
        hi,
        lo_attained,
        hi_attained,
        (),
        (),
        method,
        {
            **dict(diagnostics),
            "metadata_source": "range_formula_as_set" if interval_count is not None else "unknown",
        },
        is_interval,
        interval_count,
    )


def _range_via_optimization_bounds(
    expr: sp.Expr,
    condition: sp.Expr,
    vars_: tuple[sp.Symbol, ...],
    val_sym: sp.Symbol,
    domain: str,
) -> FunctionRangeResult:
    minimum = semialgebraic_minimize(expr, condition, vars_, domain=domain, return_result=True)
    maximum = semialgebraic_maximize(expr, condition, vars_, domain=domain, return_result=True)
    assert isinstance(minimum, OptimizationResult)
    assert isinstance(maximum, OptimizationResult)
    formula = sp.simplify(
        _range_formula(val_sym, minimum.value, maximum.value, minimum.attained, maximum.attained)
    )
    _, _, _, _, is_interval, interval_count = _range_metadata_from_formula(formula, val_sym)
    return FunctionRangeResult(
        expr,
        formula,
        val_sym,
        vars_,
        minimum.value,
        maximum.value,
        minimum.attained,
        maximum.attained,
        minimum.points,
        maximum.points,
        "optimization_bounds",
        {"constraints": sp.sstr(condition), "metadata_source": "optimization_bounds"},
        is_interval,
        interval_count,
    )


@with_computation_context
def function_range(
    expression: sp.Expr,
    constraints: FormulaLike | Iterable[FormulaLike] | None = None,
    variables: Sequence[sp.Symbol | str] | None = None,
    *,
    value_symbol: sp.Symbol | str | None = None,
    domain: str = "reals",
    method: str = "qe",
    return_result: bool = False,
    parameters: Sequence[sp.Symbol | str] | None = None,
    return_stratified: bool = False,
) -> sp.Expr | FunctionRangeResult | object:
    """Return a quantifier-free formula describing a real function range.

    The preferred direct backend uses the semialgebraic image formulation
    ``exists variables. constraints and value_symbol == expression``. It first
    applies guarded exact graph-elimination shortcuts for common univariate
    images, then tries complete CAD/QE, and may use optimization bounds for
    ``method='auto'`` or ``method='bounds'``.
    """

    if domain.lower() not in {"real", "reals", "r", "rr"}:
        raise NotImplementedError("function_range currently supports only the real domain")
    requested_method = method.lower()
    if requested_method not in {"qe", "cad", "auto", "bounds", "optimization"}:
        raise ValueError("method must be 'qe', 'cad', 'auto', 'bounds', or 'optimization'")
    expr = sp.sympify(expression)
    condition = _normalize_formula(constraints)
    vars_ = _normalize_variables(variables, sp.Tuple(condition, expr))
    val_sym = resolve_symbol(value_symbol or "t")
    if return_stratified:
        params = _normalize_parameters_for_problem(parameters or (), expr, condition)
        if not params:
            raise ValueError("return_stratified=True requires at least one parameter")
        vars_ = tuple(var for var in vars_ if var not in set(params))
        return _stratified_function_range(expr, condition, vars_, params, val_sym)

    if not vars_:
        formula = sp.Eq(val_sym, expr)
        result = FunctionRangeResult(
            expr,
            formula,
            val_sym,
            vars_,
            expr,
            expr,
            True,
            True,
            (),
            (),
            "constant",
            {"metadata_source": "constant"},
            True,
            1,
        )
        return result if return_result else result.formula

    if requested_method in {"bounds", "optimization"}:
        result = _range_via_optimization_bounds(expr, condition, vars_, val_sym, domain)
        return result if return_result else result.formula

    direct_special = _try_direct_special_image(expr, condition, vars_, val_sym)
    if direct_special is not None:
        direct_method = "direct_univariate_image"
        try:
            if len(vars_) == 1 and sp.Poly(sp.expand(expr), vars_[0]).degree() >= 2:
                direct_method = "optimization_bounds_direct_univariate_image"
        except _EXPECTED_ERRORS:
            pass
        result = _range_result_from_formula(
            expr,
            direct_special,
            val_sym,
            vars_,
            direct_method,
            {"constraints": sp.sstr(condition), "requested_method": requested_method},
        )
        return result if return_result else result.formula

    semialgebraic_graph = _try_semialgebraic_graph_image(expr, condition, vars_, val_sym)
    if semialgebraic_graph is not None:
        result = _range_result_from_formula(
            expr,
            semialgebraic_graph,
            val_sym,
            vars_,
            "qe_image_semialgebraic_graph",
            {"constraints": sp.sstr(condition), "requested_method": requested_method},
        )
        return result if return_result else result.formula

    direct = _try_solved_graph_image(expr, condition, vars_, val_sym)
    if direct is not None:
        result = _range_result_from_formula(
            expr,
            direct,
            val_sym,
            vars_,
            "qe_image_solved_graph",
            {"constraints": sp.sstr(condition), "requested_method": requested_method},
        )
        return result if return_result else result.formula

    cad_formula = (
        _try_complete_cad_image(expr, condition, vars_, val_sym)
        if requested_method == "cad"
        else None
    )
    if cad_formula is not None:
        result = _range_result_from_formula(
            expr,
            cad_formula,
            val_sym,
            vars_,
            "qe_image_complete_cad",
            {"constraints": sp.sstr(condition), "requested_method": requested_method},
        )
        return result if return_result else result.formula

    if requested_method == "cad":
        raise NotImplementedError(
            "CAD/QE-image range computation failed for this expression/domain"
        )

    result = _range_via_optimization_bounds(expr, condition, vars_, val_sym, domain)
    result = FunctionRangeResult(
        result.expression,
        result.formula,
        result.value_symbol,
        result.variables,
        result.infimum,
        result.supremum,
        result.minimum_attained,
        result.maximum_attained,
        result.minimizers,
        result.maximizers,
        "optimization_bounds_after_qe_image_fallback",
        {**dict(result.diagnostics), "requested_method": requested_method},
        result.is_interval,
        result.interval_count,
    )
    return result if return_result else result.formula


__all__ = [
    "FunctionRangeResult",
    "OptimizationResult",
    "ParametricOptimizationResult",
    "ParametricFunctionRangeResult",
    "OptimizationCertificationPolicy",
    "polynomial_locus_dimension",
    "function_range",
    "semialgebraic_maximize",
    "semialgebraic_minimize",
]
