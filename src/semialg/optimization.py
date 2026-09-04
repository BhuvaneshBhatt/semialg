from __future__ import annotations

from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass
from typing import Literal

import sympy as sp
from sympy.core.sympify import SympifyError
from sympy.logic.boolalg import BooleanFalse, BooleanTrue
from sympy.polys.polyerrors import GeneratorsNeeded, PolynomialError

from ._critical_loci import (
    LocusKey,
    ProjectedCriticalLocus,
    iter_positive_dimensional_critical_loci,
)
from ._linear_relations import safe_linear_solution
from .context import with_computation_context
from .exact_arithmetic import compare_exact_reals
from .formulas.boolean import bounded_dnf_branches, is_false_expr
from .interval_decomposition import (
    finite_real_roots as _finite_real_roots,
)
from .interval_decomposition import (
    one_dimensional_intervals as _one_dimensional_intervals,
)
from .interval_decomposition import (
    relational_polynomials as _relational_polynomials,
)
from .normalization import normalize_constraints, normalize_problem_variables, require_conjunction
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
from .structural_keys import point_key

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


def satisfies_formula(*args, **kwargs):
    """Import instance checking only for optimization candidate validation."""

    from .instances.real_fallbacks import satisfies_formula as impl

    return impl(*args, **kwargs)


@dataclass(frozen=True)
class _Candidate:
    value: sp.Expr
    point: Mapping[sp.Symbol, sp.Expr] | None
    attained: bool


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
    seen: set[tuple[tuple[sp.Symbol, sp.Expr], ...]] = set()
    for item in tied:
        if item.attained and item.point is not None:
            key = point_key(item.point)
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
    cuts: set[sp.Expr] = set()
    cut_values: list[sp.Expr] = []
    for poly in _relational_polynomials(condition):
        if poly.free_symbols <= {variable}:
            for root in _finite_real_roots(poly, variable):
                key = root
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
    for atom in require_conjunction(
        condition, message="optimization candidate enumeration supports conjunctions"
    ):
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
        key = (op if op in {"==", "!="} else "ineq", sp.factor(residual))
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
        for atom in require_conjunction(
            condition, message="optimization candidate enumeration supports conjunctions"
        ):
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

    expanded = tuple(sp.expand(eq) for eq in equations)
    equations = tuple(eq for eq in expanded if eq != 0)
    if not equations or not solve_variables:
        return ()
    points: list[dict[sp.Symbol, sp.Expr]] = []
    try:
        from .solve.zero_dimensional import is_zero_dimensional, solve_zero_dimensional_system

        if is_zero_dimensional(equations, solve_variables):
            result = solve_zero_dimensional_system(
                equations,
                vars=tuple(solve_variables),
                backend="rur",
                real=True,
                return_result=True,
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
                replacement = safe_linear_solution(eq, var)
                if replacement is None:
                    continue
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
    seen: set[tuple[tuple[sp.Symbol, sp.Expr], ...]] = set()
    subsets = _pruned_active_subsets(equalities, inequalities, variables, condition)
    for subset in subsets:
        active = tuple(equalities) + tuple(subset)
        equations, multipliers = _kkt_system(objective, variables, active)
        solve_vars = (*variables, *multipliers)
        dimension = polynomial_locus_dimension(equations, solve_vars)
        if dimension == 0:
            for point in _solve_exact_equations(equations, solve_vars, variables, condition):
                key = point_key(point)
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
                    key = point_key(point)
                    if key not in seen:
                        points.append(point)
                        seen.add(key)

        if len(active) >= len(variables) and polynomial_locus_dimension(active, variables) == 0:
            for point in _solve_exact_equations(active, variables, variables, condition):
                key = point_key(point)
                if key not in seen:
                    points.append(point)
                    seen.add(key)
    return tuple(points)


def _attained_point_on_value_locus(
    objective: sp.Expr,
    condition: sp.Expr,
    variables: tuple[sp.Symbol, ...],
    value: sp.Expr,
) -> Mapping[sp.Symbol, sp.Expr] | None:
    """Return one exact optimizer witness after an image-range endpoint is attained."""

    try:
        from .solve.find_instance import find_instance

        point = find_instance(
            sp.And(condition, sp.Eq(objective, value), evaluate=False),
            variables,
            strategy="auto",
            exact=True,
        )
    except _EXPECTED_ERRORS:
        return None
    return point if isinstance(point, Mapping) else None


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
    visited_loci: frozenset[LocusKey],
) -> tuple[_Candidate, ...]:
    """Optimize recursively over positive-dimensional critical/singular loci.

    Projected KKT loci cover ordinary positive-dimensional critical sets.  Active
    constraint Jacobian-rank-deficiency loci are included separately because a
    singular positive-dimensional boundary can carry varying objective values even
    when no KKT multiplier solves the ordinary stationarity equations.
    """

    if recursion_depth >= policy.recursion_limit:
        return ()
    out: list[_Candidate] = []
    seen_loci: set[LocusKey] = set()

    def consume_locus(locus: ProjectedCriticalLocus) -> None:
        """Optimize over one projected critical locus and reuse its structural identity in recursion."""
        if locus.key in visited_loci or locus.key in seen_loci:
            return
        seen_loci.add(locus.key)
        locus_condition = sp.And(
            condition, *(sp.Eq(eq, 0) for eq in locus.equations), evaluate=False
        )
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
                    visited_loci=visited_loci | frozenset({locus.key}),
                    allow_equality_reduction=True,
                )
                lifted = _lift_reduced_points(result.points, substitutions, variables)
                out.append(_Candidate(result.value, lifted[0] if lifted else None, result.attained))
                return

            cert = _certify_optimum_by_range(
                objective, locus_condition, variables, kind=kind, policy=policy
            )
            if cert is None:
                return
            value, attained = cert
            point = (
                _attained_point_on_value_locus(objective, locus_condition, variables, value)
                if attained
                else None
            )
            out.append(_Candidate(value, point, attained))
        except (NotImplementedError, ValueError, TypeError, ArithmeticError, PolynomialError):
            return

    for locus in iter_positive_dimensional_critical_loci(
        objective, condition, variables, equalities, inequalities, include_singular=True
    ):
        consume_locus(locus)
    return tuple(out)


def _multivariate_exact_candidates(
    objective: sp.Expr,
    condition: sp.Expr,
    variables: tuple[sp.Symbol, ...],
    *,
    kind: str,
    policy: OptimizationCertificationPolicy,
    recursion_depth: int = 0,
    visited_loci: frozenset[LocusKey] = frozenset(),
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
    visited_loci: frozenset[LocusKey] = frozenset(),
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


from ._optimization_certification import (
    _certify_candidate_by_qe,
    _certify_optimum_by_range,
    clear_optimization_range_cache,
    optimization_range_cache_info,
)
from ._optimization_specializations import (
    _try_equality_reduced_optimization,
    _try_optimization_specializations,
)


def _finalize_polynomial_optimization(
    objective: sp.Expr,
    condition: sp.Expr,
    variables: tuple[sp.Symbol, ...],
    candidates: Sequence[tuple[sp.Expr, Mapping[sp.Symbol, sp.Expr], bool]],
    *,
    kind: str,
    policy: OptimizationCertificationPolicy,
) -> OptimizationResult:
    """Choose candidates and apply the two global exact certification paths."""

    candidate_value = None
    candidate_points: tuple[Mapping[sp.Symbol, sp.Expr], ...] = ()
    candidate_attained = False
    if candidates:
        candidate_value, candidate_points, candidate_attained = _best_candidates(
            candidates, kind=kind
        )

    common_diagnostics = {
        "candidate_count": len(candidates),
        "constraints": sp.sstr(condition),
        "constraint_formula": condition,
    }
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
            {**common_diagnostics, "global_certificate": "complete_cad_no_better_point"},
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
            {**common_diagnostics, "global_certificate": "complete_cad_function_range"},
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
        {**common_diagnostics, "global_certificate": "candidate_exhaustion_without_cad_range"},
        False,
    )


def _optimize_conjunction(
    objective: sp.Expr,
    condition: sp.Expr,
    variables: tuple[sp.Symbol, ...],
    *,
    kind: str,
    policy: OptimizationCertificationPolicy,
    recursion_depth: int = 0,
    visited_loci: frozenset[LocusKey] = frozenset(),
    allow_equality_reduction: bool = True,
    allow_specializations: bool = True,
) -> OptimizationResult:
    """Optimize one conjunction with exact specializations, KKT, and certificates."""

    if is_false_expr(condition):
        raise ValueError("optimization domain is empty")
    if allow_specializations:
        specialized = _try_optimization_specializations(
            objective,
            condition,
            variables,
            kind=kind,
            policy=policy,
            recursion_depth=recursion_depth,
            visited_loci=visited_loci,
        )
        if specialized is not None:
            return specialized

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

    if allow_equality_reduction:
        reduced = _try_equality_reduced_optimization(
            objective,
            condition,
            variables,
            kind=kind,
            policy=policy,
            recursion_depth=recursion_depth,
            visited_loci=visited_loci,
        )
        if reduced is not None:
            return reduced

    candidates = _optimization_candidates(
        objective,
        condition,
        variables,
        kind=kind,
        policy=policy,
        recursion_depth=recursion_depth,
        visited_loci=visited_loci,
    )
    return _finalize_polynomial_optimization(
        objective, condition, variables, candidates, kind=kind, policy=policy
    )


from ._optimization_parametric import _normalize_parameters_for_problem, _stratified_optimization


def _combine_branch_results(
    objective: sp.Expr,
    variables: tuple[sp.Symbol, ...],
    branch_results: Sequence[OptimizationResult],
    *,
    kind: str,
) -> OptimizationResult:
    """Combine exact branch optima and preserve the weakest certification status."""

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
    certified = all(item.certified for item in branch_results)
    return OptimizationResult(
        objective,
        variables,
        best.value,
        tuple(points),
        any(result.attained for result in tied),
        kind,
        best.method if len(branch_results) == 1 else "exact_branchwise_optimization",
        {
            **dict(best.diagnostics),
            "branch_count": len(branch_results),
            "all_branches_certified": certified,
        },
        certified,
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
    max_boolean_branches: int = 32,
    parameters: Sequence[sp.Symbol | str] | None = None,
    return_stratified: bool = False,
    eliminate_quantifiers: bool = False,
    branch_solver=None,
) -> OptimizationResult | sp.Expr | object:
    """Shared exact implementation for minimization and maximization APIs.

    Boolean domains are decomposed into boundedly many conjunctions. Each
    branch is optimized independently and exact branch bounds are compared
    before the combined result and certificate status are assembled.
    """

    if max_boolean_branches < 1:
        raise ValueError("max_boolean_branches must be positive")
    if domain.lower() not in {"real", "reals", "r", "rr"}:
        raise NotImplementedError("semialgebraic optimization supports only the real domain")
    policy = OptimizationCertificationPolicy(certification, range_cost_limit, recursion_limit)
    obj = sp.sympify(objective)
    condition = normalize_constraints(constraints)
    vars_ = normalize_problem_variables(variables, sp.Tuple(condition, obj))
    if not vars_:
        value = sp.simplify(obj)
        result = OptimizationResult(obj, vars_, value, (), True, kind, "constant", {}, True)
        return result if return_result else result.value

    expansion = bounded_dnf_branches(condition, max_branches=max_boolean_branches)
    if not expansion.complete:
        raise NotImplementedError(
            f"optimization Boolean expansion exceeded max_boolean_branches={max_boolean_branches}"
        )
    branch_results: list[OptimizationResult] = []
    solve_branch = branch_solver or _optimize_conjunction
    for branch in expansion.branches:
        branch_condition = sp.And(
            *[piece for piece in branch if piece not in (True, sp.true)], evaluate=False
        )
        if any(piece in (False, sp.false) for piece in branch):
            continue
        try:
            branch_results.append(
                solve_branch(obj, branch_condition, vars_, kind=kind, policy=policy)
            )
        except ValueError:
            continue
    result = _combine_branch_results(obj, vars_, branch_results, kind=kind)
    return result if return_result else result.value


@with_computation_context
def semialgebraic_minimize(
    objective: sp.Expr,
    constraints: FormulaLike | Iterable[FormulaLike] | None = None,
    variables: Sequence[sp.Symbol | str] | None = None,
    *,
    domain: str = "reals",
    return_result: bool = False,
    certification: Literal["auto", "complete", "candidate"] = "auto",
    range_cost_limit: int = 2500,
    recursion_limit: int = 4,
    max_boolean_branches: int = 32,
    parameters: Sequence[sp.Symbol | str] | None = None,
    return_stratified: bool = False,
    eliminate_quantifiers: bool = False,
) -> list[object] | OptimizationResult | object:
    """Return ``[minimum_or_infimum, optimizer_points]`` by default.

    Set ``return_result=True`` for an :class:`OptimizationResult` with
    attainment, certification, method, and diagnostic metadata.

    Multivariate polynomial problems use exact active-set/KKT enumeration,
    singular active-locus solving, RUR-backed zero-dimensional solving, and a
    complete-CAD image certificate when available.
    """

    if eliminate_quantifiers and not return_stratified:
        raise ValueError("eliminate_quantifiers=True requires return_stratified=True")
    if return_stratified:
        obj = sp.sympify(objective)
        condition = normalize_constraints(constraints)
        params = _normalize_parameters_for_problem(parameters or (), obj, condition)
        if not params:
            raise ValueError("return_stratified=True requires at least one parameter")
        vars_ = normalize_problem_variables(variables, sp.Tuple(condition, obj))
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
            max_boolean_branches=max_boolean_branches,
            eliminate_quantifiers=eliminate_quantifiers,
        )
    result = _optimize(
        objective,
        constraints,
        variables,
        kind="min",
        domain=domain,
        return_result=True,
        certification=certification,
        range_cost_limit=range_cost_limit,
        recursion_limit=recursion_limit,
        max_boolean_branches=max_boolean_branches,
    )
    if return_result:
        return result
    if not isinstance(result, OptimizationResult):
        raise TypeError("optimization backend returned an unexpected result type")
    return [result.value, list(result.points)]


@with_computation_context
def semialgebraic_maximize(
    objective: sp.Expr,
    constraints: FormulaLike | Iterable[FormulaLike] | None = None,
    variables: Sequence[sp.Symbol | str] | None = None,
    *,
    domain: str = "reals",
    return_result: bool = False,
    certification: Literal["auto", "complete", "candidate"] = "auto",
    range_cost_limit: int = 2500,
    recursion_limit: int = 4,
    max_boolean_branches: int = 32,
    parameters: Sequence[sp.Symbol | str] | None = None,
    return_stratified: bool = False,
    eliminate_quantifiers: bool = False,
) -> list[object] | OptimizationResult | object:
    """Return ``[maximum_or_supremum, optimizer_points]`` by default.

    Set ``return_result=True`` for an :class:`OptimizationResult` with
    attainment, certification, method, and diagnostic metadata.
    """

    if eliminate_quantifiers and not return_stratified:
        raise ValueError("eliminate_quantifiers=True requires return_stratified=True")
    if return_stratified:
        obj = sp.sympify(objective)
        condition = normalize_constraints(constraints)
        params = _normalize_parameters_for_problem(parameters or (), obj, condition)
        if not params:
            raise ValueError("return_stratified=True requires at least one parameter")
        vars_ = normalize_problem_variables(variables, sp.Tuple(condition, obj))
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
            max_boolean_branches=max_boolean_branches,
            eliminate_quantifiers=eliminate_quantifiers,
        )
    result = _optimize(
        objective,
        constraints,
        variables,
        kind="max",
        domain=domain,
        return_result=True,
        certification=certification,
        range_cost_limit=range_cost_limit,
        recursion_limit=recursion_limit,
        max_boolean_branches=max_boolean_branches,
    )
    if return_result:
        return result
    if not isinstance(result, OptimizationResult):
        raise TypeError("optimization backend returned an unexpected result type")
    return [result.value, list(result.points)]


from ._optimization_range import function_range  # noqa: F401
from ._range_special_cases import _graph_formula_for_expression  # noqa: F401

__all__ = [
    "FunctionRangeResult",
    "clear_optimization_range_cache",
    "OptimizationResult",
    "optimization_range_cache_info",
    "ParametricOptimizationResult",
    "ParametricFunctionRangeResult",
    "OptimizationCertificationPolicy",
    "polynomial_locus_dimension",
    "function_range",
    "semialgebraic_maximize",
    "semialgebraic_minimize",
]
