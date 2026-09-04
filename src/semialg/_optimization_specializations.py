from __future__ import annotations

from collections.abc import Mapping
from itertools import combinations

import sympy as sp
from sympy.core.sympify import SympifyError
from sympy.polys.polyerrors import GeneratorsNeeded, PolynomialError

from ._critical_loci import LocusKey
from .exact_arithmetic import compare_exact_reals
from .interval_decomposition import (
    finite_real_roots as _finite_real_roots,
)
from .normalization import require_conjunction
from .optimization_results import (
    OptimizationCertificationPolicy,
    OptimizationResult,
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


from .optimization import (  # noqa: E402
    _best_candidates,
    _Candidate,
    _constraint_data,
    _is_feasible,
    _lift_reduced_points,
    _real_point,
    _reduce_linear_equalities,
)


def _call_optimize_conjunction(*args, **kwargs):
    from .optimization import _optimize_conjunction

    return _optimize_conjunction(*args, **kwargs)


def _affine_polynomial(expr: sp.Expr, variables: tuple[sp.Symbol, ...]) -> bool:
    try:
        return sp.Poly(sp.expand(expr), *variables).total_degree() <= 1
    except _EXPECTED_ERRORS:
        return False


def _closed_linear_polytope(condition: sp.Expr, variables: tuple[sp.Symbol, ...]) -> bool:
    """Return whether a conjunction is made only of closed affine relations."""

    for atom in require_conjunction(
        condition, message="linear-polytope specialization expects a conjunction"
    ):
        if atom in (sp.true, True):
            continue
        if isinstance(atom, (sp.StrictLessThan, sp.StrictGreaterThan, sp.Unequality)):
            return False
        if not getattr(atom, "is_Relational", False):
            return False
        residual, _ = _relation_parts(atom)
        if not _affine_polynomial(residual, variables):
            return False
    return True


def _linear_polytope_specialization(
    objective: sp.Expr,
    condition: sp.Expr,
    variables: tuple[sp.Symbol, ...],
    *,
    kind: str,
) -> OptimizationResult | None:
    """Optimize an affine objective over a bounded closed polytope by vertices.

    Boundedness is certified by the package's affine projection-bounds path.
    Once compactness is established, the fundamental theorem of linear
    programming certifies that an affine objective reaches an optimum at an
    extreme point.  Candidate extreme points are obtained by exact active
    affine systems; feasibility is checked against the original formula.
    """

    if not variables or not _affine_polynomial(objective, variables):
        return None
    if not _closed_linear_polytope(condition, variables):
        return None
    try:
        from .reasoning_regions import region_bounded

        if not region_bounded(condition, variables):
            return None
    except _EXPECTED_ERRORS:
        return None

    equalities, inequalities, _ = _constraint_data(condition, variables)
    n = len(variables)
    if len(equalities) > n:
        # Redundant equalities are handled reliably by the general path.
        return None
    active_needed = max(0, n - len(equalities))
    if active_needed > len(inequalities):
        return None
    points: list[dict[sp.Symbol, sp.Expr]] = []
    for active in combinations(inequalities, active_needed):
        equations = (*equalities, *active)
        try:
            solution = sp.linsolve(equations, variables)
        except _EXPECTED_ERRORS:
            continue
        if solution is sp.EmptySet:
            continue
        for values in solution:
            if any(value.free_symbols & set(variables) for value in values):
                continue
            point = {var: sp.simplify(value) for var, value in zip(variables, values, strict=True)}
            if _real_point(point) and _is_feasible(condition, point):
                if point_key(point) not in {point_key(p) for p in points}:
                    points.append(point)
    if not points:
        # A nonempty zero-dimensional affine region can be caught by the
        # general equality-reduction path; do not guess here.
        return None
    candidates = tuple(
        _Candidate(sp.simplify(objective.subs(point)), point, True) for point in points
    )
    value, best_points, attained = _best_candidates(candidates, kind=kind)
    return OptimizationResult(
        objective,
        variables,
        value,
        best_points,
        attained,
        kind,
        "exact_linear_polytope_vertices",
        {
            "candidate_count": len(points),
            "constraint_formula": condition,
            "global_certificate": "bounded_closed_polytope_extreme_point_theorem",
        },
        True,
    )


def _separable_groups(
    objective: sp.Expr,
    condition: sp.Expr,
    variables: tuple[sp.Symbol, ...],
) -> tuple[sp.Expr, dict[sp.Symbol, sp.Expr], dict[sp.Symbol, sp.Expr]] | None:
    """Return constant/objective/constraint pieces for a Cartesian product."""

    varset = set(variables)
    objective_parts = {var: sp.Integer(0) for var in variables}
    constant = sp.Integer(0)
    for term in sp.Add.make_args(sp.expand(objective)):
        support = term.free_symbols & varset
        if not support:
            constant += term
        elif len(support) == 1:
            objective_parts[next(iter(support))] += term
        else:
            return None
    constraint_parts: dict[sp.Symbol, list[sp.Expr]] = {var: [] for var in variables}
    for atom in require_conjunction(
        condition, message="separable optimization expects a conjunction"
    ):
        if atom in (sp.true, True):
            continue
        support = atom.free_symbols & varset
        if not support:
            try:
                if not bool(sp.simplify(atom)):
                    return None
            except (TypeError, ValueError):
                return None
        elif len(support) == 1:
            constraint_parts[next(iter(support))].append(atom)
        else:
            return None
    return (
        sp.simplify(constant),
        {var: sp.simplify(expr) for var, expr in objective_parts.items()},
        {var: sp.And(*atoms) if atoms else sp.true for var, atoms in constraint_parts.items()},
    )


def _closed_interval_piece(
    objective: sp.Expr,
    condition: sp.Expr,
    variable: sp.Symbol,
    *,
    kind: str,
) -> OptimizationResult | None:
    """Exact polynomial optimization on one explicit finite closed interval."""

    for atom in require_conjunction(
        condition, message="closed interval specialization expects a conjunction"
    ):
        if isinstance(atom, (sp.StrictLessThan, sp.StrictGreaterThan, sp.Unequality)):
            return None
    try:
        from .implicit_geometry import extract_symbolic_box_bounds

        bounds = extract_symbolic_box_bounds(condition, (variable,))
        _, lower, upper = bounds.limits[0]
    except _EXPECTED_ERRORS:
        return None
    if lower in (-sp.oo, sp.oo) or upper in (-sp.oo, sp.oo):
        return None
    if (lower.free_symbols | upper.free_symbols) & {variable}:
        return None
    try:
        poly = sp.Poly(sp.expand(objective), variable)
    except _EXPECTED_ERRORS:
        return None
    points = [sp.simplify(lower), sp.simplify(upper)]
    derivative = sp.diff(poly.as_expr(), variable)
    if derivative != 0:
        try:
            roots = _finite_real_roots(derivative, variable)
        except _EXPECTED_ERRORS:
            return None
        for root in roots:
            try:
                if compare_exact_reals(lower, root) <= 0 and compare_exact_reals(root, upper) <= 0:
                    points.append(root)
            except _EXPECTED_ERRORS:
                return None
    candidates = tuple(
        _Candidate(
            sp.simplify(objective.subs(variable, point)),
            {variable: point},
            True,
        )
        for point in dict.fromkeys(points)
        if _is_feasible(condition, {variable: point})
    )
    if not candidates:
        return None
    value, best_points, attained = _best_candidates(candidates, kind=kind)
    return OptimizationResult(
        objective,
        (variable,),
        value,
        best_points,
        attained,
        kind,
        "exact_closed_interval_calculus",
        {
            "candidate_count": len(candidates),
            "constraint_formula": condition,
            "global_certificate": "polynomial_extrema_on_compact_interval",
        },
        True,
    )


def _separable_specialization(
    objective: sp.Expr,
    condition: sp.Expr,
    variables: tuple[sp.Symbol, ...],
    *,
    kind: str,
    policy: OptimizationCertificationPolicy,
    recursion_depth: int,
    visited_loci: frozenset[LocusKey],
) -> OptimizationResult | None:
    """Optimize Cartesian-product problems as independent exact 1-D problems."""

    if len(variables) < 2:
        return None
    split = _separable_groups(objective, condition, variables)
    if split is None:
        return None
    constant, objective_parts, constraint_parts = split
    subresults: list[OptimizationResult] = []
    for var in variables:
        sub = _closed_interval_piece(objective_parts[var], constraint_parts[var], var, kind=kind)
        if sub is None:
            sub = _call_optimize_conjunction(
                objective_parts[var],
                constraint_parts[var],
                (var,),
                kind=kind,
                policy=policy,
                recursion_depth=recursion_depth + 1,
                visited_loci=visited_loci,
                allow_equality_reduction=True,
                allow_specializations=False,
            )
        subresults.append(sub)
    value = sp.simplify(constant + sum((item.value for item in subresults), sp.Integer(0)))
    attained = all(item.attained for item in subresults)
    points: tuple[Mapping[sp.Symbol, sp.Expr], ...] = ()
    if attained and all(item.points for item in subresults):
        # Cartesian product of minimizer sets can explode. Preserve one exact
        # representative point while reporting the component counts.
        point = {var: subresults[i].points[0][var] for i, var in enumerate(variables)}
        points = (point,)
    return OptimizationResult(
        objective,
        variables,
        value,
        points,
        attained,
        kind,
        "separable_exact_optimization",
        {
            "constraint_formula": condition,
            "component_methods": tuple(item.method for item in subresults),
            "component_point_counts": tuple(len(item.points) for item in subresults),
        },
        all(item.certified for item in subresults),
    )


def _try_optimization_specializations(
    objective: sp.Expr,
    condition: sp.Expr,
    variables: tuple[sp.Symbol, ...],
    *,
    kind: str,
    policy: OptimizationCertificationPolicy,
    recursion_depth: int,
    visited_loci: frozenset[LocusKey],
) -> OptimizationResult | None:
    """Run cheap exact specializations before the general KKT path."""

    linear = _linear_polytope_specialization(objective, condition, variables, kind=kind)
    if linear is not None:
        return linear
    return _separable_specialization(
        objective,
        condition,
        variables,
        kind=kind,
        policy=policy,
        recursion_depth=recursion_depth,
        visited_loci=visited_loci,
    )


def _try_equality_reduced_optimization(
    objective: sp.Expr,
    condition: sp.Expr,
    variables: tuple[sp.Symbol, ...],
    *,
    kind: str,
    policy: OptimizationCertificationPolicy,
    recursion_depth: int,
    visited_loci: frozenset[LocusKey],
) -> OptimizationResult | None:
    """Eliminate linear equalities and lift an exact reduced optimum."""

    if not variables:
        return None
    reduced_obj, reduced_cond, reduced_vars, substitutions = _reduce_linear_equalities(
        objective, condition, variables
    )
    if len(reduced_vars) >= len(variables):
        return None
    reduced_result = _call_optimize_conjunction(
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
