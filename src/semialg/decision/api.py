from __future__ import annotations

from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass, field

import sympy as sp
from sympy.logic.boolalg import Boolean
from sympy.polys.polyerrors import CoercionFailed, PolynomialError

from ..algebraic.equality_ideal import EqualityIdealContext
from ..algebraic.rational_univariate import solve_formula_with_rur
from ..decision_diagnostics import solution_capability_diagnostics
from ..formula import ParsedPrenexFormula, parse_formula
from ..formulas.boolean import is_false_expr, is_true_expr
from ..incidence import decompose_conjunctive_formula
from ..inequality_reduction import reduce_conjunctive_inequalities
from ..normalization import conjuncts, normalize_formula
from ..qe import qe_by_complete_cad
from ..simplify.boolean import simplify_boolean
from ..solve import reduce_formula
from ..solve.planner import (
    affine_presolve,
    exact_linear_feasibility,
    groebner_reduce_conjunction,
    profile_semialgebraic_system,
    reconstruct_affine_assignment,
    reconstruct_affine_solution_formula,
)
from ..structural_keys import symbol_identity_key
from ._inputs import (
    as_real_symbol as _as_real_symbol,
)
from ._metadata import (
    collect_solution_metadata as _collect_solution_metadata,
)
from ._metadata import (
    components_formula as _components_formula,
)
from ._metadata import (
    one_dim_components as _one_dim_components,
)
from .solution import (
    EquivalenceResult,
    ImplicationResult,
    IntervalComponent,
    SatisfiabilityResult,
    SemialgebraicSolution,
    TautologyResult,
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


def _safe_simplify_expr(expr: sp.Expr) -> sp.Expr:
    """Simplify algebraic expressions without sending Boolean formulas to radsimp."""

    if isinstance(expr, Boolean):
        return simplify_boolean(expr)
    return sp.simplify(expr)


def _merge_variables(
    first: Sequence[sp.Symbol],
    second: Sequence[sp.Symbol],
    *formulas: sp.Expr,
) -> tuple[sp.Symbol, ...]:
    """Merge explicit variables and free symbols while preserving order."""

    out: list[sp.Symbol] = []
    seen: set[sp.Symbol] = set()
    for source in (first, second):
        for sym in source:
            if sym not in seen:
                out.append(sym)
                seen.add(sym)
    for formula in formulas:
        for sym in sorted(getattr(formula, "free_symbols", set()), key=symbol_identity_key):
            if sym not in seen:
                out.append(sym)
                seen.add(sym)
    return tuple(out)


def _fast_parameter_conditions(
    expr: sp.Expr,
    variables: tuple[sp.Symbol, ...],
    parameters: tuple[sp.Symbol, ...],
) -> sp.Expr | None:
    """Fast parameter conditions for common one-variable polynomial atoms."""

    if len(variables) != 1 or not isinstance(expr, sp.core.relational.Relational):
        return None
    var = variables[0]
    try:
        poly = sp.Poly(sp.expand(expr.lhs - expr.rhs), var)
    except _RECOVERABLE_ERRORS:
        return None
    degree = poly.degree()
    if isinstance(expr, sp.Equality):
        if degree == 0:
            return sp.Eq(poly.as_expr(), 0)
        if degree == 1:
            return sp.Ne(poly.LC(), 0)
        if degree == 2:
            a2, a1, a0 = poly.all_coeffs()
            discriminant = sp.expand(a1**2 - 4 * a2 * a0)
            quadratic = sp.And(sp.Ne(a2, 0), discriminant >= 0)
            linear = sp.And(sp.Eq(a2, 0), sp.Ne(a1, 0))
            constant_zero = sp.And(sp.Eq(a2, 0), sp.Eq(a1, 0), sp.Eq(a0, 0))
            return _safe_simplify_expr(sp.Or(quadratic, linear, constant_zero))
    return None


def _fast_solution_formula(
    expr: sp.Expr, variables: tuple[sp.Symbol, ...]
) -> tuple[sp.Expr, str, bool | None]:
    """Fast conservative formula/satisfiability path for common solves."""

    if is_true_expr(expr):
        return sp.true, "trivial", True
    if is_false_expr(expr):
        return sp.false, "trivial", False
    if len(variables) == 1:
        components = _one_dim_components(expr, variables[0])
        if components is not None:
            reduced = _components_formula(components)
            return reduced, "one_dimensional_components", bool(components)
        reduced = reduce_conjunctive_inequalities(expr, variables[0])
        if reduced is not None:
            return (
                reduced,
                "sympy_reduce_inequalities",
                reduced is not sp.false and reduced != sp.false,
            )
    if len(variables) == 2:
        try:
            from ..implicit_geometry import decompose_cylindrical_formula_to_vertical_bounds_2d

            cells = tuple(decompose_cylindrical_formula_to_vertical_bounds_2d(expr, variables))
            if cells:
                return expr, "vertical_bounds_2d", True
        except _RECOVERABLE_ERRORS:
            pass
    return expr, "cad", None


def _point_formula(point: Mapping[sp.Symbol, sp.Expr], variables: Sequence[sp.Symbol]) -> sp.Expr:
    """Return an exact conjunction describing one finite solution point."""

    atoms = [sp.Eq(var, sp.sympify(point[var])) for var in variables]
    return sp.And(*atoms) if atoms else sp.true


def _finite_points_formula(
    points: Sequence[Mapping[sp.Symbol, sp.Expr]], variables: Sequence[sp.Symbol]
) -> sp.Expr:
    """Return an exact finite-set formula from point assignments."""

    if not points:
        return sp.false
    formulas = [_point_formula(point, variables) for point in points]
    return sp.Or(*formulas) if len(formulas) > 1 else formulas[0]


def _try_rur_formula(
    formula: sp.Expr,
    variables: Sequence[sp.Symbol],
    *,
    max_solutions: int | None = None,
):
    """Try the exact RUR finite-system backend for a Boolean formula.

    Returns ``None`` when the formula is outside the supported finite equality
    fragment. A returned object with ``status == 'unknown'`` means RUR saw at
    least one unsupported branch, so callers must not use it as an UNSAT proof.
    """

    if not variables:
        return None
    try:
        result = solve_formula_with_rur(
            formula, tuple(variables), real=True, max_solutions=max_solutions
        )
    except _RECOVERABLE_ERRORS:
        return None
    if result is None:
        return None
    if str(result.status).lower() == "unknown":
        return None
    if result.partial and not result.assignments:
        return None
    return result


@dataclass
class _PlannedSolveOutcome:
    reduced: sp.Expr
    satisfiable: bool
    method: str
    diagnostics: dict[str, object] = field(default_factory=dict)
    solved: object | None = None
    rur_result: object | None = None


def _planner_step(
    method: str, accepted: bool, reason: str, **metadata: object
) -> dict[str, object]:
    return {
        "method": method,
        "accepted": bool(accepted),
        "reason": reason,
        "metadata": metadata,
    }


def _run_auto_solve_plan(
    expr: sp.Expr,
    variables: tuple[sp.Symbol, ...],
    *,
    domain: str,
    strategy: str | None,
    sample_count: int,
    depth: int = 0,
) -> _PlannedSolveOutcome:
    """Execute the structural exact solver plan before general CAD/QE.

    Every accepted transformation is equivalence preserving.  The planner may
    prove satisfiability without CAD only when the selected backend itself gives
    an exact certificate (interval reduction, Fourier--Motzkin, RUR, explicit
    cylindrical cells, or recursively solved independent blocks).
    """

    if depth > 12:
        parsed = ParsedPrenexFormula(variables, (), parse_formula(expr), expr)
        solved = reduce_formula(parsed, domain=domain, return_result=True, strategy=strategy)
        reduced = _safe_simplify_expr(solved.result)
        sat = is_satisfiable(reduced, variables, domain=domain, strategy=strategy)
        return _PlannedSolveOutcome(
            reduced,
            bool(sat),
            getattr(solved, "method", "cad"),
            {"planner_steps": (_planner_step("cad_qe", True, "planner recursion guard"),)},
            solved=solved,
        )

    steps: list[dict[str, object]] = []
    original_vars = tuple(variables)
    profile = profile_semialgebraic_system(expr, original_vars)
    diagnostics: dict[str, object] = {
        "system_profile": {
            "atom_count": profile.atom_count,
            "equality_count": profile.equality_count,
            "inequality_count": profile.inequality_count,
            "polynomial": profile.polynomial,
            "linear": profile.linear,
            "max_total_degree": profile.max_total_degree,
            "conjunctive": profile.conjunctive,
            "variable_blocks": profile.variable_blocks,
            "suggested_order": profile.suggested_order,
            "potentially_zero_dimensional": profile.potentially_zero_dimensional,
        }
    }

    work_expr = expr
    work_vars = original_vars
    presolved = None
    if profile.polynomial and profile.conjunctive and original_vars:
        try:
            presolved = affine_presolve(expr, original_vars)
        except _RECOVERABLE_ERRORS:
            presolved = None
        if presolved is not None and presolved.changed:
            work_expr = presolved.formula
            work_vars = presolved.variables
            steps.append(
                _planner_step(
                    "affine_presolve",
                    True,
                    "globally safe affine equality substitution",
                    substitutions=presolved.substitutions,
                    remaining_variables=work_vars,
                )
            )
        else:
            steps.append(_planner_step("affine_presolve", False, "no safe affine elimination"))

    def reconstruct(formula: sp.Expr) -> sp.Expr:
        if presolved is None or not presolved.substitutions:
            return formula
        return reconstruct_affine_solution_formula(presolved, formula)

    if is_false_expr(work_expr):
        steps.append(_planner_step("trivial", True, "presolve proved inconsistency"))
        diagnostics["planner_steps"] = tuple(steps)
        return _PlannedSolveOutcome(sp.false, False, "affine_presolve", diagnostics)
    if is_true_expr(work_expr) or not work_vars:
        if not is_true_expr(work_expr):
            simplified = _safe_simplify_expr(work_expr)
            if is_false_expr(simplified):
                steps.append(_planner_step("trivial", True, "zero-variable residual is false"))
                diagnostics["planner_steps"] = tuple(steps)
                return _PlannedSolveOutcome(sp.false, False, "affine_presolve", diagnostics)
            if not is_true_expr(simplified):
                # This should only occur if undeclared parameters escaped the
                # solve-variable list.  Leave the general reducer to handle it.
                work_expr = simplified
            else:
                work_expr = sp.true
        if is_true_expr(work_expr):
            result_formula = reconstruct(sp.true)
            steps.append(_planner_step("trivial", True, "all remaining constraints eliminated"))
            diagnostics["planner_steps"] = tuple(steps)
            return _PlannedSolveOutcome(
                result_formula,
                True,
                "affine_presolve" if presolved and presolved.changed else "trivial",
                diagnostics,
            )

    # Small explicit disjunctions are solved branch-by-branch.  This mirrors
    # Solve's DNF dispatch while avoiding an exponential normalization step: we
    # only split an Or that is already present and cap the branch count.
    if isinstance(work_expr, sp.Or) and len(work_expr.args) <= 8:
        branch_outcomes: list[_PlannedSolveOutcome] = []
        for branch in work_expr.args:
            branch_outcomes.append(
                _run_auto_solve_plan(
                    branch,
                    work_vars,
                    domain=domain,
                    strategy=strategy,
                    sample_count=sample_count,
                    depth=depth + 1,
                )
            )
        sat_branches = [branch for branch in branch_outcomes if branch.satisfiable]
        reduced_union = (
            sp.Or(*(branch.reduced for branch in sat_branches)) if sat_branches else sp.false
        )
        reduced_union = reconstruct(reduced_union) if sat_branches else sp.false
        methods = tuple(branch.method for branch in branch_outcomes)
        steps.append(
            _planner_step(
                "boolean_branch_decomposition",
                True,
                "small explicit disjunction solved branch-by-branch",
                branch_count=len(branch_outcomes),
                branch_methods=methods,
            )
        )
        diagnostics["planner_steps"] = tuple(steps)
        diagnostics["branch_plans"] = tuple(branch.diagnostics for branch in branch_outcomes)
        # Preserve the familiar finite-system method when every surviving
        # branch was solved by RUR; otherwise expose the decomposition.
        method = (
            "rational_univariate"
            if sat_branches
            and all(branch.method == "rational_univariate" for branch in sat_branches)
            else "boolean_branch_decomposition[" + ",".join(methods) + "]"
        )
        return _PlannedSolveOutcome(reduced_union, bool(sat_branches), method, diagnostics)

    # Solve independent variable-incidence blocks separately.  This is the
    # semialgebraic analogue of Solve's block-diagonal system decomposition.
    components = decompose_conjunctive_formula(work_expr, work_vars)
    if components:
        child_formulas: list[sp.Expr] = []
        child_methods: list[str] = []
        child_diags: list[dict[str, object]] = []
        all_sat = True
        for component_formula, component_vars in components:
            child = _run_auto_solve_plan(
                component_formula,
                tuple(component_vars),
                domain=domain,
                strategy=strategy,
                sample_count=sample_count,
                depth=depth + 1,
            )
            child_formulas.append(child.reduced)
            child_methods.append(child.method)
            child_diags.append(child.diagnostics)
            if not child.satisfiable:
                all_sat = False
                break
        if all_sat:
            combined = sp.And(*child_formulas) if child_formulas else sp.true
            combined = reconstruct(combined)
        else:
            combined = sp.false
        method = "incidence_decomposition[" + ",".join(child_methods) + "]"
        steps.append(
            _planner_step(
                "incidence_decomposition",
                True,
                "independent variable blocks solved separately",
                block_count=len(components),
                child_methods=tuple(child_methods),
            )
        )
        diagnostics["planner_steps"] = tuple(steps)
        diagnostics["child_plans"] = tuple(child_diags)
        return _PlannedSolveOutcome(combined, all_sat, method, diagnostics)
    steps.append(_planner_step("incidence_decomposition", False, "system is not separable"))

    # Preserve the finite-system solver contract: when the original system can
    # be zero dimensional, run RUR on the affine-presolved system before the
    # generic univariate-set reducer.  Affine presolve can make this RUR much
    # smaller while the result remains an exact finite-point certificate.
    if profile.potentially_zero_dimensional:
        rur_result = _try_rur_formula(
            work_expr,
            work_vars,
            max_solutions=sample_count if sample_count else None,
        )
        if rur_result is not None and not rur_result.partial:
            assignments = [dict(point) for point in rur_result.assignments]
            if presolved is not None and presolved.substitutions:
                assignments = [
                    reconstruct_affine_assignment(presolved, point) for point in assignments
                ]
            reduced = (
                _finite_points_formula(assignments, original_vars) if assignments else sp.false
            )
            steps.append(
                _planner_step(
                    "rational_univariate",
                    True,
                    "finite system solved by RUR after affine presolve",
                    solution_count=len(assignments),
                )
            )
            diagnostics["planner_steps"] = tuple(steps)
            return _PlannedSolveOutcome(
                reduced,
                bool(assignments),
                "rational_univariate",
                diagnostics,
                rur_result=rur_result,
            )
        steps.append(
            _planner_step(
                "rational_univariate",
                False,
                "finite-system RUR attempt did not certify the branch",
            )
        )

    # Univariate semialgebraic reduction is substantially cheaper than either
    # RUR or a general CAD and already yields an exact set description.
    if len(work_vars) == 1:
        reduced, fast_method, fast_sat = _fast_solution_formula(work_expr, work_vars)
        if fast_sat is not None:
            steps.append(_planner_step(fast_method, True, "exact univariate reduction"))
            diagnostics["planner_steps"] = tuple(steps)
            return _PlannedSolveOutcome(
                reconstruct(reduced if fast_sat else sp.false),
                bool(fast_sat),
                fast_method,
                diagnostics,
            )
        steps.append(
            _planner_step("univariate_reduction", False, "univariate fast path unsupported")
        )

    work_profile = profile_semialgebraic_system(work_expr, work_vars)

    # After affine equality substitution, a linear conjunction usually consists
    # solely of inequalities.  Fourier--Motzkin can then decide feasibility
    # without constructing a CAD; the unreduced linear formula itself is an
    # exact representation of the solution set.
    if work_profile.linear and work_profile.conjunctive:
        linear = exact_linear_feasibility(work_expr, work_vars)
        if linear is not None:
            linear_sat, certificate = linear
            steps.append(
                _planner_step(
                    "linear_fourier_motzkin",
                    True,
                    "exact linear feasibility after affine presolve",
                    eliminated_condition=certificate,
                )
            )
            diagnostics["planner_steps"] = tuple(steps)
            return _PlannedSolveOutcome(
                reconstruct(work_expr) if linear_sat else sp.false,
                bool(linear_sat),
                "linear_fourier_motzkin",
                diagnostics,
            )
        steps.append(
            _planner_step("linear_fourier_motzkin", False, "linear elimination unsupported")
        )

    # Explicit triangular/cylindrical bounds are already a solved set form and
    # prove nonemptiness through their extracted cells.
    if len(work_vars) > 1:
        try:
            from ..cad_algorithms.cells import extract_explicit_cylindrical_solution

            explicit_cyl = extract_explicit_cylindrical_solution(work_expr, work_vars)
        except _RECOVERABLE_ERRORS:
            explicit_cyl = None
        if explicit_cyl is not None:
            steps.append(
                _planner_step(
                    "explicit_cylindrical_bounds",
                    True,
                    "recognized exact triangular/cylindrical bound structure",
                )
            )
            diagnostics["planner_steps"] = tuple(steps)
            return _PlannedSolveOutcome(
                reconstruct(work_expr), True, "explicit_cylindrical_bounds", diagnostics
            )
        steps.append(_planner_step("explicit_cylindrical_bounds", False, "no explicit bound form"))

    # Use the cheap equation-count test only as a gate for exact ideal
    # analysis.  Once a Groebner basis is justified, decide finite-vs-positive
    # dimensionality from the leading monomial ideal rather than guessing from
    # the number of equations.
    exact_ideal = None
    if work_profile.potentially_zero_dimensional and work_profile.conjunctive:
        equality_residuals = tuple(
            sp.expand(atom.lhs - atom.rhs)
            for atom in conjuncts(work_expr)
            if isinstance(atom, sp.Equality) and sp.expand(atom.lhs - atom.rhs) != 0
        )
        try:
            exact_ideal = EqualityIdealContext(equality_residuals, work_vars)
        except _RECOVERABLE_ERRORS:
            exact_ideal = None

    if exact_ideal is not None and exact_ideal.inconsistent:
        steps.append(
            _planner_step(
                "equality_ideal_analysis",
                True,
                "equality ideal is the unit ideal",
                dimension=-1,
                quotient_dimension=0,
            )
        )
        diagnostics["planner_steps"] = tuple(steps)
        return _PlannedSolveOutcome(sp.false, False, "equality_ideal_analysis", diagnostics)

    should_try_rur = work_profile.potentially_zero_dimensional
    if exact_ideal is not None:
        should_try_rur = exact_ideal.zero_dimensional
        steps.append(
            _planner_step(
                "equality_ideal_analysis",
                True,
                "computed exact dimension from the leading monomial ideal",
                dimension=exact_ideal.dimension,
                quotient_dimension=exact_ideal.quotient_dimension,
            )
        )

    if should_try_rur:
        rur_result = _try_rur_formula(
            work_expr,
            work_vars,
            max_solutions=sample_count if sample_count else None,
        )
        if rur_result is not None and not rur_result.partial:
            assignments = [dict(point) for point in rur_result.assignments]
            if presolved is not None and presolved.substitutions:
                assignments = [
                    reconstruct_affine_assignment(presolved, point) for point in assignments
                ]
            reduced = (
                _finite_points_formula(assignments, original_vars) if assignments else sp.false
            )
            steps.append(
                _planner_step(
                    "rational_univariate",
                    True,
                    "finite polynomial equality system solved by RUR",
                    solution_count=len(assignments),
                    quotient_dimension=(
                        exact_ideal.quotient_dimension if exact_ideal is not None else None
                    ),
                )
            )
            diagnostics["planner_steps"] = tuple(steps)
            return _PlannedSolveOutcome(
                reduced,
                bool(assignments),
                "rational_univariate",
                diagnostics,
                rur_result=rur_result,
            )
        steps.append(
            _planner_step("rational_univariate", False, "RUR did not certify a finite system")
        )
    else:
        reason = (
            f"exact equality-ideal dimension is {exact_ideal.dimension}"
            if exact_ideal is not None
            else "equation count proves a finite algebraic set is impossible or not applicable"
        )
        steps.append(_planner_step("rational_univariate", False, reason))

    # A bounded Groebner presolve simplifies equality varieties and every
    # inequality modulo the equality ideal.  For positive-dimensional systems
    # this often shrinks the later CAD without pretending the answer is a list
    # of isolated points.
    groebner = None
    if work_profile.polynomial and work_profile.conjunctive and work_profile.equality_count:
        groebner = groebner_reduce_conjunction(
            work_expr,
            work_vars,
            order=work_profile.suggested_order,
        )
    if groebner is not None:
        groebner_expr, groebner_meta = groebner
        steps.append(
            _planner_step(
                "groebner_presolve",
                True,
                "reduced polynomial system modulo equality ideal",
                **groebner_meta,
            )
        )
        if is_false_expr(groebner_expr):
            diagnostics["planner_steps"] = tuple(steps)
            return _PlannedSolveOutcome(sp.false, False, "groebner_presolve", diagnostics)
        work_expr = groebner_expr
        # If Groebner reduction exposed a one-dimensional/linear residual,
        # exploit it before the general reducer.
        post_profile = profile_semialgebraic_system(work_expr, work_vars)
        if len(work_vars) == 1:
            reduced, fast_method, fast_sat = _fast_solution_formula(work_expr, work_vars)
            if fast_sat is not None:
                steps.append(_planner_step(fast_method, True, "post-Groebner univariate reduction"))
                diagnostics["planner_steps"] = tuple(steps)
                return _PlannedSolveOutcome(
                    reconstruct(reduced if fast_sat else sp.false),
                    bool(fast_sat),
                    "groebner_presolve+" + fast_method,
                    diagnostics,
                )
        if post_profile.linear and post_profile.conjunctive:
            linear = exact_linear_feasibility(work_expr, work_vars)
            if linear is not None:
                linear_sat, certificate = linear
                steps.append(
                    _planner_step(
                        "linear_fourier_motzkin",
                        True,
                        "post-Groebner system is linear",
                        eliminated_condition=certificate,
                    )
                )
                diagnostics["planner_steps"] = tuple(steps)
                return _PlannedSolveOutcome(
                    reconstruct(work_expr) if linear_sat else sp.false,
                    bool(linear_sat),
                    "groebner_presolve+linear_fourier_motzkin",
                    diagnostics,
                )
    else:
        steps.append(_planner_step("groebner_presolve", False, "not useful/applicable"))

    parsed = ParsedPrenexFormula(work_vars, (), parse_formula(work_expr), work_expr)
    solved = reduce_formula(parsed, domain=domain, return_result=True, strategy=strategy)
    reduced = _safe_simplify_expr(solved.result)
    try:
        satisfiable = is_satisfiable(reduced, work_vars, domain=domain, strategy=strategy)
    except PolynomialError:
        satisfiable = reduced is not sp.false and reduced != sp.false
    final_reduced = reconstruct(reduced) if satisfiable else sp.false
    selected = getattr(solved, "method", "cad")
    steps.append(_planner_step("cad_qe", True, "general exact fallback", backend=selected))
    diagnostics["planner_steps"] = tuple(steps)
    return _PlannedSolveOutcome(
        final_reduced,
        bool(satisfiable),
        selected,
        diagnostics,
        solved=solved,
    )


def _make_quantified_sentence(
    formula: sp.Expr, variables: Sequence[sp.Symbol]
) -> tuple[tuple[str, sp.Symbol], ...]:
    return tuple(("exists", var) for var in variables)


def _truth_from_qe_result(result) -> bool:
    if result.is_sentence:
        return bool(result.truth_value)
    simplified = _safe_simplify_expr(result.formula)
    if is_true_expr(simplified):
        return True
    if is_false_expr(simplified):
        return False
    # A non-sentence result means parameters escaped the requested variable set.
    # Treat satisfiability existentially over remaining free symbols.
    remaining = tuple(sorted(simplified.free_symbols, key=symbol_identity_key))
    if not remaining:
        return bool(simplified)
    return bool(
        qe_by_complete_cad(
            remaining,
            _make_quantified_sentence(simplified, remaining),
            parse_formula(simplified),
            return_result=True,
        ).truth_value
    )


from ._predicates import equivalent, implies, is_satisfiable, is_tautology  # noqa: E402


def _parameter_solution_data(
    formula: sp.Expr,
    variables: tuple[sp.Symbol, ...],
    parameters: tuple[sp.Symbol, ...],
    *,
    domain: str,
) -> tuple[sp.Expr | None, object | None]:
    """Compute exact parameter solvability data when parameters are present."""

    if not parameters:
        return None, None
    pieces = conjuncts(formula)
    param_formula = pieces[0] if len(pieces) == 1 else formula
    conditions = _fast_parameter_conditions(param_formula, variables, parameters)
    if conditions is None:
        from ..parameters import solvability_conditions

        conditions = solvability_conditions(param_formula, variables, parameters, domain=domain)
    decomposition = None
    try:
        from ..parameter_stratification import parameterized_cylindrical_decomposition

        decomposition = parameterized_cylindrical_decomposition(
            formula,
            variables,
            parameters,
            domain=domain,
            specialize_fibers=True,
        )
        if conditions is None:
            conditions = decomposition.parameter_condition
    except (TypeError, ValueError, ArithmeticError, NotImplementedError, PolynomialError):
        decomposition = None
    return conditions, decomposition


def _trivial_solution(
    formula: sp.Expr,
    variables: tuple[sp.Symbol, ...],
    parameters: tuple[sp.Symbol, ...],
    *,
    satisfiable: bool,
    sample_count: int,
    parameter_conditions: sp.Expr | None,
    parameter_decomposition: object | None,
    strategy: str | None,
) -> SemialgebraicSolution:
    """Build the structured result for a constant true or false formula."""

    result_formula = sp.true if satisfiable else sp.false
    samples = (
        tuple({var: sp.Integer(0) for var in variables} for _ in range(1 if sample_count else 0))
        if satisfiable
        else ()
    )
    meta = _collect_solution_metadata(result_formula, variables)
    return SemialgebraicSolution(
        result_formula,
        variables,
        samples,
        satisfiable,
        "trivial",
        solution_capability_diagnostics(formula),
        parameters=parameters,
        simplified_constraints=(),
        parameter_conditions=(
            parameter_conditions
            if parameter_conditions is not None
            else (sp.true if satisfiable else sp.false)
        ),
        parameter_decomposition=parameter_decomposition if satisfiable else None,
        dimension=meta["dimension"],
        bounded=meta["bounded"],
        closed=meta["closed"],
        compact=meta["compact"],
        components=meta["components"],
        cells=meta["cells"],
        cylindrical_solution=meta.get("cylindrical_solution"),
        connectivity=meta.get("connectivity"),
    )


from ._solve_api import solve_semialgebraic  # noqa: E402


def canonicalize_one_dimensional_formula(
    formula: FormulaLike | Iterable[FormulaLike], variable: sp.Symbol | str
) -> sp.Expr:
    """Return a canonical interval-union formula for supported 1D systems."""

    expr = normalize_formula(formula)
    var = _as_real_symbol(variable)
    components = _one_dim_components(expr, var)
    if components is None:
        reduced = reduce_conjunctive_inequalities(expr, var)
        components = _one_dim_components(reduced, var) if reduced is not None else None
    if components is None:
        return _safe_simplify_expr(expr)
    return _components_formula(components)


__all__ = [
    "IntervalComponent",
    "SemialgebraicSolution",
    "EquivalenceResult",
    "ImplicationResult",
    "TautologyResult",
    "SatisfiabilityResult",
    "equivalent",
    "implies",
    "is_satisfiable",
    "is_tautology",
    "canonicalize_one_dimensional_formula",
    "solve_semialgebraic",
]
