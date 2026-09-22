from __future__ import annotations

from collections.abc import Iterable, Sequence
from dataclasses import dataclass, field

import sympy as sp
from sympy.logic.boolalg import Boolean
from sympy.polys.polyerrors import PolynomialError

from .._errors import EXACT_OPERATION_ERRORS as _RECOVERABLE_ERRORS
from ..algebraic.equality_ideal import EqualityIdealContext
from ..decision_diagnostics import solution_capability_diagnostics
from ..formula import ParsedPrenexFormula, parse_formula
from ..formulas.boolean import is_false_expr, is_true_expr
from ..incidence import decompose_conjunctive_formula
from ..inequality_reduction import reduce_conjunctive_inequalities
from ..normalization import conjuncts, normalize_formula
from ..qe import qe_by_complete_cad
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
from ._formula_support import (
    fast_parameter_conditions as _fast_parameter_conditions,
)
from ._formula_support import fast_solution_formula as _fast_solution_formula
from ._formula_support import safe_simplify_expr as _safe_simplify_expr
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


from ._rur import (  # noqa: E402
    finite_points_formula as _finite_points_formula,
)
from ._rur import (  # noqa: E402
    try_rur_formula as _try_rur_formula,
)


@dataclass
class _PlannedSolveOutcome:
    reduced: sp.Expr
    satisfiable: bool
    method: str
    diagnostics: dict[str, object] = field(default_factory=dict)
    solved: object | None = None
    rur_result: object | None = None


def _profile_diagnostics(profile) -> dict[str, object]:
    """Return the stable, serializable view of a solver-planning profile."""

    return {
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


def _completed_rur_outcome(
    rur_result,
    *,
    presolved,
    original_vars: tuple[sp.Symbol, ...],
    steps: list[dict[str, object]],
    diagnostics: dict[str, object],
    reason: str,
    metadata: dict[str, object] | None = None,
) -> _PlannedSolveOutcome:
    """Convert one complete finite RUR result into the planner's result model."""

    assignments = [dict(point) for point in rur_result.assignments]
    if presolved is not None and presolved.substitutions:
        assignments = [reconstruct_affine_assignment(presolved, point) for point in assignments]
    reduced = _finite_points_formula(assignments, original_vars) if assignments else sp.false
    step_metadata: dict[str, object] = {"solution_count": len(assignments)}
    if metadata:
        step_metadata.update(metadata)
    steps.append(_planner_step("rational_univariate", True, reason, **step_metadata))
    diagnostics["planner_steps"] = tuple(steps)
    return _PlannedSolveOutcome(
        reduced,
        bool(assignments),
        "rational_univariate",
        diagnostics,
        rur_result=rur_result,
    )


def _planner_step(
    method: str, accepted: bool, reason: str, **metadata: object
) -> dict[str, object]:
    return {
        "method": method,
        "accepted": bool(accepted),
        "reason": reason,
        "metadata": metadata,
    }


@dataclass
class _SolvePlanState:
    """Mutable per-call state shared by exact planner stages.

    This object is local to one solve.  It makes stage ownership
    explicit without adding caching or changing the order of mathematical work.
    """

    original_expr: sp.Expr
    original_vars: tuple[sp.Symbol, ...]
    work_expr: sp.Expr
    work_vars: tuple[sp.Symbol, ...]
    profile: object
    presolved: object | None
    steps: list[dict[str, object]]
    diagnostics: dict[str, object]

    def reconstruct(self, formula: sp.Expr) -> sp.Expr:
        if self.presolved is None or not self.presolved.substitutions:
            return formula
        return reconstruct_affine_solution_formula(self.presolved, formula)


def _prepare_solve_plan(expr: sp.Expr, variables: tuple[sp.Symbol, ...]) -> _SolvePlanState:
    """Profile and affine-presolve a planner request exactly once."""

    steps: list[dict[str, object]] = []
    profile = profile_semialgebraic_system(expr, variables)
    diagnostics: dict[str, object] = {"system_profile": _profile_diagnostics(profile)}
    work_expr = expr
    work_vars = variables
    presolved = None
    if profile.polynomial and profile.conjunctive and variables:
        try:
            presolved = affine_presolve(expr, variables)
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
    return _SolvePlanState(
        expr, variables, work_expr, work_vars, profile, presolved, steps, diagnostics
    )


def _finish_with_general_reducer(
    state: _SolvePlanState, *, domain: str, strategy: str | None
) -> _PlannedSolveOutcome:
    """Run the unchanged general CAD/QE fallback and normalize its outcome."""

    parsed = ParsedPrenexFormula(
        state.work_vars, (), parse_formula(state.work_expr), state.work_expr
    )
    solved = reduce_formula(parsed, domain=domain, return_result=True, strategy=strategy)
    reduced = _safe_simplify_expr(solved.result)
    try:
        satisfiable = is_satisfiable(reduced, state.work_vars, domain=domain, strategy=strategy)
    except PolynomialError:
        satisfiable = reduced is not sp.false and reduced != sp.false
    final_reduced = state.reconstruct(reduced) if satisfiable else sp.false
    selected = getattr(solved, "method", "cad")
    state.steps.append(_planner_step("cad_qe", True, "general exact fallback", backend=selected))
    state.diagnostics["planner_steps"] = tuple(state.steps)
    return _PlannedSolveOutcome(
        final_reduced, bool(satisfiable), selected, state.diagnostics, solved=solved
    )


def _try_structural_solve_decomposition(
    state: _SolvePlanState,
    *,
    domain: str,
    strategy: str | None,
    sample_count: int,
    depth: int,
) -> _PlannedSolveOutcome | None:
    """Solve explicit Boolean branches or independent incidence blocks exactly."""

    work_expr, work_vars = state.work_expr, state.work_vars
    steps, diagnostics, reconstruct = state.steps, state.diagnostics, state.reconstruct
    if isinstance(work_expr, sp.Or) and len(work_expr.args) <= 8:
        branch_outcomes = [
            _run_auto_solve_plan(
                branch,
                work_vars,
                domain=domain,
                strategy=strategy,
                sample_count=sample_count,
                depth=depth + 1,
            )
            for branch in work_expr.args
        ]
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
        method = (
            "rational_univariate"
            if sat_branches
            and all(branch.method == "rational_univariate" for branch in sat_branches)
            else "boolean_branch_decomposition[" + ",".join(methods) + "]"
        )
        return _PlannedSolveOutcome(reduced_union, bool(sat_branches), method, diagnostics)

    components = decompose_conjunctive_formula(work_expr, work_vars)
    if not components:
        steps.append(_planner_step("incidence_decomposition", False, "system is not separable"))
        return None
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
    combined = (
        reconstruct(sp.And(*child_formulas) if child_formulas else sp.true) if all_sat else sp.false
    )
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

    state = _prepare_solve_plan(expr, tuple(variables))
    steps = state.steps
    original_vars = state.original_vars
    profile = state.profile
    diagnostics = state.diagnostics
    work_expr = state.work_expr
    work_vars = state.work_vars
    presolved = state.presolved
    reconstruct = state.reconstruct

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

    decomposed = _try_structural_solve_decomposition(
        state, domain=domain, strategy=strategy, sample_count=sample_count, depth=depth
    )
    if decomposed is not None:
        return decomposed

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
            return _completed_rur_outcome(
                rur_result,
                presolved=presolved,
                original_vars=original_vars,
                steps=steps,
                diagnostics=diagnostics,
                reason="finite system solved by RUR after affine presolve",
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
            return _completed_rur_outcome(
                rur_result,
                presolved=presolved,
                original_vars=original_vars,
                steps=steps,
                diagnostics=diagnostics,
                reason="finite polynomial equality system solved by RUR",
                metadata={
                    "quotient_dimension": (
                        exact_ideal.quotient_dimension if exact_ideal is not None else None
                    )
                },
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

    state.work_expr = work_expr
    state.work_vars = work_vars
    return _finish_with_general_reducer(state, domain=domain, strategy=strategy)


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
        from ..decomposition import parametric_cad

        decomposition = parametric_cad(
            formula,
            variables,
            parameters=parameters,
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
