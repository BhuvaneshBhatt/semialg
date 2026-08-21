from __future__ import annotations

from collections.abc import Sequence

import sympy as sp
from sympy import S
from sympy.logic.boolalg import And as SymAnd
from sympy.solvers.inequalities import reduce_inequalities

from ..algebraic.rational_univariate import RationalUnivariateError, solve_formula_with_rur
from ..cad.constants import (
    PROJECTION_COLLINS,
    PROJECTION_LAZARD,
    PROJECTION_MCCALLUM,
    PROJECTION_TTICAD,
)
from ..cad.reduced import decomp_form_reduced_safe
from ..context import with_computation_context
from ..formula import ParsedPrenexFormula, parse_formula, parse_quant_form_text, to_sympy
from ..model import ProjectionConfig
from ..planner.select import select_strat_for_form
from ..qe.complete import qe_by_complete_cad, qe_from_cad
from ..qe.virtual_substitution import try_quadratic_virtual_substitution_qe
from ..simplify.boolean import simplify_boolean
from ..status import SolverStatus
from ..tticad.safe import decompose_tticad_safe
from .domains import SolveDomain, apply_assumptions, normalize_assumptions, normalize_domain
from .integer.diophantine import solve_int_methods
from .integer.linear_divisibility import detect_lin_reduction
from .preprocess import semialgebraicize
from .result import SolveResult


def parse_and_maybe_prep(
    text: str, *, symbols=None, variable_order=None, use_preprocess: bool = True
):
    parsed = parse_quant_form_text(text, symbols=symbols, variable_order=variable_order)
    preprocess_changed = False
    preprocess_meta = None
    if use_preprocess:
        prep = semialgebraicize(parsed.matrix)
        preprocess_changed = prep.changed
        preprocess_meta = prep
        if prep.changed:
            parsed = ParsedPrenexFormula(
                tuple(parsed.vars) + tuple(prep.aux_vars),
                tuple(parsed.quantifiers) + tuple(("exists", aux) for aux in prep.aux_vars),
                prep.formula,
                prep.sympy_expr,
            )
    return parsed, preprocess_changed, preprocess_meta


def _finite_assignments_to_formula(assignments, free_variables: Sequence[sp.Symbol]) -> sp.Expr:
    free_variables = tuple(free_variables)
    if not assignments:
        return sp.false
    if not free_variables:
        return sp.true
    pieces = []
    seen: set[tuple[str, ...]] = set()
    for assignment in assignments:
        key = tuple(sp.sstr(sp.simplify(assignment[var])) for var in free_variables)
        if key in seen:
            continue
        seen.add(key)
        pieces.append(
            sp.And(
                *(sp.Eq(var, sp.simplify(assignment[var])) for var in free_variables),
                evaluate=False,
            )
        )
    if not pieces:
        return sp.false
    return pieces[0] if len(pieces) == 1 else sp.Or(*pieces, evaluate=False)


def _try_rational_univariate_reduction(parsed: ParsedPrenexFormula):
    """Use RUR when the full formula has a finite algebraic solution set.

    For existential formulas this projects finite points to the free variables;
    for quantifier-free formulas it returns an equivalent finite disjunction of
    point equalities. Returning ``None`` means the formula was outside the RUR
    fragment and should continue to VS/CAD.
    """

    if parsed.quantifiers and any(str(q).lower() != "exists" for q, _ in parsed.quantifiers):
        return None
    quantified = {sym for _, sym in parsed.quantifiers}
    matrix_symbols = tuple(
        sorted(getattr(parsed.matrix_expr, "free_symbols", set()), key=lambda sym: sym.name)
    )
    ordered_symbols = tuple(dict.fromkeys(tuple(parsed.vars) + matrix_symbols))
    all_symbols = tuple(
        sym for sym in ordered_symbols if sym in matrix_symbols or sym in quantified
    )
    if not all_symbols:
        return None
    try:
        rur = solve_formula_with_rur(parsed.matrix_expr, all_symbols, real=True)
    except RationalUnivariateError:
        return None
    if rur is None or not rur.complete or rur.status == SolverStatus.UNKNOWN:
        return None
    free_variables = tuple(sym for sym in all_symbols if sym not in quantified)
    formula = _finite_assignments_to_formula(rur.assignments, free_variables)
    try:
        formula = sp.simplify_logic(formula, form="dnf")
    except (TypeError, ValueError, sp.SympifyError):
        pass
    return formula, rur


def _planner_ready_parsed(parsed: ParsedPrenexFormula) -> ParsedPrenexFormula:
    """Rebuild Formula nodes before planner analysis.

    Some test runners reload ``semialg.formula`` while keeping imported parser
    objects alive. Re-parsing from the SymPy expression avoids class-identity
    mismatches and gives the planner a normalized formula tree.
    """

    matrix_expr = to_sympy(parsed.matrix)
    return ParsedPrenexFormula(
        tuple(parsed.vars),
        tuple(parsed.quantifiers),
        parse_formula(matrix_expr),
        matrix_expr,
    )


def _safe_strategy_selection(parsed: ParsedPrenexFormula, fallback_vars: Sequence[sp.Symbol]):
    """Select a CAD backend, falling back to conservative Collins on planner issues."""

    from ..planner.select import StrategySelection

    try:
        planner_parsed = _planner_ready_parsed(parsed)
        return select_strat_for_form(planner_parsed.matrix, parsed=planner_parsed)
    except Exception as exc:
        return StrategySelection(
            backend=PROJECTION_COLLINS,
            variable_order=tuple(fallback_vars),
            partial=False,
            projection=ProjectionConfig(operator=PROJECTION_COLLINS),
            notes=(f"Planner selection failed; using conservative Collins fallback: {exc}",),
        )


def _reduce_reals(parsed: ParsedPrenexFormula, config=None, *, strategy: str | None = None):
    """Dispatch real quantifier elimination through the configured exact strategy stack."""
    selection = None
    vars_ = parsed.vars
    strategy_name = (strategy or "collins").lower()
    if strategy_name in {"auto", "planner", "rur", "rational-univariate", "rational_univariate"}:
        rur_reduction = _try_rational_univariate_reduction(parsed)
        if rur_reduction is not None:
            formula, rur_result = rur_reduction
            return formula, rur_result, selection
        if strategy_name in {"rur", "rational-univariate", "rational_univariate"}:
            raise NotImplementedError(
                "formula is outside the rational-univariate finite-solver fragment"
            )
    use_vs_prepass = strategy_name in {
        "auto",
        "planner",
        "virtual-substitution",
        "virtual_substitution",
        "vs",
    }
    if use_vs_prepass and parsed.quantifiers:
        vs_result = try_quadratic_virtual_substitution_qe(
            vars_,
            parsed.quantifiers,
            to_sympy(parsed.matrix),
            full=True,
        )
        if vs_result is not None and not vs_result.remaining_quantifiers:
            return vs_result.formula, vs_result, selection
        if vs_result is not None and vs_result.remaining_quantifiers:
            reduced_matrix = parse_formula(vs_result.formula)
            remaining_symbols = tuple(
                sym for sym in getattr(vs_result.formula, "free_symbols", set())
            )
            remaining_quantified = tuple(sym for _, sym in vs_result.remaining_quantifiers)
            free_after_vs = tuple(
                sym for sym in vars_ if sym in remaining_symbols and sym not in remaining_quantified
            )
            extra_free = tuple(
                sorted(
                    (set(remaining_symbols) - set(free_after_vs) - set(remaining_quantified)),
                    key=lambda s: s.name,
                )
            )
            parsed = ParsedPrenexFormula(
                free_after_vs + extra_free + remaining_quantified,
                vs_result.remaining_quantifiers,
                reduced_matrix,
                vs_result.formula,
            )
            vars_ = parsed.vars
    if strategy_name in {"auto", "planner"}:
        selection = _safe_strategy_selection(parsed, vars_)
        vars_ = selection.variable_order
    elif strategy_name in {"mccallum", "lazard", "tticad"}:
        base = _safe_strategy_selection(parsed, vars_)
        vars_ = base.variable_order
        from ..planner.select import StrategySelection

        backend = {
            "mccallum": PROJECTION_MCCALLUM,
            "lazard": PROJECTION_LAZARD,
            "tticad": PROJECTION_TTICAD,
        }[strategy_name]
        selection = StrategySelection(
            backend=backend,
            variable_order=vars_,
            partial=base.partial,
            projection=base.projection,
            notes=(*base.notes, f"Explicit strategy requested: {strategy_name}."),
        )
    if selection is not None and selection.backend != PROJECTION_COLLINS:
        if selection.backend == PROJECTION_TTICAD:
            safe = decompose_tticad_safe(parsed.matrix, vars_)
        else:
            backend_name = "lazard" if selection.backend == PROJECTION_LAZARD else "mccallum"
            safe = decomp_form_reduced_safe(parsed.matrix, vars_, backend=backend_name)
        result = qe_from_cad(
            safe.cad,
            vars_,
            parsed.quantifiers,
            parsed.matrix,
            backend=getattr(safe, "effective_backend", safe.cad.backend),
        )
        from ..planner.select import StrategySelection

        selection = StrategySelection(
            backend=getattr(safe, "effective_backend", selection.backend),
            variable_order=vars_,
            partial=selection.partial,
            projection=selection.projection,
            notes=(
                *selection.notes,
                f"Reduced backend used_fallback={getattr(safe, 'used_fallback', None)}.",
            ),
        )
    else:
        result = qe_by_complete_cad(vars_, parsed.quantifiers, parsed.matrix)
    if result.is_sentence:
        formula = sp.true if result.truth_value else sp.false
    else:
        formula = result.formula
    if "vs_result" in locals() and vs_result is not None and vs_result.remaining_quantifiers:
        if selection is not None:
            from ..planner.select import StrategySelection

            selection = StrategySelection(
                backend=selection.backend,
                variable_order=vars_,
                partial=selection.partial,
                projection=selection.projection,
                notes=(*selection.notes, *vs_result.notes),
            )
    return formula, result, selection


def _complex_formula_support(expr: sp.Expr) -> bool:
    return not bool(
        expr.atoms(sp.StrictLessThan, sp.StrictGreaterThan, sp.LessThan, sp.GreaterThan)
    )


def _finite_set_to_formula(vars_, solset):
    if solset is S.EmptySet:
        return sp.false
    if isinstance(solset, sp.FiniteSet):
        pieces = []
        for item in solset:
            if len(vars_) == 1 and not isinstance(item, tuple):
                item = (item,)
            elif not isinstance(item, tuple):
                continue
            pieces.append(sp.And(*[sp.Eq(v, val) for v, val in zip(vars_, item, strict=True)]))
        if not pieces:
            return sp.false
        return sp.Or(*pieces) if len(pieces) > 1 else pieces[0]
    return None


def _reduce_complexes(parsed: ParsedPrenexFormula):
    expr = to_sympy(parsed.matrix)
    vars_ = tuple(parsed.vars)
    if parsed.quantifiers:
        # Limited experimental support: existential quantifier over a quantifier-free equality system.
        if all(q == "exists" for q, _ in parsed.quantifiers) and _complex_formula_support(expr):
            if len(vars_) == 1:
                sol = sp.solveset(expr, vars_[0], domain=S.Complexes)
                return sp.Ne(sol, S.EmptySet), {"set": sol}
        raise NotImplementedError(
            "Complex quantified reduction is implemented only for very small existential cases"
        )
    if not _complex_formula_support(expr):
        raise NotImplementedError("Complex-domain inequalities are not supported")
    if len(vars_) == 1:
        sol = sp.solveset(expr, vars_[0], domain=S.Complexes)
        formula = _finite_set_to_formula(vars_, sol)
        if formula is None:
            formula = sp.Contains(vars_[0], sol)
        return formula, {"set": sol}
    try:
        sol = sp.nonlinsolve(
            [a.lhs - a.rhs if isinstance(a, sp.Equality) else a for a in sp.And.make_args(expr)],
            vars_,
        )
    except (sp.PolynomialError, ValueError, TypeError):
        sol = None
    if sol is not None:
        formula = _finite_set_to_formula(vars_, sol)
        if formula is not None:
            return formula, {"set": sol}
    return sp.simplify(expr), {"note": "complex symbolic passthrough"}


def _integer_set_to_formula(var, solset):
    if solset is S.EmptySet:
        return sp.false
    if isinstance(solset, sp.FiniteSet):
        vals = sorted(solset, key=sp.default_sort_key)
        pieces = [sp.Eq(var, v) for v in vals]
        return sp.Or(*pieces) if len(pieces) > 1 else pieces[0]
    return sp.And(sp.Contains(var, S.Integers), sp.Contains(var, solset))


def reduce_int_univar(expr: sp.Expr, var: sp.Symbol):
    args = sp.And.make_args(expr) if isinstance(expr, SymAnd) else (expr,)
    sets = []
    for arg in args:
        if isinstance(arg, sp.Equality):
            sets.append(sp.solveset(arg, var, domain=S.Integers))
        elif isinstance(
            arg, (sp.StrictLessThan, sp.StrictGreaterThan, sp.LessThan, sp.GreaterThan)
        ):
            r = reduce_inequalities([arg], var)
            sets.append(sp.solveset(r, var, domain=S.Integers))
        else:
            sets.append(sp.solveset(arg, var, domain=S.Integers))
    sol = sets[0]
    for s in sets[1:]:
        sol = sol.intersect(s)
    return _integer_set_to_formula(var, sol), {"set": sol}


def _reduce_integers(parsed: ParsedPrenexFormula):
    expr = to_sympy(parsed.matrix)
    vars_ = tuple(parsed.vars)
    if parsed.quantifiers:
        if len(vars_) == 1:
            var = vars_[0]
            formula, meta = reduce_int_univar(expr, var)
            if all(q == "exists" for q, _ in parsed.quantifiers):
                solset = meta["set"]
                return sp.true if solset is not S.EmptySet else sp.false, meta
            raise NotImplementedError(
                "Only existential univariate integer quantification is supported"
            )
        raise NotImplementedError(
            "Integer-domain quantified reduction is not implemented beyond simple univariate existential cases"
        )

    specialized = solve_int_methods(expr, vars_)
    if specialized is not None:
        return specialized.formula, {"specialized_result": specialized}

    reduction = detect_lin_reduction(expr, vars_)
    if reduction is not None:
        reduced_vars = tuple(v for v in vars_ if v != reduction.solved_variable)
        meta = {
            "linear_divisibility_reduction": {
                "solved_variable": sp.sstr(reduction.solved_variable),
                "replacement": sp.sstr(reduction.replacement),
                "divisibility_condition": sp.sstr(reduction.divisibility_condition),
            }
        }
        if len(reduced_vars) == 0:
            return sp.simplify(
                sp.And(
                    sp.Eq(reduction.solved_variable, reduction.replacement),
                    reduction.reduced_formula,
                )
            ), meta
        if len(reduced_vars) == 1:
            subformula, submeta = reduce_int_univar(reduction.reduced_formula, reduced_vars[0])
            rebuilt = sp.And(sp.Eq(reduction.solved_variable, reduction.replacement), subformula)
            meta["subproblem"] = submeta
            return sp.simplify(rebuilt), meta
        rebuilt = sp.And(
            sp.Eq(reduction.solved_variable, reduction.replacement), reduction.reduced_formula
        )
        return sp.simplify(rebuilt), meta

    if len(vars_) == 1:
        return reduce_int_univar(expr, vars_[0])

    args = sp.And.make_args(expr) if isinstance(expr, SymAnd) else (expr,)
    equalities = [a for a in args if isinstance(a, sp.Equality)]
    if len(equalities) == 1 and len(args) == 1:
        try:
            solset = sp.diophantine(equalities[0].lhs - equalities[0].rhs)
            formula = (
                _finite_set_to_formula(vars_, sp.FiniteSet(*list(solset)[:10]))
                if solset
                else sp.false
            )
            return formula if formula is not None else expr, {"set": solset}
        except (ValueError, TypeError, sp.SympifyError):
            pass
    return sp.And(sp.Contains(sp.Tuple(*vars_), S.Integers ** len(vars_)), sp.simplify(expr)), {
        "note": "integer symbolic restriction"
    }


def _drop_trivial_real_bounds(formula: sp.Expr) -> sp.Expr:
    if not isinstance(formula, sp.And):
        return formula
    kept = []
    for arg in formula.args:
        if isinstance(arg, sp.StrictLessThan) and arg.rhs is sp.oo:
            continue
        if isinstance(arg, sp.StrictGreaterThan) and arg.rhs is -sp.oo:
            continue
        kept.append(arg)
    if not kept:
        return sp.true
    return sp.And(*kept)


def _normalize_reduced_real_formula(
    formula: sp.Expr, free_variables: Sequence[sp.Symbol]
) -> sp.Expr:
    """Lightweight cleanup for formulas produced by specialized QE backends."""

    if getattr(formula, "is_Boolean", False):
        simplified = simplify_boolean(formula)
    else:
        try:
            simplified = sp.simplify(formula)
        except Exception:
            simplified = formula
    try:
        symbols = tuple(
            sym for sym in free_variables if sym in getattr(simplified, "free_symbols", set())
        )
        if len(symbols) == 1:
            reduced = reduce_inequalities(
                list(sp.And(simplified).args) if isinstance(simplified, sp.And) else [simplified],
                symbols[0],
            )
            return simplify_boolean(_drop_trivial_real_bounds(reduced))
    except Exception:
        pass
    dropped = _drop_trivial_real_bounds(simplified)
    return (
        simplify_boolean(dropped) if getattr(dropped, "is_Boolean", False) else sp.simplify(dropped)
    )


@with_computation_context
def reduce_formula(
    parsed: ParsedPrenexFormula,
    config=None,
    *,
    domain: str | SolveDomain | None = None,
    assumptions=None,
    return_result: bool = False,
    strategy: str | None = None,
):
    dom = normalize_domain(domain)
    if assumptions:
        matrix = parse_formula(apply_assumptions(to_sympy(parsed.matrix), assumptions))
        parsed = ParsedPrenexFormula(
            parsed.vars,
            parsed.quantifiers,
            matrix,
            apply_assumptions(parsed.sympy_expr, assumptions),
        )
    if dom is SolveDomain.REALS:
        reduced, qe_result, selection = _reduce_reals(parsed, config=config, strategy=strategy)
        metadata = {"qe_result": qe_result}
        if selection is not None:
            metadata["strategy_selection"] = selection
        elif (strategy or "collins").lower() in {"auto", "planner"}:
            metadata["strategy_selection"] = _safe_strategy_selection(parsed, parsed.vars)
        method = (
            "quadratic_virtual_substitution"
            if getattr(qe_result, "backend", "") == "quadratic-virtual-substitution-qe"
            else (
                "rational_univariate"
                if getattr(qe_result, "backend", "") == "rational-univariate-formula-solver"
                else "cad"
            )
        )
        if method != "rational_univariate":
            reduced = _normalize_reduced_real_formula(
                reduced,
                tuple(sym for sym in parsed.vars if sym not in {v for _, v in parsed.quantifiers}),
            )
        solve_result = SolveResult(method=method, domain=dom, result=reduced, metadata=metadata)
    elif dom is SolveDomain.COMPLEXES:
        reduced, meta = _reduce_complexes(parsed)
        solve_result = SolveResult(
            method="symbolic_complex", domain=dom, result=reduced, metadata=meta
        )
    elif dom is SolveDomain.INTEGERS:
        reduced, meta = _reduce_integers(parsed)
        solve_result = SolveResult(
            method="symbolic_integer", domain=dom, result=reduced, metadata=meta
        )
    elif dom is SolveDomain.RATIONALS:
        # Rational QE is not a real-closed-field problem. Use the real result as
        # a conservative symbolic approximation and mark the domain in metadata.
        reduced, qe_result, selection = _reduce_reals(parsed, config=config, strategy=strategy)
        solve_result = SolveResult(
            method="rational_via_real_relaxation",
            domain=dom,
            result=reduced,
            metadata={"qe_result": qe_result, "domain_note": "computed over the real relaxation"},
        )
    else:
        solve_result = SolveResult(
            method="unsupported_domain",
            domain=dom,
            result=sp.false,
            metadata={"reason": f"unsupported domain {dom.value}"},
        )
    if assumptions:
        solve_result.metadata["assumptions"] = tuple(
            map(sp.sstr, normalize_assumptions(assumptions))
        )
    return solve_result if return_result else solve_result.result


def reduce_text(
    text: str,
    *,
    symbols=None,
    variable_order=None,
    config=None,
    domain: str | SolveDomain | None = None,
    assumptions=None,
    use_preprocess: bool = True,
    return_result: bool = False,
    strategy: str | None = None,
):
    if variable_order is not None:
        variable_order = [
            sp.Symbol(v, real=True) if isinstance(v, str) else v for v in variable_order
        ]
    parsed, preprocess_changed, preprocess_meta = parse_and_maybe_prep(
        text, symbols=symbols, variable_order=variable_order, use_preprocess=use_preprocess
    )
    solved = reduce_formula(
        parsed,
        config=config,
        domain=domain,
        assumptions=assumptions,
        return_result=True,
        strategy=strategy,
    )
    if hasattr(solved.result, "xreplace"):
        prefer_plain_symbols = (strategy or "").lower() == "mccallum"
        canonical_by_name = {
            symbol.name: (
                sp.Symbol(symbol.name)
                if prefer_plain_symbols
                else (symbol if symbol.is_real is True else sp.Symbol(symbol.name, real=True))
            )
            for symbol in parsed.vars
        }
        replacements = {
            symbol: canonical_by_name[symbol.name]
            for symbol in getattr(solved.result, "free_symbols", set())
            if symbol.name in canonical_by_name and symbol != canonical_by_name[symbol.name]
        }
        if replacements:
            solved.result = solved.result.xreplace(replacements)
    solved.normalized_text = text
    solved.preprocess_changed = preprocess_changed
    if preprocess_meta is not None:
        solved.metadata["preprocess"] = {
            "changed": preprocess_meta.changed,
            "notes": preprocess_meta.notes,
            "assumptions": tuple(map(sp.sstr, preprocess_meta.assumptions)),
        }
    return solved if return_result else solved.result


__all__ = ["reduce_formula", "reduce_text"]
