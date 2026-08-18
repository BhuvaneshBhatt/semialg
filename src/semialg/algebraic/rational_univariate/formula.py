from __future__ import annotations

from collections.abc import Mapping, Sequence

import sympy as sp

from ...formulas.boolean import bounded_dnf_branches, is_false_expr, is_true_expr
from .representation import RationalUnivariateError, RationalUnivariateFormulaResult
from .signs import evaluate_boolean_formula_at_point, solve_rur_semialgebraic_system


def _formula_disjuncts(
    formula: sp.Expr | bool,
    *,
    max_branches: int = 64,
) -> tuple[tuple[sp.Expr | bool, ...], ...]:
    """Return bounded DNF-style conjunction branches."""

    expansion = bounded_dnf_branches(formula, max_branches=max_branches)
    if not expansion.complete:
        raise RationalUnivariateError(
            "; ".join(expansion.notes) or "bounded Boolean branch expansion exceeded the RUR limit"
        )
    return expansion.branches


def _branch_equalities_and_constraints(
    branch: Sequence[sp.Expr | bool],
    variables: Sequence[sp.Symbol],
) -> tuple[tuple[sp.Expr, ...], sp.Expr | bool] | None:
    equalities: list[sp.Expr] = []
    constraints: list[sp.Expr | bool] = []
    variable_set = set(variables)
    for atom in branch:
        if is_true_expr(atom):
            continue
        if is_false_expr(atom):
            return tuple(), sp.false
        if isinstance(atom, sp.Equality):
            residual = sp.expand(atom.lhs - atom.rhs)
            if residual == 0:
                continue
            if residual.is_number:
                return tuple(), sp.false
            if residual.free_symbols <= variable_set:
                try:
                    sp.Poly(residual, *variables, domain=sp.QQ)
                except (sp.PolynomialError, ValueError, TypeError):
                    constraints.append(atom)
                else:
                    equalities.append(residual)
                    continue
        constraints.append(atom)
    constraint_formula: sp.Expr | bool
    if not constraints:
        constraint_formula = sp.true
    elif len(constraints) == 1:
        constraint_formula = constraints[0]
    else:
        constraint_formula = sp.And(*constraints, evaluate=False)
    return tuple(equalities), constraint_formula


def solve_formula_with_rur(
    formula: sp.Expr | bool,
    variables: Sequence[sp.Symbol],
    *,
    real: bool = True,
    max_solutions: int | None = None,
) -> RationalUnivariateFormulaResult | None:
    """Try to solve a finite Boolean formula by RUR branch enumeration.

    Each disjunctive branch must contain enough rational polynomial equalities
    in ``variables`` to define a zero-dimensional candidate set. Remaining
    relations are evaluated exactly at the algebraic candidate points. The
    function returns ``None`` when no branch is in the supported RUR fragment;
    an empty result with status ``"unsat"`` means at least one branch was
    supported and all supported branches were unsatisfiable.
    """

    variable_tuple = tuple(variables)
    if not variable_tuple:
        truth = evaluate_boolean_formula_at_point(formula, {})
        return RationalUnivariateFormulaResult(
            variables=variable_tuple,
            assignments=(({},) if truth else ()),
            status="satisfied" if truth else "unsat",
            solved_branches=1,
        )

    assignments: list[Mapping[sp.Symbol, sp.Expr]] = []
    solved_branches = 0
    skipped_branches = 0
    notes: list[str] = []
    seen: set[tuple[str, ...]] = set()
    try:
        branches = _formula_disjuncts(formula)
    except RationalUnivariateError as exc:
        return RationalUnivariateFormulaResult(
            variables=variable_tuple,
            assignments=tuple(),
            status="unknown",
            solved_branches=0,
            skipped_branches=1,
            notes=(str(exc),),
        )

    for branch in branches:
        parsed = _branch_equalities_and_constraints(branch, variable_tuple)
        if parsed is None:
            skipped_branches += 1
            continue
        equalities, constraints = parsed
        if is_false_expr(constraints):
            solved_branches += 1
            continue
        if len(equalities) < len(variable_tuple):
            skipped_branches += 1
            notes.append("skipped branch without enough rational equalities for RUR")
            continue
        try:
            branch_solutions = solve_rur_semialgebraic_system(
                equalities,
                variable_tuple,
                constraints,
                real=real,
                as_assignments=True,
            )
        except RationalUnivariateError as exc:
            skipped_branches += 1
            notes.append(f"skipped branch outside RUR fragment: {exc}")
            continue
        solved_branches += 1
        for assignment in branch_solutions:  # type: ignore[assignment]
            key = tuple(sp.sstr(sp.simplify(assignment[var])) for var in variable_tuple)
            if key not in seen:
                seen.add(key)
                assignments.append(dict(assignment))
                if max_solutions is not None and len(assignments) >= max_solutions:
                    return RationalUnivariateFormulaResult(
                        variables=variable_tuple,
                        assignments=tuple(assignments),
                        status="satisfied",
                        solved_branches=solved_branches,
                        skipped_branches=skipped_branches,
                        notes=tuple(notes),
                    )
    if solved_branches == 0:
        return None
    if skipped_branches and not assignments:
        return RationalUnivariateFormulaResult(
            variables=variable_tuple,
            assignments=tuple(),
            status="unknown",
            solved_branches=solved_branches,
            skipped_branches=skipped_branches,
            notes=tuple(notes) + ("RUR coverage was partial; skipped branches may be satisfiable",),
        )
    return RationalUnivariateFormulaResult(
        variables=variable_tuple,
        assignments=tuple(assignments),
        status="partial" if skipped_branches else ("satisfied" if assignments else "unsat"),
        solved_branches=solved_branches,
        skipped_branches=skipped_branches,
        notes=tuple(notes),
    )
