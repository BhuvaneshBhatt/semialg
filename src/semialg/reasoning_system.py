from __future__ import annotations

from collections.abc import Iterable, Sequence

import sympy as sp
from sympy.logic.boolalg import Boolean

from ._errors import EXACT_OPERATION_ERRORS as _RECOVERABLE_ERRORS
from .decision import equivalent, implies, is_satisfiable
from .inequality_reduction import reduce_conjunctive_inequalities
from .normalization import conjuncts, normalize_formula, normalize_variables
from .reasoning_results import SimplifiedSystem
from .simplify.boolean import simplify_boolean
from .symbolic_simplify import simplify_boole

FormulaLike = sp.Expr | Boolean | bool


def _ordered_unique_exprs(exprs: Iterable[sp.Expr]) -> tuple[sp.Expr, ...]:
    out: list[sp.Expr] = []
    seen: set[sp.Expr] = set()
    for expr in exprs:
        simplified = sp.simplify(expr) if not getattr(expr, "is_Relational", False) else expr
        key = simplified
        if key not in seen:
            out.append(simplified)
            seen.add(key)
    return tuple(out)


def _try_reduce_inequalities_1d(expr: sp.Expr, variables: Sequence[sp.Symbol]) -> sp.Expr | None:
    if len(variables) != 1:
        return None
    var = variables[0]
    return reduce_conjunctive_inequalities(expr, var)


def _simple_equality_substitutions(
    atoms: Sequence[sp.Expr], variables: Sequence[sp.Symbol]
) -> dict[sp.Symbol, sp.Expr]:
    substitutions: dict[sp.Symbol, sp.Expr] = {}
    variable_set = set(variables)
    for atom in atoms:
        if not isinstance(atom, sp.Equality):
            continue
        pairs = ((atom.lhs, atom.rhs), (atom.rhs, atom.lhs))
        for lhs, rhs in pairs:
            if (
                isinstance(lhs, sp.Symbol)
                and lhs in variable_set
                and lhs not in getattr(rhs, "free_symbols", set())
            ):
                if lhs not in substitutions:
                    substitutions[lhs] = sp.simplify(rhs)
                break
    # Keep only acyclic substitutions in user-variable order.
    clean: dict[sp.Symbol, sp.Expr] = {}
    for sym in variables:
        if sym not in substitutions:
            continue
        rhs = substitutions[sym].xreplace(clean)
        if sym not in getattr(rhs, "free_symbols", set()):
            clean[sym] = rhs
    return clean


def _apply_substitutions_to_constraints(
    atoms: Sequence[sp.Expr], substitutions: dict[sp.Symbol, sp.Expr], *, keep_definitions: bool
) -> tuple[sp.Expr, ...]:
    out: list[sp.Expr] = []
    for atom in atoms:
        if isinstance(atom, sp.Equality):
            if atom.lhs in substitutions and sp.simplify(atom.rhs - substitutions[atom.lhs]) == 0:
                if keep_definitions:
                    out.append(atom)
                continue
            if atom.rhs in substitutions and sp.simplify(atom.lhs - substitutions[atom.rhs]) == 0:
                if keep_definitions:
                    out.append(atom)
                continue
        try:
            new_atom = atom.subs(substitutions)
            if new_atom in (True, sp.true):
                continue
            if new_atom in (False, sp.false):
                out.append(sp.false)
            else:
                out.append(new_atom)
        except _RECOVERABLE_ERRORS:
            out.append(atom)
    return tuple(out)


def simplify_system(
    constraints: FormulaLike | Iterable[FormulaLike],
    variables: Sequence[sp.Symbol | str] | None = None,
    *,
    assumptions: FormulaLike | Iterable[FormulaLike] = True,
    return_result: bool = False,
    strategy: str | None = None,
    eliminate_equalities: bool = False,
    output: str = "formula",
) -> sp.Expr | tuple[sp.Expr, ...] | SimplifiedSystem:
    """Simplify a real semialgebraic system with CAD/QE-backed checks.

    The routine focuses on dependable semantic simplifications:
    contradiction detection, duplicate removal, redundancy removal by implication,
    and one-dimensional inequality normalization via SymPy's real inequality
    reducer. It intentionally returns a mathematically equivalent formula rather
    than promising a globally minimal or prettiest CAD cell union.
    """

    if output not in {"formula", "constraints", "result"}:
        raise ValueError("output must be 'formula', 'constraints', or 'result'")
    expr = normalize_formula(constraints)
    asm = normalize_formula(assumptions)
    vars_ = normalize_variables(variables, sp.And(asm, expr))
    initial_atoms = _ordered_unique_exprs(conjuncts(expr))
    substitutions = _simple_equality_substitutions(initial_atoms, vars_)
    if substitutions:
        substituted_atoms = _apply_substitutions_to_constraints(
            initial_atoms, substitutions, keep_definitions=not eliminate_equalities
        )
        expr = sp.And(*substituted_atoms) if substituted_atoms else sp.true
    combined = sp.And(asm, expr) if asm is not sp.true and asm != sp.true else expr

    explicit_cyl = None
    try:
        from .cad_algorithms.cells import extract_explicit_cylindrical_solution

        explicit_cyl = extract_explicit_cylindrical_solution(combined, vars_)
    except _RECOVERABLE_ERRORS:
        explicit_cyl = None

    if combined is sp.false or combined == sp.false:
        result = SimplifiedSystem(
            sp.false,
            (),
            vars_,
            True,
            substitutions=substitutions,
            diagnostics={"reason": "inconsistent"},
        )
        return (
            result
            if (return_result or output == "result")
            else (() if output == "constraints" else sp.false)
        )
    if explicit_cyl is None and not is_satisfiable(combined, vars_, strategy=strategy):
        result = SimplifiedSystem(
            sp.false,
            (),
            vars_,
            True,
            substitutions=substitutions,
            diagnostics={"reason": "inconsistent"},
        )
        return (
            result
            if (return_result or output == "result")
            else (() if output == "constraints" else sp.false)
        )
    if explicit_cyl is not None and asm in (sp.true, True):
        formula = explicit_cyl.as_formula(closed=False)
        atoms = _ordered_unique_exprs(conjuncts(formula))
        result = SimplifiedSystem(
            formula=formula,
            constraints=atoms,
            variables=vars_,
            inconsistent=False,
            removed_redundant=(),
            substitutions=substitutions,
            method="explicit_cylindrical_cell_formula",
            diagnostics={
                "used_cylindrical_cell_formula": True,
                "cell_count": len(explicit_cyl.cells),
            },
        )
        return (
            result
            if (return_result or output == "result")
            else (atoms if output == "constraints" else formula)
        )

    normalized_1d = _try_reduce_inequalities_1d(expr, vars_)
    if normalized_1d is not None and equivalent(expr, normalized_1d, vars_, strategy=strategy):
        expr = normalized_1d

    try:
        expr = simplify_boole(expr, vars_, assumptions=asm, strategy=strategy)
    except _RECOVERABLE_ERRORS:
        try:
            expr = simplify_boolean(expr)
        except _RECOVERABLE_ERRORS:
            pass

    atoms = _ordered_unique_exprs(conjuncts(expr))
    if not atoms:
        result = SimplifiedSystem(
            sp.true, (), vars_, False, substitutions=substitutions, method="trivial"
        )
        if return_result or output == "result":
            return result
        return () if output == "constraints" else sp.true

    kept: list[sp.Expr] = list(atoms)
    removed: list[sp.Expr] = []
    changed = True
    while changed:
        changed = False
        for atom in list(kept):
            others = [other for other in kept if other is not atom]
            if not others:
                continue
            premise = sp.And(*(others + ([] if asm is sp.true or asm == sp.true else [asm])))
            try:
                if implies(premise, atom, vars_, strategy=strategy):
                    kept.remove(atom)
                    removed.append(atom)
                    changed = True
                    break
            except _RECOVERABLE_ERRORS:
                continue

    formula = sp.And(*kept) if kept else sp.true
    try:
        formula = simplify_boole(formula, vars_, assumptions=asm, strategy=strategy)
    except _RECOVERABLE_ERRORS:
        try:
            formula = simplify_boolean(formula)
        except _RECOVERABLE_ERRORS:
            pass
    cylindrical_formula_used = False
    if len(vars_) > 1:
        try:
            from .cad_algorithms.cells import extract_cylindrical_solution

            cyl = extract_cylindrical_solution(formula, vars_, selected_only=True)
            # Use a finite cylindrical formula only when it is small enough to
            # improve readability rather than explode the output. This lets the
            # simplifier exploit arbitrary-dimensional nested CAD cells while
            # staying conservative for large decompositions.
            if cyl.cells and len(cyl.cells) <= 12:
                cyl_formula = cyl.as_formula(closed=False)
                if explicit_cyl is not None or equivalent(
                    formula, cyl_formula, vars_, strategy=strategy
                ):
                    formula = cyl_formula
                    kept = list(conjuncts(formula))
                    cylindrical_formula_used = True
        except _RECOVERABLE_ERRORS:
            pass

    result = SimplifiedSystem(
        formula=formula,
        constraints=tuple(kept),
        variables=vars_,
        inconsistent=False,
        removed_redundant=tuple(removed),
        substitutions=substitutions,
        diagnostics={
            "input_constraint_count": len(atoms),
            "kept_constraint_count": len(kept),
            "used_cylindrical_cell_formula": cylindrical_formula_used,
        },
    )
    if return_result or output == "result":
        return result
    if output == "constraints":
        return tuple(kept)
    return formula


__all__ = ["SimplifiedSystem", "simplify_system"]
