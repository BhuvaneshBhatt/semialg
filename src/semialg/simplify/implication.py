from __future__ import annotations

from collections.abc import Iterable
from functools import lru_cache

import sympy as sp
from sympy.core.relational import Relational
from sympy.logic.boolalg import And as SymAnd

from ..formulas.boolean import bounded_dnf_branches, make_and


def _dnf_branches(expr: sp.Expr) -> tuple[tuple[sp.Expr, ...], ...]:
    """Return a DNF branch representation for a SymPy Boolean expression."""

    simplified = sp.simplify_logic(expr, form="dnf")
    expansion = bounded_dnf_branches(simplified, max_branches=4096)
    if not expansion.complete:
        return (_branch_atoms(simplified),)
    branches = []
    for branch in expansion.branches:
        expr_branch = [item for item in branch if item is not sp.true and item is not True]
        if any(item is sp.false or item is False for item in branch):
            continue
        branches.append(_branch_atoms(make_and(*expr_branch)))
    return tuple(branches)


def _branch_atoms(expr: sp.Expr) -> tuple[sp.Expr, ...]:
    if expr is sp.true:
        return tuple()
    if isinstance(expr, SymAnd):
        atoms = expr.args
    else:
        atoms = (expr,)
    return tuple(sorted(dict.fromkeys(atoms), key=sp.sstr))


def _branch_expr(branch: Iterable[sp.Expr]) -> sp.Expr:
    items = tuple(branch)
    if not items:
        return sp.true
    return sp.And(*items)


def _safe_not(expr: sp.Expr) -> sp.Expr:
    return sp.Not(expr)


def _complexity_ok(expr: sp.Expr, *, max_atoms: int, max_vars: int, max_degree: int) -> bool:
    atoms = tuple(expr.atoms(Relational))
    if len(atoms) > max_atoms:
        return False
    vars_ = sorted(expr.free_symbols, key=lambda s: s.name)
    if len(vars_) > max_vars:
        return False
    for atom in atoms:
        try:
            poly = sp.Poly(sp.expand(atom.lhs - atom.rhs), *vars_)
        except Exception:
            return False
        if poly.total_degree() > max_degree:
            return False
    return True


@lru_cache(maxsize=512)
def _is_unsatisfiable_cached(expr_text: str, symbol_names: tuple[str, ...]) -> bool | None:
    symbols = {name: sp.Symbol(name, real=True) for name in symbol_names}
    try:
        expr = sp.sympify(expr_text, locals=symbols)
    except Exception:
        return None
    try:
        from ..cad.decomposition import decomp_collins_complete
        from ..formula import parse_formula
        from ..qe.complete import evaluate_formula_on_cell

        formula = parse_formula(expr)
        vars_ = tuple(symbols[name] for name in symbol_names)
        polys = []
        from ..formula import formula_polynomials

        polys = tuple(formula_polynomials(formula)) or (sp.Integer(1),)
        cad = decomp_collins_complete(polys, vars_)
        level = len(vars_)
        return not any(
            evaluate_formula_on_cell(formula, cell, vars_)
            for cell in cad.cells_by_level.get(level, ())
        )
    except Exception:
        return None


def is_unsatisfiable_by_cad(
    expr: sp.Expr, *, max_atoms: int = 8, max_vars: int = 3, max_degree: int = 4
) -> bool | None:
    """Conservatively decide small real-polynomial unsatisfiability by CAD.

    ``None`` means the implication minimizer should leave the formula alone.
    The size guard prevents the pretty-printer from accidentally launching an
    expensive CAD on every large output branch.
    """

    if expr is sp.false:
        return True
    if expr is sp.true:
        return False
    if not _complexity_ok(expr, max_atoms=max_atoms, max_vars=max_vars, max_degree=max_degree):
        return None
    names = tuple(sorted(sym.name for sym in expr.free_symbols))
    return _is_unsatisfiable_cached(sp.sstr(expr), names)


def implies_by_cad(antecedent: sp.Expr, consequent: sp.Expr, **kwargs) -> bool | None:
    """Return whether ``antecedent => consequent`` for guarded small formulas."""

    return is_unsatisfiable_by_cad(sp.And(antecedent, _safe_not(consequent)), **kwargs)


def minimize_conj_by_impl(branch: tuple[sp.Expr, ...], **kwargs) -> tuple[sp.Expr, ...] | None:
    """Remove atoms implied by the other atoms in a conjunction."""

    atoms = list(branch)
    changed = False
    idx = 0
    while idx < len(atoms):
        atom = atoms[idx]
        rest = atoms[:idx] + atoms[idx + 1 :]
        if not rest:
            idx += 1
            continue
        result = implies_by_cad(_branch_expr(rest), atom, **kwargs)
        if result is True:
            atoms.pop(idx)
            changed = True
            continue
        idx += 1
    return tuple(atoms) if changed else None


def minimize_disj_by_impl(expr: sp.Expr, **kwargs) -> sp.Expr:
    """Remove CAD-provably redundant DNF branches and redundant branch atoms.

    This is a genuine implication-based minimizer for small semialgebraic
    formulas. It is intentionally guarded and conservative: when CAD cannot
    cheaply prove an implication, the original branch is kept.
    """

    branches = [tuple(branch) for branch in _dnf_branches(expr)]
    if not branches:
        return sp.false
    # First minimize atoms inside each branch.
    for pos, branch in enumerate(tuple(branches)):
        minimized = minimize_conj_by_impl(branch, **kwargs)
        if minimized is not None:
            branches[pos] = minimized
    # Remove duplicate branches after atom minimization.
    unique: list[tuple[sp.Expr, ...]] = []
    for branch in branches:
        if branch not in unique:
            unique.append(branch)
    branches = unique
    keep = [True] * len(branches)
    for i, left in enumerate(branches):
        if not keep[i]:
            continue
        left_expr = _branch_expr(left)
        for j, right in enumerate(branches):
            if i == j or not keep[j]:
                continue
            right_expr = _branch_expr(right)
            result = implies_by_cad(left_expr, right_expr, **kwargs)
            if result is True:
                keep[i] = False
                break
    kept = [_branch_expr(branch) for branch, flag in zip(branches, keep, strict=True) if flag]
    if not kept:
        return sp.false
    return sp.simplify_logic(sp.Or(*kept), form="dnf")


__all__ = [
    "implies_by_cad",
    "is_unsatisfiable_by_cad",
    "minimize_conj_by_impl",
    "minimize_disj_by_impl",
]
