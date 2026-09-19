from __future__ import annotations

from collections.abc import Iterable
from dataclasses import dataclass

import sympy as sp
from sympy.core.relational import (
    Equality,
    GreaterThan,
    LessThan,
    Relational,
    StrictGreaterThan,
    StrictLessThan,
    Unequality,
)
from sympy.logic.boolalg import And as SymAnd
from sympy.logic.boolalg import BooleanFalse, BooleanTrue
from sympy.logic.boolalg import Not as SymNot
from sympy.logic.boolalg import Or as SymOr

RELATION_TYPES = (
    sp.Equality,
    sp.Unequality,
    sp.StrictLessThan,
    sp.LessThan,
    sp.StrictGreaterThan,
    sp.GreaterThan,
)


@dataclass(frozen=True)
class BooleanBranchExpansion:
    branches: tuple[tuple[sp.Expr | bool, ...], ...]
    complete: bool = True
    truncated: bool = False
    notes: tuple[str, ...] = ()


def is_true_expr(expr: object) -> bool:
    return expr is True or expr is sp.true or expr == sp.true or expr is sp.S.true


def is_false_expr(expr: object) -> bool:
    return expr is False or expr is sp.false or expr == sp.false or expr is sp.S.false


def make_and(*args: sp.Expr) -> sp.Expr:
    filtered = [arg for arg in args if not is_true_expr(arg)]
    if any(is_false_expr(arg) for arg in args):
        return sp.false
    if not filtered:
        return sp.true
    return sp.And(*filtered, evaluate=False)


def make_or(*args: sp.Expr) -> sp.Expr:
    filtered = [arg for arg in args if not is_false_expr(arg)]
    if any(is_true_expr(arg) for arg in args):
        return sp.true
    if not filtered:
        return sp.false
    return sp.Or(*filtered, evaluate=False)


def relation_residual(relation: sp.Expr) -> sp.Expr:
    return sp.expand(relation.lhs - relation.rhs)  # type: ignore[attr-defined]


def canonical_relation(atom: Relational) -> tuple[sp.Expr, str]:
    if isinstance(atom, Equality):
        return sp.expand(atom.lhs - atom.rhs), "="
    if isinstance(atom, Unequality):
        return sp.expand(atom.lhs - atom.rhs), "!="
    if isinstance(atom, StrictLessThan):
        return sp.expand(atom.lhs - atom.rhs), "<"
    if isinstance(atom, LessThan):
        return sp.expand(atom.lhs - atom.rhs), "<="
    if isinstance(atom, StrictGreaterThan):
        return sp.expand(atom.rhs - atom.lhs), "<"
    if isinstance(atom, GreaterThan):
        return sp.expand(atom.rhs - atom.lhs), "<="
    raise ValueError(f"unsupported relational atom: {atom!r}")


def relation_from_residual(expr: sp.Expr, operator: str) -> sp.Expr:
    expr = sp.expand(expr)
    if operator == "=":
        return sp.Eq(expr, 0)
    if operator == "!=":
        return sp.Ne(expr, 0)
    if operator == "<":
        return expr < 0
    if operator == "<=":
        return expr <= 0
    raise ValueError(f"unsupported relation operator: {operator!r}")


def negate_relation(atom: Relational) -> sp.Expr:
    expr, operator = canonical_relation(atom)
    if operator == "=":
        return relation_from_residual(expr, "!=")
    if operator == "!=":
        return relation_from_residual(expr, "=")
    if operator == "<":
        return relation_from_residual(-expr, "<=")
    if operator == "<=":
        return relation_from_residual(-expr, "<")
    raise ValueError(f"unsupported relation operator: {operator!r}")


def to_negation_normal_form(formula: sp.Expr, negate: bool = False) -> sp.Expr:
    if is_true_expr(formula) or isinstance(formula, BooleanTrue):
        return sp.false if negate else sp.true
    if is_false_expr(formula) or isinstance(formula, BooleanFalse):
        return sp.true if negate else sp.false
    if isinstance(formula, Relational):
        return negate_relation(formula) if negate else formula
    if isinstance(formula, SymNot):
        return to_negation_normal_form(formula.args[0], not negate)
    if isinstance(formula, SymAnd):
        mapped = tuple(to_negation_normal_form(arg, negate) for arg in formula.args)
        return make_or(*mapped) if negate else make_and(*mapped)
    if isinstance(formula, SymOr):
        mapped = tuple(to_negation_normal_form(arg, negate) for arg in formula.args)
        return make_and(*mapped) if negate else make_or(*mapped)
    raise ValueError(f"unsupported Boolean formula node: {formula!r}")


def iter_relational_atoms(formula: sp.Expr) -> Iterable[Relational]:
    if is_true_expr(formula) or is_false_expr(formula):
        return
    if isinstance(formula, Relational):
        yield formula
        return
    if isinstance(formula, (SymAnd, SymOr)):
        for arg in formula.args:
            yield from iter_relational_atoms(arg)
        return
    raise ValueError(f"unsupported Boolean formula node: {formula!r}")


def bounded_dnf_branches(
    formula: sp.Expr | bool, *, max_branches: int = 64
) -> BooleanBranchExpansion:
    notes: list[str] = []

    def expand(node: sp.Expr | bool) -> list[tuple[sp.Expr | bool, ...]]:
        if is_true_expr(node) or is_false_expr(node) or isinstance(node, bool):
            return [(node,)]
        if isinstance(node, sp.Or):
            branches: list[tuple[sp.Expr | bool, ...]] = []
            for arg in node.args:
                branches.extend(expand(arg))
                if len(branches) > max_branches:
                    notes.append("bounded Boolean branch expansion exceeded the limit")
                    raise OverflowError(notes[-1])
            return branches
        if isinstance(node, sp.And):
            branches: list[tuple[sp.Expr | bool, ...]] = [tuple()]
            for arg in node.args:
                arg_branches = expand(arg)
                expanded: list[tuple[sp.Expr | bool, ...]] = []
                for prefix in branches:
                    for suffix in arg_branches:
                        expanded.append(prefix + suffix)
                        if len(expanded) > max_branches:
                            notes.append("bounded Boolean branch expansion exceeded the limit")
                            raise OverflowError(notes[-1])
                branches = expanded
            return branches
        return [(node,)]

    try:
        return BooleanBranchExpansion(
            tuple(expand(formula)), complete=True, truncated=False, notes=tuple(notes)
        )
    except OverflowError:
        return BooleanBranchExpansion(tuple(), complete=False, truncated=True, notes=tuple(notes))
