from __future__ import annotations

from dataclasses import dataclass

import sympy as sp

from ..formula import And, Atom, BoolConst, Formula, Not, Or


def _formula_truth_status(formula: Formula, subs: dict[sp.Symbol, sp.Expr]):
    if isinstance(formula, BoolConst):
        return formula.value
    if isinstance(formula, Atom):
        value = sp.expand(formula.expr).subs(subs)
        if value.free_symbols:
            return None
        value = sp.simplify(value)
        if formula.op == "=":
            return bool(value == 0)
        if formula.op == "!=":
            return bool(value != 0)
        numeric = sp.N(value, 50)
        if formula.op == "<":
            return bool(numeric < 0)
        if formula.op == "<=":
            return bool(numeric <= 0)
        if formula.op == ">":
            return bool(numeric > 0)
        if formula.op == ">=":
            return bool(numeric >= 0)
        raise ValueError(f"Unsupported operator: {formula.op}")
    if isinstance(formula, And):
        saw_unknown = False
        for arg in formula.args:
            value = _formula_truth_status(arg, subs)
            if value is False:
                return False
            if value is None:
                saw_unknown = True
        return None if saw_unknown else True
    if isinstance(formula, Or):
        saw_unknown = False
        for arg in formula.args:
            value = _formula_truth_status(arg, subs)
            if value is True:
                return True
            if value is None:
                saw_unknown = True
        return None if saw_unknown else False
    if isinstance(formula, Not):
        value = _formula_truth_status(formula.arg, subs)
        return None if value is None else (not value)
    raise TypeError(f"Unsupported formula node: {type(formula)}")


@dataclass(frozen=True)
class PruningDecision:
    status: bool | None
    should_prune: bool


def evaluate_pruning_status(formula: Formula, subs: dict[sp.Symbol, sp.Expr]) -> PruningDecision:
    status = _formula_truth_status(formula, subs)
    return PruningDecision(status=status, should_prune=status is not None)


__all__ = ["PruningDecision", "evaluate_pruning_status"]
