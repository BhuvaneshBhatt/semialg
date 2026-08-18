from __future__ import annotations

from dataclasses import dataclass

from ..formula import And, Formula, Or
from .formula_normalize import normalize_formula


@dataclass(frozen=True)
class FormulaBranches:
    normalized: Formula
    top_level_branches: tuple[Formula, ...]
    branch_count: int


def split_top_level_branches(formula: Formula) -> FormulaBranches:
    normalized = normalize_formula(formula)
    if isinstance(normalized, Or):
        return FormulaBranches(
            normalized=normalized,
            top_level_branches=normalized.args,
            branch_count=len(normalized.args),
        )
    return FormulaBranches(normalized=normalized, top_level_branches=(normalized,), branch_count=1)


def conjunctive_branches(formula: Formula) -> tuple[Formula, ...]:
    branches = split_top_level_branches(formula).top_level_branches
    out = []
    for branch in branches:
        if isinstance(branch, And):
            out.extend(branch.args)
        else:
            out.append(branch)
    return tuple(out)
