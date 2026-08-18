from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass

import sympy as sp

from .witness_generation import sample_free_assignments


@dataclass(frozen=True)
class PartialAssignExtension:
    base_assignment: dict[sp.Symbol, object]
    extensions: tuple[dict[sp.Symbol, object], ...]
    constraint: sp.Expr


def complete_partial_assign(
    partial_assignment: Mapping[sp.Symbol, object],
    variables: Sequence[sp.Symbol],
    *,
    domain_rules: Mapping[sp.Symbol, Sequence[object]] | None = None,
    sample_count: int = 5,
    modulus: int | None = None,
    seed: int | None = None,
) -> list[dict[sp.Symbol, object]]:
    partial_assignment = dict(partial_assignment)
    missing = [v for v in variables if v not in partial_assignment]
    if not missing:
        return [partial_assignment]
    fills = sample_free_assignments(
        missing, domain_rules=domain_rules, sample_count=sample_count, modulus=modulus, seed=seed
    )
    return [{**partial_assignment, **fill} for fill in fills]


def extend_partial_cons(
    partial_assignment: Mapping[sp.Symbol, object],
    constraint: sp.Expr,
    variables: Sequence[sp.Symbol],
    *,
    domain_rules: Mapping[sp.Symbol, Sequence[object]] | None = None,
    sample_count: int = 5,
    modulus: int | None = None,
    seed: int | None = None,
) -> PartialAssignExtension:
    partial_assignment = dict(partial_assignment)
    missing = [v for v in variables if v not in partial_assignment]
    valid: list[dict[sp.Symbol, object]] = []
    if len(missing) == 1 and modulus is None:
        try:
            solved = sp.solve(constraint.subs(partial_assignment), missing[0], dict=True)
        except Exception:
            solved = []
        for sol in solved:
            valid.append({**partial_assignment, missing[0]: sp.simplify(sol[missing[0]])})
    if not valid:
        candidates = complete_partial_assign(
            partial_assignment,
            variables,
            domain_rules=domain_rules,
            sample_count=sample_count,
            modulus=modulus,
            seed=seed,
        )
        for assignment in candidates:
            try:
                reduced = sp.simplify(constraint.subs(assignment))
            except Exception:
                continue
            if reduced is sp.true or reduced is True:
                valid.append(assignment)
    return PartialAssignExtension(dict(partial_assignment), tuple(valid), constraint)


__all__ = ["PartialAssignExtension", "complete_partial_assign", "extend_partial_cons"]
