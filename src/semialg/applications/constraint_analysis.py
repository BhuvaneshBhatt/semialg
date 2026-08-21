"""Exact redundancy and feasibility diagnostics for semialgebraic constraints."""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field

import sympy as sp

from ..decision import implies, is_satisfiable
from ..normalization import normalize_formula, normalize_variables

FormulaLike = sp.Expr | sp.logic.boolalg.Boolean | bool


@dataclass(frozen=True)
class ConstraintRedundancyResult:
    """Exact redundancy classification for an explicit constraint list."""

    constraints: tuple[sp.Expr, ...]
    variables: tuple[sp.Symbol, ...]
    feasible: bool
    redundant_indices: tuple[int, ...]
    essential_indices: tuple[int, ...]
    witnesses: Mapping[int, Mapping[sp.Symbol, sp.Expr] | None] = field(default_factory=dict)
    method: str = "exact_implication_redundancy"
    certified: bool = True


@dataclass(frozen=True)
class FeasibleSetDiagnosticResult:
    """Exact feasibility diagnostics with a witness or irreducible conflict."""

    constraints: tuple[sp.Expr, ...]
    variables: tuple[sp.Symbol, ...]
    feasible: bool
    witness: Mapping[sp.Symbol, sp.Expr] | None
    redundant_indices: tuple[int, ...] = ()
    essential_indices: tuple[int, ...] = ()
    conflict_indices: tuple[int, ...] = ()
    method: str = "exact_semialgebraic_feasibility_diagnostics"
    certified: bool = True


def _normalize_constraints(
    constraints: Sequence[FormulaLike], variables: Sequence[sp.Symbol | str] | None
) -> tuple[tuple[sp.Expr, ...], tuple[sp.Symbol, ...]]:
    if isinstance(constraints, (str, bytes)):
        raise TypeError("constraints must be a sequence of formulas")
    pieces = tuple(normalize_formula(item) for item in constraints)
    if not pieces:
        raise ValueError("at least one constraint is required")
    vars_ = normalize_variables(variables, sp.And(*pieces))
    return pieces, vars_


def analyze_constraint_redundancy(
    constraints: Sequence[FormulaLike],
    variables: Sequence[sp.Symbol | str] | None = None,
) -> ConstraintRedundancyResult:
    """Determine which constraints are implied by all the others exactly."""

    pieces, vars_ = _normalize_constraints(constraints, variables)
    full = sp.And(*pieces)
    sat = is_satisfiable(full, vars_, return_result=True)
    redundant: list[int] = []
    essential: list[int] = []
    witnesses: dict[int, Mapping[sp.Symbol, sp.Expr] | None] = {}
    for index, constraint in enumerate(pieces):
        others = sp.And(*(piece for pos, piece in enumerate(pieces) if pos != index))
        result = implies(others, constraint, vars_, return_result=True)
        if bool(result):
            redundant.append(index)
        else:
            essential.append(index)
            witnesses[index] = getattr(result, "counterexample", None)
    return ConstraintRedundancyResult(
        constraints=pieces,
        variables=vars_,
        feasible=bool(sat),
        redundant_indices=tuple(redundant),
        essential_indices=tuple(essential),
        witnesses=witnesses,
    )


def _irreducible_conflict(
    constraints: tuple[sp.Expr, ...], variables: tuple[sp.Symbol, ...]
) -> tuple[int, ...]:
    active = list(range(len(constraints)))
    changed = True
    while changed:
        changed = False
        for index in tuple(active):
            candidate = [pos for pos in active if pos != index]
            formula = sp.And(*(constraints[pos] for pos in candidate))
            if not bool(is_satisfiable(formula, variables)):
                active = candidate
                changed = True
    return tuple(active)


def diagnose_feasible_set(
    constraints: Sequence[FormulaLike],
    variables: Sequence[sp.Symbol | str] | None = None,
    *,
    find_conflict: bool = True,
) -> FeasibleSetDiagnosticResult:
    """Diagnose a semialgebraic constraint set exactly.

    Feasible systems return an exact witness plus redundancy information.
    Infeasible systems optionally return an inclusion-minimal conflicting set
    of constraint indices.  The conflict is irreducible, not guaranteed to
    have minimum cardinality.
    """

    pieces, vars_ = _normalize_constraints(constraints, variables)
    sat = is_satisfiable(sp.And(*pieces), vars_, return_result=True)
    if bool(sat):
        redundancy = analyze_constraint_redundancy(pieces, vars_)
        return FeasibleSetDiagnosticResult(
            constraints=pieces,
            variables=vars_,
            feasible=True,
            witness=sat.witness,
            redundant_indices=redundancy.redundant_indices,
            essential_indices=redundancy.essential_indices,
        )
    conflict = _irreducible_conflict(pieces, vars_) if find_conflict else ()
    return FeasibleSetDiagnosticResult(
        constraints=pieces,
        variables=vars_,
        feasible=False,
        witness=None,
        conflict_indices=conflict,
    )


__all__ = [
    "ConstraintRedundancyResult",
    "FeasibleSetDiagnosticResult",
    "analyze_constraint_redundancy",
    "diagnose_feasible_set",
]
