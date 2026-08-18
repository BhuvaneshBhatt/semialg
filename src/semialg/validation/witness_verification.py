from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass

import sympy as sp

from .solution_checking import form_sat_by_assign


@dataclass(frozen=True)
class CandidateWitnessVerdict:
    assignment: dict[sp.Symbol, object]
    is_valid: bool
    details: str = ""


def verify_candidate_witness(
    formula: sp.Expr,
    assignment: Mapping[sp.Symbol, object],
    *,
    domain=sp.Complexes,
    modulus: int | None = None,
    check_numeric_equalities: bool = True,
) -> CandidateWitnessVerdict:
    ok = form_sat_by_assign(
        formula,
        assignment,
        domain=domain,
        modulus=modulus,
        check_numeric_equalities=check_numeric_equalities,
    )
    return CandidateWitnessVerdict(
        dict(assignment), ok, "" if ok else "assignment does not satisfy formula"
    )


def verify_cand_wits(
    formula: sp.Expr,
    assignments: Sequence[Mapping[sp.Symbol, object]],
    *,
    domain=sp.Complexes,
    modulus: int | None = None,
    check_numeric_equalities: bool = True,
) -> list[CandidateWitnessVerdict]:
    return [
        verify_candidate_witness(
            formula,
            a,
            domain=domain,
            modulus=modulus,
            check_numeric_equalities=check_numeric_equalities,
        )
        for a in assignments
    ]


__all__ = ["CandidateWitnessVerdict", "verify_candidate_witness", "verify_cand_wits"]
