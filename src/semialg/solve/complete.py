from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field
from typing import Any

import sympy as sp

from ..formula import ParsedPrenexFormula, parse_quant_form_text
from ..qe.complete import CompleteQEResult, qe_by_complete_cad
from .domains import SolveDomain, normalize_domain


@dataclass
class CompleteSolveResult:
    """Public result object for the complete real-polynomial QE path."""

    method: str
    domain: SolveDomain
    result: sp.Expr
    status: str
    backend: str
    metadata: dict[str, Any] = field(default_factory=dict)
    qe_result: CompleteQEResult | None = None


def _result_metadata(qe: CompleteQEResult) -> dict[str, Any]:
    diag = qe.diagnostics
    return {
        "variables": tuple(map(sp.sstr, qe.variables)),
        "free_variables": tuple(map(sp.sstr, qe.free_variables)),
        "quantified_variables": tuple(map(sp.sstr, qe.quantified_variables)),
        "quantifier_blocks": tuple(
            (block.quantifier, tuple(map(sp.sstr, block.variables)))
            for block in qe.quantifier_blocks
        ),
        "is_sentence": qe.is_sentence,
        "truth_value": qe.truth_value,
        "satisfying_cells": qe.satisfying_cell_indices,
        "cell_count": len(qe.cad.cells),
        "cell_union_cell_count": None if qe.cell_union is None else len(qe.cell_union.cells),
        "witness_samples": qe.witness_samples,
        "variable_reordered": None if diag is None else diag.variable_reordered,
        "diagnostic_notes": () if diag is None else diag.notes,
    }


def reduce_complete_formula(
    parsed: ParsedPrenexFormula,
    *,
    free_variables: Sequence[sp.Symbol] | None = None,
) -> CompleteSolveResult:
    qe = qe_by_complete_cad(
        parsed.vars, parsed.quantifiers, parsed.matrix, free_variables=free_variables
    )
    return CompleteSolveResult(
        method="complete_cad_qe",
        domain=normalize_domain("reals"),
        result=qe.formula,
        status=qe.status,
        backend=qe.backend,
        qe_result=qe,
        metadata=_result_metadata(qe),
    )


def reduce_complete_text(
    text: str,
    *,
    symbols: Mapping[str, sp.Symbol] | None = None,
    variable_order: Sequence[sp.Symbol | str] | None = None,
    free_variables: Sequence[sp.Symbol | str] | None = None,
    return_result: bool = False,
):
    converted_order = None
    if variable_order is not None:
        converted_order = tuple(
            sp.Symbol(v, real=True) if isinstance(v, str) else v for v in variable_order
        )
    converted_free = None
    if free_variables is not None:
        converted_free = tuple(
            sp.Symbol(v, real=True) if isinstance(v, str) else v for v in free_variables
        )
    parsed = parse_quant_form_text(
        text, symbols=dict(symbols or {}), variable_order=converted_order
    )
    solved = reduce_complete_formula(parsed, free_variables=converted_free)
    return solved if return_result else solved.result


def resolve_complete_text(text: str, **kwargs):
    result = reduce_complete_text(text, return_result=True, **kwargs)
    if result.result == sp.true:
        return True
    if result.result == sp.false:
        return False
    return result.result


__all__ = [
    "CompleteSolveResult",
    "reduce_complete_formula",
    "reduce_complete_text",
    "resolve_complete_text",
]
