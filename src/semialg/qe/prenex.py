from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass

import sympy as sp

from ..formula import Formula, ParsedPrenexFormula, parse_quant_form_text, to_sympy
from ..model import QEConfig
from .blocks import QuantifierBlock, blocks_to_quantifiers, norm_quant_blocks
from .complete import CompleteQEResult, qe_by_complete_cad


@dataclass(frozen=True)
class PrenexQERequest:
    vars: tuple[sp.Symbol, ...]
    quantifiers: tuple[tuple[str, sp.Symbol], ...]
    matrix: Formula


def eval_quantifier_free(formula: Formula, point: Mapping[sp.Symbol, sp.Expr]) -> bool:
    value = sp.simplify(to_sympy(formula).subs(dict(point)))
    if value is sp.true or value == sp.true:
        return True
    if value is sp.false or value == sp.false:
        return False
    if getattr(value, "free_symbols", None):
        raise ValueError(f"formula still has free symbols after substitution: {value!r}")
    return bool(value)


def qe_prenex(
    vars_: Sequence[sp.Symbol],
    quantifiers: Sequence[tuple[str, sp.Symbol]],
    matrix: Formula,
    config: QEConfig | None = None,
) -> CompleteQEResult:
    # Configuration fields let reduced and partial CAD drivers share this call
    # signature while preserving complete-CAD semantics for prenex QE.
    _ = config
    return qe_by_complete_cad(tuple(vars_), tuple(quantifiers), matrix)


def qe_prenex_suffix(
    vars_: Sequence[sp.Symbol],
    quantifiers: Sequence[tuple[str, sp.Symbol]],
    matrix: Formula,
    *,
    config: QEConfig | None = None,
    **_ignored: object,
) -> CompleteQEResult:
    return qe_prenex(vars_, quantifiers, matrix, config=config)


def qe_blocks(
    vars_: Sequence[sp.Symbol],
    quantifier_blocks: Sequence[QuantifierBlock],
    matrix: Formula,
    config: QEConfig | None = None,
) -> CompleteQEResult:
    quantifiers = blocks_to_quantifiers(norm_quant_blocks(quantifier_blocks))
    return qe_prenex(vars_=vars_, quantifiers=quantifiers, matrix=matrix, config=config)


def qe_parsed(parsed: ParsedPrenexFormula, config: QEConfig | None = None) -> CompleteQEResult:
    return qe_prenex(
        vars_=parsed.vars, quantifiers=parsed.quantifiers, matrix=parsed.matrix, config=config
    )


def qe_text(
    text: str, *, symbols=None, variable_order=None, config: QEConfig | None = None
) -> CompleteQEResult:
    parsed = parse_quant_form_text(text, symbols=symbols, variable_order=variable_order)
    return qe_parsed(parsed, config=config)


__all__ = [
    "PrenexQERequest",
    "eval_quantifier_free",
    "qe_prenex",
    "qe_prenex_suffix",
    "qe_blocks",
    "qe_parsed",
    "qe_text",
]
