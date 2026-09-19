from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass

import sympy as sp


@dataclass(frozen=True)
class QuantifierBlock:
    quantifier: str
    variables: tuple[sp.Symbol, ...]


def norm_quant_blocks(blocks: Sequence[QuantifierBlock]) -> tuple[QuantifierBlock, ...]:
    normalized = []
    for block in blocks:
        q = block.quantifier.lower()
        if q not in {"exists", "forall"}:
            raise ValueError(f"Unsupported quantifier: {block.quantifier}")
        vars_ = tuple(block.variables)
        if not vars_:
            continue
        if normalized and normalized[-1].quantifier == q:
            normalized[-1] = QuantifierBlock(q, normalized[-1].variables + vars_)
        else:
            normalized.append(QuantifierBlock(q, vars_))
    return tuple(normalized)


def blocks_to_quantifiers(blocks: Sequence[QuantifierBlock]) -> tuple[tuple[str, sp.Symbol], ...]:
    flat: list[tuple[str, sp.Symbol]] = []
    for block in blocks:
        for var in block.variables:
            flat.append((block.quantifier.lower(), var))
    return tuple(flat)


def quantifiers_to_blocks(
    quantifiers: Sequence[tuple[str, sp.Symbol]],
) -> tuple[QuantifierBlock, ...]:
    blocks: list[QuantifierBlock] = []
    for q, var in quantifiers:
        q = q.lower()
        if q not in {"exists", "forall"}:
            raise ValueError(f"Unsupported quantifier: {q}")
        if blocks and blocks[-1].quantifier == q:
            blocks[-1] = QuantifierBlock(q, blocks[-1].variables + (var,))
        else:
            blocks.append(QuantifierBlock(q, (var,)))
    return tuple(blocks)


__all__ = ["QuantifierBlock", "norm_quant_blocks", "blocks_to_quantifiers", "quantifiers_to_blocks"]
