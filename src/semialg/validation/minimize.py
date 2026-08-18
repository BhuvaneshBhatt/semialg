from __future__ import annotations

from collections.abc import Callable
from dataclasses import replace

from .corpus import ValidationCase

Predicate = Callable[[ValidationCase], bool]


def minimize_case_by_atoms(case: ValidationCase, still_fails: Predicate) -> ValidationCase:
    """Greedy text-level shrinker for Boolean formulas made of top-level atoms.

    This is deliberately modest: it removes parenthesized conjunct/disjunct
    chunks while preserving any reduced case that still triggers ``still_fails``.
    More sophisticated AST-level delta debugging can build on the same API.
    """

    current = case
    chunks = _split_loose_chunks(case.formula_text)
    changed = True
    while changed and len(chunks) > 1:
        changed = False
        for idx in range(len(chunks)):
            trial_chunks = chunks[:idx] + chunks[idx + 1 :]
            if not trial_chunks:
                continue
            trial_text = " and ".join(trial_chunks)
            trial = replace(
                current,
                formula_text=trial_text,
                metadata={**dict(current.metadata), "minimized_from": case.name},
            )
            if still_fails(trial):
                current = trial
                chunks = trial_chunks
                changed = True
                break
    return current


def _split_loose_chunks(text: str) -> list[str]:
    normalized = text.replace(" OR ", " or ").replace(" AND ", " and ")
    for sep in (" or ", " and "):
        if sep in normalized:
            return [piece.strip() for piece in normalized.split(sep) if piece.strip()]
    return [text.strip()]


__all__ = ["minimize_case_by_atoms"]
