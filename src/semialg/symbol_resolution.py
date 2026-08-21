from __future__ import annotations

from collections.abc import Iterable, Sequence

import sympy as sp


def _symbols_from(objects: Iterable[object]) -> tuple[sp.Symbol, ...]:
    out: list[sp.Symbol] = []
    seen: set[sp.Symbol] = set()
    for obj in objects:
        if obj is None:
            continue
        if isinstance(obj, sp.Symbol):
            symbols = (obj,)
        elif isinstance(obj, sp.Basic):
            symbols = tuple(obj.free_symbols)
        else:
            try:
                expr = sp.sympify(obj)
            except (TypeError, ValueError, sp.SympifyError):
                continue
            symbols = tuple(expr.free_symbols)
        for symbol in symbols:
            if symbol not in seen:
                out.append(symbol)
                seen.add(symbol)
    return tuple(out)


def resolve_symbol(
    value: sp.Symbol | str,
    *,
    context: Iterable[object] = (),
    known_symbols: Iterable[sp.Symbol] = (),
    create_real: bool = True,
) -> sp.Symbol:
    """Resolve a symbol name against symbols already present in a problem.

    SymPy symbols with the same printed name but different assumptions are
    distinct objects.  Public APIs that accept string variable names therefore
    must prefer the exact symbol objects already present in the input formula,
    objective, bounds, or other contextual expressions.  Ambiguous same-name
    matches are rejected instead of guessed.
    """

    if isinstance(value, sp.Symbol):
        return value
    if not isinstance(value, str):
        raise TypeError(f"expected a Symbol or string name, got {type(value).__name__}")

    candidates: list[sp.Symbol] = []
    seen: set[sp.Symbol] = set()
    for symbol in tuple(known_symbols) + _symbols_from(context):
        if symbol.name == value and symbol not in seen:
            candidates.append(symbol)
            seen.add(symbol)
    if len(candidates) > 1:
        raise ValueError(
            f"symbol name {value!r} is ambiguous across symbols with different assumptions"
        )
    if candidates:
        return candidates[0]
    return sp.Symbol(value, real=True) if create_real else sp.Symbol(value)


def normalize_variables(
    variables: Sequence[sp.Symbol | str] | None,
    *,
    context: Iterable[object] = (),
    append_context_symbols: bool = True,
    exclude: Iterable[sp.Symbol] = (),
    create_real: bool = True,
) -> tuple[sp.Symbol, ...]:
    """Normalize a public variable list while preserving contextual symbols.

    Explicit variables keep their requested order.  When
    ``append_context_symbols`` is true, remaining free symbols from ``context``
    are appended deterministically.  Symbols in ``exclude`` are omitted.
    """

    context_symbols = _symbols_from(context)
    excluded = set(exclude)
    out: list[sp.Symbol] = []
    seen: set[sp.Symbol] = set()
    for value in variables or ():
        symbol = resolve_symbol(
            value,
            context=context,
            known_symbols=context_symbols,
            create_real=create_real,
        )
        if symbol not in excluded and symbol not in seen:
            out.append(symbol)
            seen.add(symbol)
    if append_context_symbols:
        for symbol in sorted(context_symbols, key=lambda item: (item.name, sp.srepr(item))):
            if symbol not in excluded and symbol not in seen:
                out.append(symbol)
                seen.add(symbol)
    return tuple(out)


__all__ = ["normalize_variables", "resolve_symbol"]
