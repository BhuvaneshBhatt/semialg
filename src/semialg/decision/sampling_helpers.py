from __future__ import annotations

from collections.abc import Iterable, Mapping, Sequence

import sympy as sp

from ..instances.real_fallbacks import satisfies_formula
from ..sampling import sample_points


def _normalize_sample_request(
    count: int, samples: int | str | None, sample_mode: str | None
) -> tuple[int, str]:
    """Normalize sampling controls."""

    resolved_count = int(count)
    resolved_mode = sample_mode or "auto"
    if isinstance(samples, int):
        resolved_count = int(samples)
    elif isinstance(samples, str):
        resolved_mode = samples
    if resolved_count < 0:
        raise ValueError("sample count must be nonnegative")
    key = resolved_mode.lower().replace("-", "_")
    aliases = {
        "component": "per_component",
        "components": "per_component",
        "per_components": "per_component",
        "cell": "per_cell",
        "cells": "per_cell",
        "per_cells": "per_cell",
        "point": "auto",
        "points": "auto",
        "sample": "auto",
        "samples": "auto",
    }
    key = aliases.get(key, key)
    if key not in {"auto", "per_component", "per_cell"}:
        raise ValueError(f"unsupported sample mode: {resolved_mode!r}")
    return resolved_count, key


def _dedupe_samples(
    samples: Iterable[Mapping[sp.Symbol, sp.Expr]],
) -> tuple[Mapping[sp.Symbol, sp.Expr], ...]:
    seen: set[tuple[tuple[str, str], ...]] = set()
    out: list[Mapping[sp.Symbol, sp.Expr]] = []
    for sample in samples:
        key = tuple(sorted((sym.name, sp.sstr(value)) for sym, value in sample.items()))
        if key not in seen:
            seen.add(key)
            out.append(dict(sample))
    return tuple(out)


def _samples_from_components(
    components: Sequence[object],
) -> tuple[Mapping[sp.Symbol, sp.Expr], ...]:
    out: list[Mapping[sp.Symbol, sp.Expr]] = []
    for component in components:
        variable = getattr(component, "variable", None)
        sample_fn = getattr(component, "sample_point", None)
        if sample_fn is None:
            continue
        sample = sample_fn()
        if isinstance(sample, Mapping):
            out.append({sym: sp.simplify(value) for sym, value in sample.items()})
        elif variable is not None:
            out.append({variable: sp.simplify(sample)})
    return _dedupe_samples(out)


def _samples_from_cells(cells: Sequence[object]) -> tuple[Mapping[sp.Symbol, sp.Expr], ...]:
    out: list[Mapping[sp.Symbol, sp.Expr]] = []
    for cell in cells:
        sample_fn = getattr(cell, "sample_point", None)
        if sample_fn is None:
            continue
        sample = sample_fn()
        if isinstance(sample, Mapping):
            out.append({sym: sp.simplify(value) for sym, value in sample.items()})
    return _dedupe_samples(out)


def _collect_structural_samples(
    formula: sp.Expr,
    original_formula: sp.Expr,
    variables: tuple[sp.Symbol, ...],
    metadata: Mapping[str, object],
    *,
    count: int,
    mode: str,
    strategy: str | None,
) -> tuple[Mapping[sp.Symbol, sp.Expr], ...]:
    """Return improved solve samples from exact components/cells when available."""

    if count == 0 and mode == "auto":
        return ()

    connectivity = metadata.get("connectivity")
    components = tuple(metadata.get("components", ()) or ())
    cells = tuple(metadata.get("cells", ()) or ())

    strategy_key = (strategy or "auto").lower().replace("-", "_")
    force_sampling_strategy = strategy_key in {"rational", "grid", "random", "representative"}

    structural: tuple[Mapping[sp.Symbol, sp.Expr], ...] = ()
    if force_sampling_strategy and mode == "auto":
        structural = ()
    elif mode == "per_component":
        if connectivity is not None and getattr(connectivity, "components", None):
            structural = _samples_from_components(tuple(connectivity.components))
        if not structural:
            structural = _samples_from_components(components)
        if not structural and cells:
            structural = _samples_from_cells(cells)
    elif mode == "per_cell":
        structural = _samples_from_cells(cells)
        if not structural and components:
            structural = _samples_from_components(components)
    elif components:
        structural = _samples_from_components(components)
    elif cells:
        structural = _samples_from_cells(cells)

    if structural:
        if mode in {"per_component", "per_cell"}:
            return structural
        return structural[:count]

    if count <= 0:
        return ()

    generated = tuple(
        point
        for point in sample_points(formula, variables, count=count, strategy=strategy or "auto")
        if satisfies_formula(formula, point, strict=False)
    )
    if not generated:
        generated = tuple(
            point
            for point in sample_points(
                original_formula, variables, count=count, strategy=strategy or "auto"
            )
            if satisfies_formula(original_formula, point, strict=False)
        )
    return _dedupe_samples(generated)
