#!/usr/bin/env python3
"""Deterministic FLINT/Arb benchmark on semialg's actual CAD/root pipelines."""

from __future__ import annotations

import json
import statistics
from time import perf_counter

import sympy as sp

from semialg.algebraic.cache import clear_algebraic_caches
from semialg.algebraic.roots import isolate_real_roots
from semialg.cad_algorithms.performance_cache import clear_cad_caches
from semialg.cad_algorithms.projection.collins import build_collins_proj_set


def _timed(call, repeats: int) -> dict[str, object]:
    values = []
    result = None
    for _ in range(repeats):
        start = perf_counter()
        result = call()
        values.append(perf_counter() - start)
    return {
        "min": min(values),
        "median": statistics.median(values),
        "samples": values,
        "result": result,
    }


def main() -> None:
    x, y, z = sp.symbols("x y z")
    projection_input = (
        sp.Poly(z**3 + (x + y + 1) * z + x * y + 2, x, y, z, domain=sp.QQ),
        sp.Poly(z**2 + (x - y) * z + x**2 + y - 1, x, y, z, domain=sp.QQ),
    )
    clustered = sp.Poly(
        (x**2 - 2) * ((10000 * x - 14143) ** 2 - 2) * ((10000 * x - 14144) ** 2 - 2),
        x,
        domain=sp.QQ,
    )

    def projection() -> object:
        clear_cad_caches()
        return build_collins_proj_set(projection_input, (x, y, z))

    def isolation() -> object:
        clear_algebraic_caches()
        return isolate_real_roots(clustered)

    accelerated_projection = _timed(projection, 3)
    accelerated_isolation = _timed(isolation, 3)

    # Compare against the exact fallback *inside the same production pipeline*.
    import semialg.algebraic.roots as roots_module
    import semialg.cad_algorithms.polynomial_utils as polynomial_utils

    saved_resultant = polynomial_utils.flint_resultant
    saved_eval = roots_module.flint_eval
    saved_gcd = roots_module.flint_gcd
    saved_rational_roots = roots_module.flint_rational_roots
    try:
        polynomial_utils.flint_resultant = lambda *args, **kwargs: None
        roots_module.flint_eval = lambda *args, **kwargs: None
        roots_module.flint_gcd = lambda *args, **kwargs: None
        roots_module.flint_rational_roots = lambda *args, **kwargs: None
        fallback_projection = _timed(projection, 2)
        fallback_isolation = _timed(isolation, 2)
    finally:
        polynomial_utils.flint_resultant = saved_resultant
        roots_module.flint_eval = saved_eval
        roots_module.flint_gcd = saved_gcd
        roots_module.flint_rational_roots = saved_rational_roots

    acc_tower = accelerated_projection.pop("result")
    fall_tower = fallback_projection.pop("result")
    acc_roots = accelerated_isolation.pop("result")
    fall_roots = fallback_isolation.pop("result")
    payload = {
        "projection": {
            "flint": accelerated_projection,
            "sympy_fallback": fallback_projection,
            "same_counts": acc_tower.poly_count_by_level() == fall_tower.poly_count_by_level(),
        },
        "clustered_initial_isolation": {
            "flint_arb": accelerated_isolation,
            "exact_fallback": fallback_isolation,
            "same_root_count": len(acc_roots) == len(fall_roots),
            "multiplicities": [root.multiplicity for root in acc_roots],
        },
    }
    print(json.dumps(payload, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
