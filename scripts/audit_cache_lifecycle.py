"""Exercise repeated exact workloads and report cache/heap lifecycle behavior."""

from __future__ import annotations

import gc
import json
import tracemalloc

import sympy as sp

import semialg
from semialg.algebraic_decomposition import primary_decomposition
from semialg.cache_control import cache_report, clear_caches


def main(rounds: int = 30) -> None:
    x, y = sp.symbols("x y")
    clear_caches(include_sympy=True)
    gc.collect()
    tracemalloc.start()
    start_current, start_peak = tracemalloc.get_traced_memory()

    for k in range(rounds):
        f = sp.expand((x + (k % 7) * y + 1) * (x - y + (k % 5)))
        result = primary_decomposition((f,), (x, y))
        if not result.complete or not semialg.replay_certificate(result).verified:
            raise RuntimeError(f"workload failed at round {k}")

    loaded_current, loaded_peak = tracemalloc.get_traced_memory()
    before = cache_report()
    clear_caches(include_sympy=True)
    gc.collect()
    cleared_current, cleared_peak = tracemalloc.get_traced_memory()
    after = cache_report()
    tracemalloc.stop()

    print(
        json.dumps(
            {
                "rounds": rounds,
                "bytes": {
                    "start_current": start_current,
                    "start_peak": start_peak,
                    "loaded_current": loaded_current,
                    "loaded_peak": loaded_peak,
                    "cleared_current": cleared_current,
                    "cleared_peak": cleared_peak,
                },
                "cache_before_clear": before,
                "cache_after_clear": after,
            },
            indent=2,
            default=str,
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
