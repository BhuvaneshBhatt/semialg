"""Deterministic benchmarks for cache- and certificate-sensitive paths.

Run with ``PYTHONPATH=src python scripts/benchmark_package.py``.  The script
prints JSON so CI or external profiling jobs can archive comparable records.
Wall-clock values are informational; correctness is checked by exact replay or
by direct-vs-modular differential equality before a result is emitted.
"""

from __future__ import annotations

import json
from time import perf_counter

import sympy as sp

import semialg
from semialg.algebraic.groebner_utils import compute_groebner_basis
from semialg.algebraic_decomposition import primary_decomposition
from semialg.cache_control import cache_report, clear_caches


def _time(call):
    start = perf_counter()
    value = call()
    return value, perf_counter() - start


def main() -> None:
    x, y, z = sp.symbols("x y z")
    records: list[dict[str, object]] = []

    clear_caches(include_sympy=True)
    gtz_case = ((x + y) ** 2, y * (x + y))
    cold, cold_s = _time(lambda: primary_decomposition(gtz_case, (x, y)))
    warm, warm_s = _time(lambda: primary_decomposition(gtz_case, (x, y)))
    assert cold.complete and warm.complete
    assert semialg.replay_certificate(cold).verified
    assert semialg.replay_certificate(warm).verified
    records.append(
        {
            "name": "gtz_embedded_nonmonomial",
            "cold_seconds": cold_s,
            "warm_seconds": warm_s,
            "components": len(cold.components),
        }
    )

    generators = (x * y - z, x * z - y, y**2 - z**2)
    direct, direct_s = _time(
        lambda: compute_groebner_basis(
            generators, (x, y, z), order="grevlex", domain=sp.QQ, modular=False
        )
    )
    modular, modular_s = _time(
        lambda: compute_groebner_basis(
            generators, (x, y, z), order="grevlex", domain=sp.QQ, modular=True
        )
    )
    assert tuple(p.as_expr() for p in direct.polys) == tuple(p.as_expr() for p in modular.polys)
    records.append(
        {
            "name": "groebner_direct_vs_modular",
            "direct_seconds": direct_s,
            "modular_seconds": modular_s,
            "basis_size": len(direct.polys),
        }
    )

    formula = sp.And(x**2 + y**2 <= 4, x + y >= 1, x - y <= 2)
    sat, cad_s = _time(lambda: semialg.is_satisfiable(formula))
    assert sat is True
    records.append({"name": "cad_feasibility_2d", "seconds": cad_s, "satisfiable": True})

    cache_before_clear = cache_report()
    clear_caches(include_sympy=True)
    cache_after_clear = cache_report()
    assert all(info["currsize"] == 0 for info in cache_after_clear["gtz"].values())
    assert all(size == 0 for size in cache_after_clear["algebraic"]["sizes"].values())
    assert all(size == 0 for size in cache_after_clear["cad"]["sizes"].values())

    print(
        json.dumps(
            {
                "benchmarks": records,
                "cache_before_clear": cache_before_clear,
                "cache_after_clear": cache_after_clear,
            },
            indent=2,
            default=str,
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
