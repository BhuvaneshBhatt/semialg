"""Small deterministic GTZ benchmark corpus.

Run with ``PYTHONPATH=src python scripts/benchmark_gtz.py``.  Timings are
informational only; every completed case is independently certificate-replayed.
The printed ideals are also convenient inputs for external Singular primdecGTZ
comparisons.
"""

from __future__ import annotations

from time import perf_counter

import sympy as sp

from semialg.algebraic.gtz_primary import (
    clear_gtz_caches,
    gtz_cache_info,
    gtz_primary_decomposition,
    verify_gtz_primary_decomposition_certificate,
)


def main() -> None:
    x, y, z = sp.symbols("x y z")
    cases = (
        ("two_components", (x * y,), (x, y)),
        ("embedded_nonmonomial", ((x + y) ** 2, y * (x + y)), (x, y)),
        ("localized_split", (x * (x - 1),), (x, y)),
        ("space_curve", (x * y - z, x * z - y), (x, y, z)),
    )
    clear_gtz_caches()
    for name, generators, variables in cases:
        start = perf_counter()
        result = gtz_primary_decomposition(generators, variables)
        elapsed = perf_counter() - start
        verified = bool(
            result.complete
            and result.certificate is not None
            and verify_gtz_primary_decomposition_certificate(result.certificate)
        )
        print(
            f"{name}: complete={result.complete} verified={verified} "
            f"components={len(result.components)} seconds={elapsed:.6f}"
        )
    print("cache:", gtz_cache_info())


if __name__ == "__main__":
    main()
