from __future__ import annotations

import json
from pathlib import Path
from time import perf_counter

import sympy as sp

from semialg.cad_algorithms.decomposition import decomp_collins_complete
from semialg.formula import parse_formula
from semialg.partial.qe import lazy_find_inst_form
from semialg.sdp_backends import external_sdp_search
from semialg.sos_certificates import plan_sos_search, verify_psd_exact


def timed(fn):
    start = perf_counter()
    value = fn()
    return value, perf_counter() - start


def main() -> None:
    x, y, z, w = sp.symbols("x y z w")
    sparse_poly = (x**2 + y**2 + z**2 + w**2) ** 3
    plan, plan_seconds = timed(lambda: plan_sos_search(sparse_poly, (x, y, z, w)))

    matrix = sp.diag(*([sp.Integer(1)] * 18 + [sp.Integer(0)] * 6))
    psd, psd_seconds = timed(lambda: verify_psd_exact(matrix))

    formula = sp.And(x**2 < 1, y**2 < 1)
    lazy, lazy_seconds = timed(lambda: lazy_find_inst_form((x, y), parse_formula(formula)))
    full, full_seconds = timed(
        lambda: decomp_collins_complete((x**2 - 1, y**2 - 1), lazy.stats.variables)
    )

    external, external_seconds = timed(
        lambda: external_sdp_search(x**4 + y**4, (x, y), backend="clarabel")
    )

    payload = {
        "sos": {
            "dense_gram_dimension": plan.dense_gram_dimension,
            "sparse_gram_dimension": plan.gram_dimension,
            "gram_variables": plan.gram_variables,
            "basis_method": plan.basis_method,
            "planner_seconds": plan_seconds,
        },
        "psd": {
            "dimension": matrix.rows,
            "verified": psd.verified,
            "method": psd.method,
            "seconds": psd_seconds,
        },
        "partial_cad": {
            "lazy_seconds": lazy_seconds,
            "full_seconds": full_seconds,
            "lazy_leaf_cells": lazy.stats.evaluated_leaf_cells,
            "full_leaf_cells": len(full.cells_by_level[2]),
            "stopped_early": lazy.stats.stopped_early,
            "variable_order": [str(v) for v in lazy.stats.variables],
        },
        "external_sdp": {
            "backend": external.backend,
            "status": external.status,
            "success": external.success,
            "seconds": external_seconds,
        },
    }
    output = Path("semialg_p3_benchmark.json")
    output.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(output)


if __name__ == "__main__":
    main()
