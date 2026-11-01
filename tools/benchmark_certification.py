"""Record cold/warm timings and exact work counters for fixed mathematical inputs."""

from __future__ import annotations

import argparse
import json
import statistics
import time
from pathlib import Path

import sympy as sp


def workload(name):
    """Return a checked callable; the JSON manifest contains data, never executable code."""
    from semialg import function_range, quantifier_eliminate
    from semialg.algebraic import isolate_real_roots
    from semialg.algebraic.samples import FiberRootContext
    from semialg.cad_algorithms.decomposition import _roots_in_parent_tower
    from semialg.simplify.implication import _is_unsatisfiable_cached

    x, y = sp.symbols("x y", real=True)
    if name == "closed-qe":

        def run():
            value = quantifier_eliminate(sp.And(sp.Eq(x * x - 2, 0), x > 0), (("exists", x),))
            if value is not sp.true:
                raise AssertionError(value)
            return value
    elif name == "lazy-qe":

        def run():
            from semialg.formula import parse_formula
            from semialg.partial.qe import lazy_resolve_formula

            result = lazy_resolve_formula(
                (x, y), (("forall", x), ("exists", y)), parse_formula(sp.Eq(y * y, x))
            )
            if result.truth_value is not False:
                raise AssertionError(result)
            return {
                "truth": result.truth_value,
                "lifted_stacks": result.stats.lifted_stacks,
                "evaluated_leaf_cells": result.stats.evaluated_leaf_cells,
                "visited_cells_by_level": dict(result.stats.visited_cells_by_level),
            }
    elif name == "free-parameter-qe":

        def run():
            a = sp.Symbol("a", real=True)
            result = quantifier_eliminate(
                sp.Eq(y * y, x * x + a),
                (("forall", x), ("exists", y)),
                strategy="partial",
                return_result=True,
            )
            if result.formula != (a >= 0):
                raise AssertionError(result.formula)
            stats = result.backend_result.stats
            return {
                "visited_cells_by_level": dict(stats.visited_cells_by_level),
                "evaluated_leaf_cells": stats.evaluated_leaf_cells,
                "lifted_stacks": stats.lifted_stacks,
            }
    elif name == "linear-boolean":

        def run():
            a, b = sp.symbols("a b", real=True)
            matrix = sp.Or(sp.And(y > x + a, y < x + b), sp.And(y > x + b, y < x + a))
            result = quantifier_eliminate(matrix, (("forall", x), ("exists", y)))
            for av, bv in ((-1, 0), (0, 0), (2, -1)):
                if bool(result.subs({a: av, b: bv})) != (av != bv):
                    raise AssertionError(result)
            return result
    elif name == "native-multiplicity":

        def run():
            parent = isolate_real_roots(sp.Poly(x * x - 2, x))[1]
            roots = _roots_in_parent_tower(
                FiberRootContext(sp.Poly(y * y - x * x + 2, x, y), (x, y), (parent,), y)
            )
            if len(roots) != 1 or roots[0].multiplicity != 2:
                raise AssertionError(roots)
            return roots
    elif name == "shared-image":

        def run():
            t = sp.Symbol("t", real=True)
            value = function_range(x * x, sp.And(x >= -1, x <= 2), (x,), value_symbol=t)
            if sp.simplify_logic(sp.Equivalent(value, sp.And(t >= 0, t <= 4))) is not sp.true:
                raise AssertionError(value)
            return value
    elif name == "semantic-implication":

        def run():
            value = _is_unsatisfiable_cached(sp.And(x > 1, x <= 0), (x,))
            if not value:
                raise AssertionError(value)
            return value
    else:
        raise ValueError(f"unknown benchmark: {name}")
    return run


def measure(name, repeats):
    """Check every result, keeping cold outputs alive through the corresponding warm call."""
    from semialg.cache_control import cache_report, clear_caches

    run = workload(name)
    pairs = []
    for _ in range(repeats):
        clear_caches(include_sympy=True, collect=False)
        before = cache_report()
        started = time.perf_counter()
        cold_result = run()
        cold = time.perf_counter() - started
        after_cold = cache_report()
        started = time.perf_counter()
        warm_result = run()
        warm = time.perf_counter() - started
        after_warm = cache_report()
        if cold_result != warm_result:
            raise AssertionError("cache state changed the mathematical result")
        pairs.append(
            {
                "cold_seconds": cold,
                "work": cold_result if isinstance(cold_result, dict) else {},
                "warm_seconds": warm,
                "before": before,
                "after_cold": after_cold,
                "after_warm": after_warm,
            }
        )
    return {
        "id": name,
        "cold_median": statistics.median(p["cold_seconds"] for p in pairs),
        "warm_median": statistics.median(p["warm_seconds"] for p in pairs),
        "runs": pairs,
    }


def main():
    from quality_metadata import measurement_identity

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=Path("certification-benchmarks.json"))
    parser.add_argument("--repeats", type=int, default=3)
    args = parser.parse_args()
    if args.repeats < 1:
        parser.error("--repeats must be positive")
    root = Path(__file__).resolve().parents[1]
    manifest = json.loads((root / "tests/data/performance_cases.json").read_text())
    result = {
        "environment": measurement_identity(),
        "inputs": manifest,
        "results": [measure(case["id"], args.repeats) for case in manifest],
    }
    args.output.write_text(json.dumps(result, indent=2) + "\n")


if __name__ == "__main__":
    main()
