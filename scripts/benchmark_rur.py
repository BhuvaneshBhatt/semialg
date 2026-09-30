"""Benchmark the public RUR constructor with time and peak-memory measurements."""

from __future__ import annotations

import argparse
import gc
import json
import resource
import statistics
import time
import tracemalloc

from semialg.algebraic.rational_univariate import compute_rur
from semialg.benchmarks.rur_corpus import rur_benchmark_cases


def run_case(case, repeats: int):
    samples = []
    for _ in range(repeats):
        gc.collect()
        tracemalloc.start()
        start = time.perf_counter()
        result = compute_rur(case.equations, case.variables)
        elapsed = time.perf_counter() - start
        _current, peak = tracemalloc.get_traced_memory()
        tracemalloc.stop()
        samples.append((elapsed, peak, result.quotient_dimension))
    return {
        "name": case.name,
        "median_seconds": statistics.median(x[0] for x in samples),
        "min_seconds": min(x[0] for x in samples),
        "peak_tracemalloc_bytes": max(x[1] for x in samples),
        "peak_rss_kib": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
        "quotient_dimension": samples[-1][2],
        "repeats": repeats,
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--repeats", type=int, default=3)
    parser.add_argument("--case")
    args = parser.parse_args()
    cases = rur_benchmark_cases()
    if args.case:
        cases = tuple(c for c in cases if c.name == args.case)
    print(json.dumps([run_case(c, args.repeats) for c in cases], indent=2))


if __name__ == "__main__":
    main()
