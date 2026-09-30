"""Benchmark packed RUR elimination against domain characteristic polynomials."""

from __future__ import annotations

import argparse
import json
import resource
import statistics
import time
import tracemalloc
from functools import partial

from sympy.polys.domains import GF, QQ
from sympy.polys.matrices import DomainMatrix

from semialg.algebraic.rational_univariate.linear_algebra import packed_krylov_coefficients


def _matrix(dimension, domain):
    zero, one = domain.zero, domain.one
    matrix = [[zero for _ in range(dimension)] for _ in range(dimension)]
    for column in range(dimension - 1):
        matrix[column + 1][column] = one
    matrix[0][-1] = -one
    matrix[1][-1] = -one
    # Exact rank-one similarity transform keeps e_0 cyclic while avoiding a
    # benchmark that depends on the literal companion layout.
    vector = [domain.convert((index % 5) + 1) for index in range(dimension)]
    vector[-1] = zero
    transform = [[one if i == j else zero for j in range(dimension)] for i in range(dimension)]
    inverse = [[one if i == j else zero for j in range(dimension)] for i in range(dimension)]
    for i in range(dimension - 1):
        transform[i][-1] += vector[i]
        inverse[i][-1] -= vector[i]

    def multiply(left, right):
        return [
            [
                sum(
                    (left[i][k] * right[k][j] for k in range(dimension)),
                    zero,
                )
                for j in range(dimension)
            ]
            for i in range(dimension)
        ]

    return multiply(multiply(transform, matrix), inverse)


def _measure(function, repeats):
    samples = []
    peak_python = 0
    for _ in range(repeats):
        tracemalloc.start()
        start = time.perf_counter()
        function()
        samples.append(time.perf_counter() - start)
        _, peak = tracemalloc.get_traced_memory()
        tracemalloc.stop()
        peak_python = max(peak_python, peak)
    return {
        "median_seconds": statistics.median(samples),
        "samples_seconds": samples,
        "peak_python_bytes": peak_python,
        "peak_rss_kib": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--repeats", type=int, default=3)
    parser.add_argument("--dimensions", default="16,24,32,48,64,96,128")
    args = parser.parse_args()
    dimensions = tuple(int(item) for item in args.dimensions.split(","))
    report = {}
    for name, domain in (("QQ", QQ), ("GF(65521)", GF(65521))):
        report[name] = {}
        for dimension in dimensions:
            matrix = _matrix(dimension, domain)
            packed = partial(packed_krylov_coefficients, matrix, domain)
            characteristic = DomainMatrix.from_list(matrix, domain).charpoly
            report[name][str(dimension)] = {
                "packed": _measure(packed, args.repeats),
                "charpoly": _measure(characteristic, args.repeats),
            }
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
