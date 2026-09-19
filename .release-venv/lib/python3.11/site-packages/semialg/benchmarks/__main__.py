"""Command-line entry point for semialg's benchmark/conformance suite."""

from __future__ import annotations

import argparse

from .cases import cad_cases, groebner_family, tticad_cases
from .runner import run_cad, run_groebner, run_tticad


def _show(result) -> None:
    m = result.metrics
    print(
        f"{result.case.name}: correct={result.correct} wall={m.wall_seconds:.6g}s cpu={m.cpu_seconds:.6g}s"
    )
    for key, value in m.counters.items():
        print(f"  {key}: {value}")


def main() -> None:
    parser = argparse.ArgumentParser(description="Run semialg conformance/performance benchmarks")
    sub = parser.add_subparsers(dest="command", required=True)
    cad_p = sub.add_parser("cad", help="run Wilson CAD Example Bank cases")
    cad_p.add_argument("--case", help="case name or source identifier")
    tti_p = sub.add_parser("tticad", help="run the published TTICAD corpus")
    tti_p.add_argument("--case", help="case name or source identifier")
    gb_p = sub.add_parser("groebner", help="run a parameterized Gröbner family")
    gb_p.add_argument("family", choices=("cyclic", "katsura", "eco", "noon"))
    gb_p.add_argument("n", type=int)
    args = parser.parse_args()

    if args.command == "groebner":
        _show(run_groebner(groebner_family(args.family, args.n)))
        return
    cases = cad_cases() if args.command == "cad" else tticad_cases()
    if args.case:
        cases = tuple(
            c
            for c in cases
            if c.name == args.case or (c.provenance and c.provenance.identifier == args.case)
        )
        if not cases:
            parser.error(f"unknown benchmark case: {args.case}")
    runner = run_cad if args.command == "cad" else run_tticad
    for case in cases:
        _show(runner(case))


if __name__ == "__main__":
    main()
