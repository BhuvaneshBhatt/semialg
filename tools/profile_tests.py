"""Profile unchanged pytest cases, preserving failures and per-phase outcomes."""

from __future__ import annotations

import argparse
import cProfile
import hashlib
import json
import os
import subprocess
import sys
import time
from pathlib import Path

import pytest


@pytest.hookimpl(hookwrapper=True)
def pytest_runtest_call(item):
    profiler = cProfile.Profile()
    profiler.enable()
    try:
        yield
    finally:
        profiler.disable()
        folder = Path(os.environ["SEMIALG_PROFILE_OUTPUT"])
        name = hashlib.sha256(item.nodeid.encode()).hexdigest()[:16]
        profiler.dump_stats(str(folder / (name + ".prof")))


def pytest_runtest_logreport(report):
    row = {
        "nodeid": report.nodeid,
        "phase": report.when,
        "outcome": report.outcome,
        "seconds": report.duration,
        "failure": str(report.longrepr) if report.failed else None,
        "profile": hashlib.sha256(report.nodeid.encode()).hexdigest()[:16] + ".prof",
    }
    with (Path(os.environ["SEMIALG_PROFILE_OUTPUT"]) / "phases.jsonl").open("a") as stream:
        stream.write(json.dumps(row) + "\n")


def main():
    from quality_metadata import measurement_identity

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=Path("build/profile"))
    parser.add_argument("--timeout", type=float, default=1200)
    parser.add_argument("--seed", type=int, default=28419)
    parser.add_argument("pytest_args", nargs=argparse.REMAINDER)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    root = Path(__file__).resolve().parents[1]
    pytest_args = args.pytest_args or ["-m", "slow", "tests"]
    if pytest_args[0] == "--":
        pytest_args = pytest_args[1:]
    command = [
        sys.executable,
        "-m",
        "pytest",
        "-p",
        "profile_tests",
        "-o",
        "addopts=",
        "-o",
        "pythonpath=",
        f"--hypothesis-seed={args.seed}",
        *pytest_args,
    ]
    env = {
        **os.environ,
        "SEMIALG_PROFILE_OUTPUT": str(args.output.resolve()),
        "PYTHONPATH": os.pathsep.join([str(root / "tools"), *sys.path]),
    }
    (args.output / "phases.jsonl").write_text("")
    started = time.perf_counter()
    timed_out = False
    with (args.output / "pytest.log").open("w") as stream:
        try:
            result = subprocess.run(
                command,
                env=env,
                stdout=stream,
                stderr=subprocess.STDOUT,
                timeout=args.timeout,
                check=False,
            )
            code = result.returncode
        except subprocess.TimeoutExpired:
            timed_out = True
            code = 124
    metadata = {
        **measurement_identity(),
        "seed": args.seed,
        "command": command,
        "seconds": time.perf_counter() - started,
        "returncode": code,
        "timed_out": timed_out,
    }
    (args.output / "run.json").write_text(json.dumps(metadata, indent=2) + "\n")
    return code


if __name__ == "__main__":
    raise SystemExit(main())
