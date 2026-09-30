"""Run slow pytest cases serially in fresh processes and record diagnostics."""

from __future__ import annotations

import argparse
import json
import os
import shutil
import signal
import subprocess
import sys
import time
from dataclasses import asdict, dataclass
from datetime import UTC, datetime
from pathlib import Path


@dataclass(frozen=True)
class RunResult:
    """Outcome and timing for one isolated pytest node."""

    nodeid: str
    status: str
    seconds: float
    returncode: int | None
    exceeded_soft_timeout: bool
    command: list[str]
    log: str


def _pytest_command(nodeid: str, extra_args: list[str]) -> list[str]:
    return [
        sys.executable,
        "-m",
        "pytest",
        "-o",
        "addopts=",
        "-m",
        "slow",
        nodeid,
        "-ra",
        *extra_args,
    ]


def collect_slow_tests(pytest_args: list[str]) -> list[str]:
    """Collect slow-test node IDs without running test bodies."""
    command = [
        sys.executable,
        "-m",
        "pytest",
        "-o",
        "addopts=",
        "-m",
        "slow",
        "--collect-only",
        "-q",
        *pytest_args,
    ]
    completed = subprocess.run(command, text=True, capture_output=True, check=False)
    if completed.returncode not in (0, 5):
        raise RuntimeError(
            "slow-test collection failed\n"
            f"command: {' '.join(command)}\n{completed.stdout}{completed.stderr}"
        )
    return sorted(
        line.strip()
        for line in completed.stdout.splitlines()
        if "::" in line and not line.lstrip().startswith(("=", "ERROR"))
    )


def run_isolated_test(
    nodeid: str,
    *,
    log_dir: Path,
    soft_timeout: float,
    hard_timeout: float,
    pytest_args: list[str],
) -> RunResult:
    """Run one pytest node in a fresh process and enforce a hard wall timeout."""
    command = _pytest_command(nodeid, pytest_args)
    safe_name = nodeid.replace("/", "__").replace("::", "--").replace("[", "_").replace("]", "_")
    log_path = log_dir / f"{safe_name}.log"
    started = time.monotonic()
    timed_out = False
    with log_path.open("w", encoding="utf-8") as log_file:
        timeout_program = shutil.which("timeout") if os.name == "posix" else None
        if timeout_program:
            bounded_command = [
                timeout_program,
                "--signal=TERM",
                "--kill-after=5",
                str(hard_timeout),
                *command,
            ]
            completed = subprocess.run(
                bounded_command, stdout=log_file, stderr=subprocess.STDOUT, check=False
            )
            returncode = completed.returncode
            timed_out = returncode in (124, 137)
        else:
            process = subprocess.Popen(
                command, stdout=log_file, stderr=subprocess.STDOUT, start_new_session=True
            )
            try:
                process.wait(timeout=hard_timeout)
            except subprocess.TimeoutExpired:
                timed_out = True
                os.killpg(process.pid, signal.SIGTERM)
                try:
                    process.wait(timeout=5)
                except subprocess.TimeoutExpired:
                    os.killpg(process.pid, signal.SIGKILL)
                    process.wait()
            returncode = process.returncode
    seconds = time.monotonic() - started
    status = "timeout" if timed_out else ("passed" if returncode == 0 else "failed")
    return RunResult(
        nodeid=nodeid,
        status=status,
        seconds=seconds,
        returncode=returncode,
        exceeded_soft_timeout=seconds > soft_timeout,
        command=command,
        log=str(log_path),
    )


def write_reports(results: list[RunResult], report_dir: Path) -> None:
    """Write JSON and Markdown reports ordered by execution and by runtime."""
    payload = {
        "generated_at": datetime.now(UTC).isoformat(),
        "summary": {
            "total": len(results),
            "passed": sum(result.status == "passed" for result in results),
            "failed": sum(result.status == "failed" for result in results),
            "timed_out": sum(result.status == "timeout" for result in results),
            "soft_timeout_exceeded": sum(result.exceeded_soft_timeout for result in results),
            "total_seconds": sum(result.seconds for result in results),
        },
        "results": [asdict(result) for result in results],
    }
    (report_dir / "slow-tests.json").write_text(
        json.dumps(payload, indent=2) + "\n", encoding="utf-8"
    )
    lines = [
        "# Isolated slow-test report",
        "",
        f"Tests: {payload['summary']['total']}",
        f"Passed: {payload['summary']['passed']}",
        f"Failed: {payload['summary']['failed']}",
        f"Timed out: {payload['summary']['timed_out']}",
        f"Total wall time: {payload['summary']['total_seconds']:.1f} s",
        "",
        "## Slowest tests",
        "",
        "| Seconds | Status | Test |",
        "| ---: | --- | --- |",
    ]
    for result in sorted(results, key=lambda item: item.seconds, reverse=True):
        lines.append(f"| {result.seconds:.2f} | {result.status} | `{result.nodeid}` |")
    lines.extend(
        [
            "",
            "## Reproduction",
            "",
            "Each JSON result contains the exact pytest command and log path.",
        ]
    )
    (report_dir / "slow-tests.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("paths", nargs="*", default=["tests"], help="pytest collection paths")
    parser.add_argument(
        "--soft-timeout", type=float, default=60.0, help="flag tests slower than this many seconds"
    )
    parser.add_argument(
        "--hard-timeout", type=float, default=300.0, help="terminate a test after this many seconds"
    )
    parser.add_argument("--report-dir", type=Path, default=Path(".test-reports/slow"))
    parser.add_argument("--start-index", type=int, default=0, help="skip this many collected tests")
    parser.add_argument("--max-tests", type=int, help="run at most N tests after --start-index")
    parser.add_argument(
        "--resume", action="store_true", help="skip node IDs already present in the JSON report"
    )
    parser.add_argument(
        "--pytest-arg",
        action="append",
        default=[],
        help="extra argument passed to each pytest process",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    if args.soft_timeout <= 0 or args.hard_timeout <= 0:
        raise SystemExit("timeouts must be positive")
    if args.soft_timeout > args.hard_timeout:
        raise SystemExit("soft timeout cannot exceed hard timeout")
    report_dir = args.report_dir
    log_dir = report_dir / "logs"
    log_dir.mkdir(parents=True, exist_ok=True)
    nodeids = collect_slow_tests(args.paths)
    nodeids = nodeids[args.start_index :]
    if args.max_tests is not None:
        nodeids = nodeids[: args.max_tests]
    results: list[RunResult] = []
    report_path = report_dir / "slow-tests.json"
    if args.resume and report_path.exists():
        previous = json.loads(report_path.read_text(encoding="utf-8"))
        results = [RunResult(**item) for item in previous.get("results", [])]
        completed = {result.nodeid for result in results}
        nodeids = [nodeid for nodeid in nodeids if nodeid not in completed]
    print(f"Scheduled {len(nodeids)} slow tests", flush=True)
    for index, nodeid in enumerate(nodeids, 1):
        print(f"[{index}/{len(nodeids)}] {nodeid}", flush=True)
        result = run_isolated_test(
            nodeid,
            log_dir=log_dir,
            soft_timeout=args.soft_timeout,
            hard_timeout=args.hard_timeout,
            pytest_args=args.pytest_arg,
        )
        results.append(result)
        marker = " SLOW" if result.exceeded_soft_timeout else ""
        print(f"  {result.status} {result.seconds:.2f}s{marker}", flush=True)
        write_reports(results, report_dir)
    return 1 if any(result.status != "passed" for result in results) else 0


if __name__ == "__main__":
    raise SystemExit(main())
