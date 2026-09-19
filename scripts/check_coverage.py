"""Validate repository and subsystem coverage from a coverage.py JSON report."""

from __future__ import annotations

import argparse
import json
from dataclasses import dataclass
from pathlib import Path, PurePosixPath
from typing import Any

try:
    import tomllib
except ModuleNotFoundError:
    import tomli as tomllib

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_POLICY = ROOT / "tests" / "coverage_policy.toml"


@dataclass(frozen=True)
class Metrics:
    statements: int
    covered_statements: int
    branches: int
    covered_branches: int

    @property
    def statement_percent(self) -> float:
        return _percent(self.covered_statements, self.statements)

    @property
    def branch_percent(self) -> float:
        return _percent(self.covered_branches, self.branches)

    @property
    def combined_percent(self) -> float:
        return _percent(
            self.covered_statements + self.covered_branches,
            self.statements + self.branches,
        )

    def __add__(self, other: Metrics) -> Metrics:
        return Metrics(
            self.statements + other.statements,
            self.covered_statements + other.covered_statements,
            self.branches + other.branches,
            self.covered_branches + other.covered_branches,
        )


ZERO_METRICS = Metrics(0, 0, 0, 0)


def _percent(covered: int, total: int) -> float:
    return 100.0 if total == 0 else covered * 100.0 / total


def _metric_from_summary(summary: dict[str, Any]) -> Metrics:
    statements = int(summary.get("num_statements", 0))
    covered_statements = int(summary.get("covered_lines", 0))
    branches = int(summary.get("num_branches", 0))
    covered_branches = int(summary.get("covered_branches", 0))
    if min(statements, covered_statements, branches, covered_branches) < 0:
        raise ValueError("coverage counts must be nonnegative")
    if covered_statements > statements:
        raise ValueError("covered statements exceed the statement denominator")
    if covered_branches > branches:
        raise ValueError("covered branches exceed the branch denominator")
    return Metrics(statements, covered_statements, branches, covered_branches)


def _canonical_source_path(path: str, source_root: str) -> str:
    normalized = path.replace("\\", "/").lstrip("./")
    root = source_root.rstrip("/")
    if normalized == root or normalized.startswith(root + "/"):
        return normalized
    package = root.split("/")[-1]
    marker = package + "/"
    index = normalized.find(marker)
    if index >= 0:
        return root.rsplit("/", 1)[0] + "/" + normalized[index:]
    return normalized


def _source_files(repo_root: Path, source_root: str) -> set[str]:
    root = repo_root / source_root
    if not root.is_dir():
        raise ValueError(f"coverage source root does not exist: {source_root}")
    return {
        path.relative_to(repo_root).as_posix()
        for path in root.rglob("*.py")
        if "__pycache__" not in path.parts
    }


def _matches(entry: dict[str, Any], path: str, source_root: str) -> bool:
    if path in entry.get("files", ()):
        return True
    for directory in entry.get("directories", ()):
        prefix = str(directory).rstrip("/") + "/"
        if path.startswith(prefix):
            return True
    pure = PurePosixPath(path)
    if any(pure.match(pattern) for pattern in entry.get("patterns", ())):
        return True
    if entry.get("top_level", False):
        relative = PurePosixPath(path).relative_to(PurePosixPath(source_root))
        return len(relative.parts) == 1
    return False


def _assign_subsystems(
    files: set[str], policy: dict[str, Any]
) -> tuple[dict[str, set[str]], dict[str, list[str]]]:
    source_root = policy["source_root"]
    subsystem_policy = policy.get("subsystems", {})
    assignments: dict[str, set[str]] = {name: set() for name in subsystem_policy}
    overlaps: dict[str, list[str]] = {}
    for path in sorted(files):
        explicit = [
            name
            for name, entry in subsystem_policy.items()
            if not entry.get("top_level", False) and _matches(entry, path, source_root)
        ]
        if len(explicit) > 1:
            overlaps[path] = explicit
            continue
        if explicit:
            assignments[explicit[0]].add(path)
            continue
        fallback = [
            name
            for name, entry in subsystem_policy.items()
            if entry.get("top_level", False) and _matches(entry, path, source_root)
        ]
        if len(fallback) > 1:
            overlaps[path] = fallback
        elif fallback:
            assignments[fallback[0]].add(path)
    return assignments, overlaps


def _floor_failures(label: str, metrics: Metrics, policy: dict[str, Any]) -> list[str]:
    failures: list[str] = []
    for metric_name, value in (
        ("statement", metrics.statement_percent),
        ("branch", metrics.branch_percent),
        ("combined", metrics.combined_percent),
    ):
        key = f"floor_{metric_name}"
        if key not in policy:
            continue
        floor = float(policy[key])
        if not 0 <= floor <= 100:
            failures.append(f"{label}: invalid {key}={floor:g}")
        elif value + 1e-12 < floor:
            failures.append(f"{label}: {metric_name} coverage {value:.2f}% is below {floor:.2f}%")
    return failures


def evaluate_report(
    report: dict[str, Any], policy: dict[str, Any], *, repo_root: Path = ROOT
) -> dict[str, Any]:
    if int(policy.get("schema_version", 0)) != 1:
        raise ValueError("coverage policy schema_version must be 1")
    source_root = str(policy.get("source_root", "src/semialg"))
    report_files = report.get("files")
    if not isinstance(report_files, dict) or not report_files:
        raise ValueError("coverage report contains no files")

    expected = _source_files(repo_root, source_root)
    canonical: dict[str, dict[str, Any]] = {}
    for raw_path, record in report_files.items():
        path = _canonical_source_path(str(raw_path), source_root)
        if path in canonical:
            raise ValueError(f"coverage report contains duplicate canonical path: {path}")
        canonical[path] = record

    measured_source = set(canonical) & expected
    missing = sorted(expected - measured_source)
    if missing:
        preview = ", ".join(missing[:8])
        suffix = " ..." if len(missing) > 8 else ""
        raise ValueError(
            f"coverage report is incomplete; missing {len(missing)} source files: {preview}{suffix}"
        )

    metrics_by_file = {
        path: _metric_from_summary(canonical[path].get("summary", {})) for path in expected
    }
    total = ZERO_METRICS
    for metrics in metrics_by_file.values():
        total += metrics
    if total.statements == 0:
        raise ValueError("coverage report has a zero statement denominator")
    if bool(policy.get("require_branches", True)) and total.branches == 0:
        raise ValueError("coverage report has no branch measurements")

    assignments, overlaps = _assign_subsystems(expected, policy)
    if overlaps:
        first_path = sorted(overlaps)[0]
        raise ValueError(
            f"coverage subsystem rules overlap for {first_path}: {', '.join(overlaps[first_path])}"
        )
    assigned = set().union(*assignments.values()) if assignments else set()
    unassigned = sorted(expected - assigned)
    if unassigned:
        preview = ", ".join(unassigned[:8])
        suffix = " ..." if len(unassigned) > 8 else ""
        raise ValueError(
            f"coverage subsystem rules leave {len(unassigned)} files unassigned: {preview}{suffix}"
        )

    failures = _floor_failures("overall", total, policy.get("overall", {}))
    subsystem_results: dict[str, dict[str, Any]] = {}
    for name, paths in assignments.items():
        if not paths:
            raise ValueError(f"coverage subsystem {name!r} matches no source files")
        metrics = ZERO_METRICS
        for path in paths:
            metrics += metrics_by_file[path]
        entry = policy["subsystems"][name]
        failures.extend(_floor_failures(name, metrics, entry))
        subsystem_results[name] = {
            "metrics": metrics,
            "files": tuple(sorted(paths)),
            "target_combined": entry.get("target_combined"),
        }

    return {
        "metrics": total,
        "subsystems": subsystem_results,
        "failures": failures,
    }


def render_markdown(result: dict[str, Any], policy: dict[str, Any]) -> str:
    total: Metrics = result["metrics"]
    lines = [
        "# Coverage summary",
        "",
        "| Scope | Statements | Branches | Combined | Floor | Target |",
        "| --- | ---: | ---: | ---: | ---: | ---: |",
    ]
    overall_policy = policy.get("overall", {})
    overall_floor = overall_policy.get("floor_combined")
    overall_target = overall_policy.get("next_combined_target")
    floor_text = "—" if overall_floor is None else f"{float(overall_floor):.2f}%"
    target_text = "—" if overall_target is None else f"{float(overall_target):.2f}%"
    lines.append(
        f"| Overall | {total.statement_percent:.2f}% | {total.branch_percent:.2f}% | "
        f"{total.combined_percent:.2f}% | {floor_text} | {target_text} |"
    )
    for name, entry in result["subsystems"].items():
        metrics: Metrics = entry["metrics"]
        subsystem_policy = policy["subsystems"][name]
        floor = subsystem_policy.get("floor_combined")
        target = entry["target_combined"]
        floor_text = "—" if floor is None else f"{float(floor):.2f}%"
        target_text = "—" if target is None else f"{float(target):.2f}%"
        lines.append(
            f"| `{name}` | {metrics.statement_percent:.2f}% | {metrics.branch_percent:.2f}% | "
            f"{metrics.combined_percent:.2f}% | {floor_text} | {target_text} |"
        )
    lines.extend(
        [
            "",
            f"Measured source files: {sum(len(v['files']) for v in result['subsystems'].values())}.",
        ]
    )
    if result["failures"]:
        lines.extend(["", "## Gate failures", ""])
        lines.extend(f"- {failure}" for failure in result["failures"])
    else:
        lines.extend(["", "All enforced coverage floors pass."])
    lines.append("")
    return "\n".join(lines)


def load_policy(path: Path) -> dict[str, Any]:
    return tomllib.loads(path.read_text(encoding="utf-8"))


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("report", type=Path, help="coverage.py JSON report")
    parser.add_argument("--policy", type=Path, default=DEFAULT_POLICY)
    parser.add_argument("--markdown", type=Path)
    args = parser.parse_args(argv)

    policy = load_policy(args.policy)
    report = json.loads(args.report.read_text(encoding="utf-8"))
    result = evaluate_report(report, policy)
    summary = render_markdown(result, policy)
    print(summary, end="")
    if args.markdown is not None:
        args.markdown.parent.mkdir(parents=True, exist_ok=True)
        args.markdown.write_text(summary, encoding="utf-8")
    return 1 if result["failures"] else 0


if __name__ == "__main__":
    raise SystemExit(main())
