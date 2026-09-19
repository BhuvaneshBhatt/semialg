"""Contracts for coverage measurement, completeness, and subsystem gates."""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "check_coverage.py"


def _checker():
    spec = importlib.util.spec_from_file_location("semialg_coverage_checker", SCRIPT)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def _repo(tmp_path: Path) -> Path:
    root = tmp_path / "repo"
    package = root / "src" / "semialg"
    package.mkdir(parents=True)
    (package / "a.py").write_text("a = 1\n", encoding="utf-8")
    (package / "b.py").write_text("b = 2\n", encoding="utf-8")
    return root


def _policy(**overall):
    return {
        "schema_version": 1,
        "source_root": "src/semialg",
        "require_branches": True,
        "overall": overall,
        "subsystems": {
            "left": {"files": ["src/semialg/a.py"]},
            "right": {"files": ["src/semialg/b.py"]},
        },
    }


def _record(statements: int, covered: int, branches: int, covered_branches: int):
    return {
        "summary": {
            "num_statements": statements,
            "covered_lines": covered,
            "num_branches": branches,
            "covered_branches": covered_branches,
        }
    }


def _report(a=None, b=None):
    return {
        "files": {
            "src/semialg/a.py": a or _record(10, 8, 4, 2),
            "src/semialg/b.py": b or _record(10, 10, 4, 4),
        }
    }


def test_metrics_keep_statement_branch_and_combined_percentages_distinct(tmp_path: Path) -> None:
    checker = _checker()
    result = checker.evaluate_report(_report(), _policy(), repo_root=_repo(tmp_path))

    metrics = result["metrics"]
    assert metrics.statement_percent == pytest.approx(90.0)
    assert metrics.branch_percent == pytest.approx(75.0)
    assert metrics.combined_percent == pytest.approx(24 / 28 * 100)


def test_overall_combined_floor_is_enforced_without_conflating_branch_coverage(
    tmp_path: Path,
) -> None:
    checker = _checker()
    result = checker.evaluate_report(
        _report(), _policy(floor_combined=90.0), repo_root=_repo(tmp_path)
    )

    assert len(result["failures"]) == 1
    assert "combined coverage" in result["failures"][0]
    assert "90.00%" in result["failures"][0]


def test_subsystem_floor_catches_small_weak_area_hidden_by_strong_neighbor(tmp_path: Path) -> None:
    checker = _checker()
    policy = _policy()
    policy["subsystems"]["left"]["floor_combined"] = 80.0
    report = _report(a=_record(10, 5, 4, 1), b=_record(90, 90, 20, 20))

    result = checker.evaluate_report(report, policy, repo_root=_repo(tmp_path))

    assert any(failure.startswith("left:") for failure in result["failures"])


def test_incomplete_report_is_rejected_instead_of_inflating_coverage(tmp_path: Path) -> None:
    checker = _checker()
    report = {"files": {"src/semialg/a.py": _record(10, 10, 2, 2)}}

    with pytest.raises(ValueError, match="incomplete"):
        checker.evaluate_report(report, _policy(), repo_root=_repo(tmp_path))


def test_unassigned_source_file_is_rejected(tmp_path: Path) -> None:
    checker = _checker()
    policy = _policy()
    del policy["subsystems"]["right"]

    with pytest.raises(ValueError, match="unassigned"):
        checker.evaluate_report(_report(), policy, repo_root=_repo(tmp_path))


def test_overlapping_subsystem_rules_are_rejected(tmp_path: Path) -> None:
    checker = _checker()
    policy = _policy()
    policy["subsystems"]["right"]["files"].append("src/semialg/a.py")

    with pytest.raises(ValueError, match="overlap"):
        checker.evaluate_report(_report(), policy, repo_root=_repo(tmp_path))


def test_branch_measurement_is_required_when_policy_requests_it(tmp_path: Path) -> None:
    checker = _checker()
    branchless = _report(
        a=_record(10, 10, 0, 0),
        b=_record(10, 10, 0, 0),
    )

    with pytest.raises(ValueError, match="no branch measurements"):
        checker.evaluate_report(branchless, _policy(), repo_root=_repo(tmp_path))


def test_invalid_coverage_counts_are_rejected(tmp_path: Path) -> None:
    checker = _checker()
    bad = _report(a=_record(10, 11, 4, 4))

    with pytest.raises(ValueError, match="covered statements exceed"):
        checker.evaluate_report(bad, _policy(), repo_root=_repo(tmp_path))


def test_zero_branch_subsystem_is_valid_when_repository_has_branch_data(tmp_path: Path) -> None:
    checker = _checker()
    report = _report(a=_record(10, 10, 0, 0), b=_record(10, 10, 4, 4))

    result = checker.evaluate_report(report, _policy(), repo_root=_repo(tmp_path))

    assert result["subsystems"]["left"]["metrics"].branch_percent == 100.0
    assert not result["failures"]


def test_markdown_summary_surfaces_every_subsystem_and_next_target(tmp_path: Path) -> None:
    checker = _checker()
    policy = _policy(next_combined_target=76.0)
    policy["subsystems"]["left"]["target_combined"] = 85.0
    result = checker.evaluate_report(_report(), policy, repo_root=_repo(tmp_path))

    markdown = checker.render_markdown(result, policy)

    assert "| Overall |" in markdown
    assert "`left`" in markdown
    assert "`right`" in markdown
    assert "76.00%" in markdown
    assert "85.00%" in markdown


def test_repository_policy_assigns_every_semialg_source_file_exactly_once() -> None:
    checker = _checker()
    policy = checker.load_policy(ROOT / "tests" / "coverage_policy.toml")
    files = checker._source_files(ROOT, policy["source_root"])
    assignments, overlaps = checker._assign_subsystems(files, policy)

    assert not overlaps
    assigned = set().union(*assignments.values())
    assert assigned == files
    assert len(assignments) >= 10


def test_repository_policy_preserves_current_gate_and_review_targets() -> None:
    checker = _checker()
    policy = checker.load_policy(ROOT / "tests" / "coverage_policy.toml")

    assert policy["overall"]["floor_combined"] == 74.5
    assert policy["overall"]["next_combined_target"] == 76.0
    decision = policy["subsystems"]["decision_parameters_roots"]
    algebraic = policy["subsystems"]["algebraic_certification"]
    geometry = policy["subsystems"]["geometry_topology"]
    assert decision["floor_combined"] == decision["target_combined"] == 85.0
    assert algebraic["floor_combined"] == algebraic["target_combined"] == 80.0
    assert geometry["floor_combined"] == geometry["target_combined"] == 80.0
