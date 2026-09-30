from pathlib import Path

from tools.run_slow_tests import RunResult, write_reports


def test_reports_preserve_diagnostics(tmp_path: Path) -> None:
    result = RunResult(
        nodeid="tests/test_example.py::test_case",
        status="passed",
        seconds=1.25,
        returncode=0,
        exceeded_soft_timeout=False,
        command=["python", "-m", "pytest"],
        log="case.log",
    )
    write_reports([result], tmp_path)
    report = (tmp_path / "slow-tests.json").read_text(encoding="utf-8")
    markdown = (tmp_path / "slow-tests.md").read_text(encoding="utf-8")
    assert '"passed": 1' in report
    assert '"nodeid": "tests/test_example.py::test_case"' in report
    assert "1.25 | passed | `tests/test_example.py::test_case`" in markdown
