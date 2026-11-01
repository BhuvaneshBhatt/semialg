"""A profiled test failure remains a failing run with usable diagnostic artifacts."""

import importlib.util
import json
from pathlib import Path

_spec = importlib.util.spec_from_file_location(
    "profile_runner", Path(__file__).parents[1] / "tools/profile_tests.py"
)
_runner = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_runner)


def test_profile_failure(tmp_path, monkeypatch):
    case = tmp_path / "test_failure.py"
    case.write_text("def test_failure():\n    assert 1 == 2\n")
    output = tmp_path / "profile"
    monkeypatch.syspath_prepend(str(Path(__file__).parents[1] / "tools"))
    monkeypatch.setattr(
        "sys.argv", ["profile_tests.py", "--output", str(output), "--", "-q", str(case)]
    )
    assert _runner.main() == 1
    report = json.loads((output / "run.json").read_text())
    assert report["returncode"] == 1 and report["timed_out"] is False
    rows = [json.loads(line) for line in (output / "phases.jsonl").read_text().splitlines()]
    (call,) = [row for row in rows if row["phase"] == "call"]
    assert call["outcome"] == "failed" and "assert 1 == 2" in call["failure"]
    assert (output / call["profile"]).is_file()
    assert len(report["source_sha256"]) == 64
