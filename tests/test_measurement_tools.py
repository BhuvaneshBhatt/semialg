"""Measurement tools retain mathematical checks, raw counters and failures."""

import importlib.util
import json
from pathlib import Path

import pytest

ROOT = Path(__file__).parents[1]
_spec = importlib.util.spec_from_file_location(
    "certification_benchmark", ROOT / "tools/benchmark_certification.py"
)
_benchmark = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_benchmark)


@pytest.mark.parametrize(
    "case",
    json.loads((ROOT / "tests/data/performance_cases.json").read_text()),
    ids=lambda case: case["id"],
)
def test_benchmark_inputs(case):
    measured = _benchmark.measure(case["id"], 1)
    assert measured["id"] == case["id"]
    (run,) = measured["runs"]
    assert run["before"]["cad"]["stats"]["complete_cad_hits"] == 0
    assert (
        run["after_warm"]["cad"]["stats"]["complete_cad_hits"]
        >= run["after_cold"]["cad"]["stats"]["complete_cad_hits"]
    )
    assert run["cold_seconds"] >= 0 and run["warm_seconds"] >= 0


def test_unknown_workload():
    with pytest.raises(ValueError, match="unknown benchmark"):
        _benchmark.workload("unregistered")
