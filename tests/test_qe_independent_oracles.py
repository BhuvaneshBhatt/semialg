"""Expected truth comes from known factor roots, not a second solver call."""

import importlib.util
import sys
from pathlib import Path

import pytest

_spec = importlib.util.spec_from_file_location(
    "qe_oracles", Path(__file__).parents[1] / "tools/qe_oracles.py"
)
_oracles = importlib.util.module_from_spec(_spec)
sys.modules[_spec.name] = _oracles
_spec.loader.exec_module(_oracles)


@pytest.mark.parametrize("case", tuple(_oracles.make_cases()), ids=lambda case: str(case.index))
def test_sign_cell_oracle(case):
    rows = _oracles.check_case(case)
    assert len(rows) == 2
    assert all(len(row["results"]) == 4 for row in rows)
