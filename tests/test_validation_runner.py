from __future__ import annotations

import sympy as sp

from semialg.validation import (
    RandomFormulaConfig,
    ValidationCase,
    built_in_smoke_cases,
    find_grid_witness,
    random_validation_cases,
    read_jsonl_cases,
    run_validation_cases,
    write_jsonl_cases,
)
from semialg.validation.minimize import minimize_case_by_atoms


def test_validation_run_01():
    report = run_validation_cases(built_in_smoke_cases())
    assert report.passed, report.to_json()
    assert len(report.results) == 3


def test_validation_run_02():
    config = RandomFormulaConfig(seed=17, variables=("x", "y"), atom_count=2, quantifier_count=1)
    first = random_validation_cases(3, config)
    second = random_validation_cases(3, config)
    assert first == second
    assert all(case.tags == ("random",) for case in first)


def test_jsonl_round_trip(tmp_path):
    cases = built_in_smoke_cases()
    path = tmp_path / "cases.jsonl"
    write_jsonl_cases(path, cases)
    assert read_jsonl_cases(path) == cases


def test_validation_run_03():
    x = sp.Symbol("x", real=True)
    check = find_grid_witness(x > 0, x >= 0, (x,))
    assert not check.equivalent_on_grid
    assert check.witness is not None
    assert check.witness.assignment[x] == 0


def test_validation_run_04():
    case = ValidationCase(
        name="bad",
        formula_text="(x > 0) and (x = 0) and (x < 1)",
        variables=("x",),
        tags=("unit",),
    )

    def still_fails(candidate: ValidationCase) -> bool:
        return "x = 0" in candidate.formula_text

    minimized = minimize_case_by_atoms(case, still_fails)
    assert "x = 0" in minimized.formula_text
    assert len(minimized.formula_text) <= len(case.formula_text)
