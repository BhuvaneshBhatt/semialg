from __future__ import annotations

import json
import time
from collections.abc import Iterable, Sequence
from dataclasses import dataclass, field
from pathlib import Path

import sympy as sp

from ..parser import parse_quantified_formula
from ..qe.complete import qe_by_complete_cad
from .checkers import CheckResult, FormulaChecker, SymPyInequalityChecker
from .corpus import ValidationCase
from .symmetric_difference import SymmetricDifferenceCheck, find_grid_witness


@dataclass(frozen=True)
class CaseValidationResult:
    case: ValidationCase
    solver_status: str
    solver_formula: str | None
    elapsed_seconds: float
    passed: bool
    checker_results: tuple[CheckResult, ...] = ()
    equivalence_checks: tuple[SymmetricDifferenceCheck, ...] = ()
    diagnostics: tuple[str, ...] = ()

    def to_json_dict(self) -> dict[str, object]:
        return {
            "case": self.case.to_json_dict(),
            "solver_status": self.solver_status,
            "solver_formula": self.solver_formula,
            "elapsed_seconds": self.elapsed_seconds,
            "passed": self.passed,
            "checker_results": [
                {
                    "checker_name": o.checker_name,
                    "available": o.available,
                    "formula": None if o.formula is None else sp.sstr(o.formula),
                    "truth_value": o.truth_value,
                    "status": o.status,
                    "diagnostics": list(o.diagnostics),
                }
                for o in self.checker_results
            ],
            "equivalence_checks": [
                {
                    "equivalent_on_grid": check.equivalent_on_grid,
                    "checked_points": check.checked_points,
                    "witness": None
                    if check.witness is None
                    else {
                        "assignment": {
                            str(k): sp.sstr(v) for k, v in check.witness.assignment.items()
                        },
                        "left_value": check.witness.left_value,
                        "right_value": check.witness.right_value,
                    },
                }
                for check in self.equivalence_checks
            ],
            "diagnostics": list(self.diagnostics),
        }


@dataclass(frozen=True)
class ValidationRunReport:
    results: tuple[CaseValidationResult, ...]
    started_at: float = field(default_factory=time.time)

    @property
    def passed(self) -> bool:
        return all(result.passed for result in self.results)

    def to_json(self) -> str:
        return json.dumps(
            {"passed": self.passed, "results": [r.to_json_dict() for r in self.results]},
            indent=2,
            sort_keys=True,
        )

    def write_json(self, path: str | Path) -> None:
        out = Path(path)
        out.parent.mkdir(parents=True, exist_ok=True)
        out.write_text(self.to_json(), encoding="utf-8")


def validate_case(
    case: ValidationCase, *, checkers: Sequence[FormulaChecker] | None = None
) -> CaseValidationResult:
    start = time.perf_counter()
    diagnostics: list[str] = []
    checker_results: list[CheckResult] = []
    equivalence_checks: list[SymmetricDifferenceCheck] = []
    try:
        parsed = parse_quantified_formula(
            case.formula_text, symbols=case.symbols(), variable_order=case.sympy_variables()
        )
        case_variables = case.sympy_variables()
        case_quantifiers = case.sympy_quantifiers() or parsed.quantifiers
        result = qe_by_complete_cad(case_variables, case_quantifiers, parsed.matrix)
        solver_formula = result.formula
        solver_status = result.status
    except Exception as exc:
        return CaseValidationResult(
            case, "error", None, time.perf_counter() - start, False, diagnostics=(repr(exc),)
        )

    passed = solver_status == "complete"
    if case.expected_text is not None:
        try:
            expected_expr = _parse_expected(case.expected_text, case.symbols())
            check = find_grid_witness(expected_expr, solver_formula, tuple(result.free_variables))
            equivalence_checks.append(check)
            passed = passed and check.equivalent_on_grid
        except Exception as exc:
            passed = False
            diagnostics.append(f"expected comparison failed: {exc!r}")

    selected_checkers = tuple(checkers or (SymPyInequalityChecker(),))
    for checker in selected_checkers:
        check_result = checker.eliminate(parsed.matrix_expr, case_variables, case_quantifiers)
        checker_results.append(check_result)
        if check_result.formula is not None and not case_quantifiers:
            check = find_grid_witness(
                check_result.formula, solver_formula, tuple(result.free_variables)
            )
            equivalence_checks.append(check)
            passed = passed and check.equivalent_on_grid

    return CaseValidationResult(
        case=case,
        solver_status=solver_status,
        solver_formula=sp.sstr(solver_formula),
        elapsed_seconds=time.perf_counter() - start,
        passed=passed,
        checker_results=tuple(checker_results),
        equivalence_checks=tuple(equivalence_checks),
        diagnostics=tuple(diagnostics),
    )


def run_validation_cases(
    cases: Iterable[ValidationCase], *, checkers: Sequence[FormulaChecker] | None = None
) -> ValidationRunReport:
    return ValidationRunReport(tuple(validate_case(case, checkers=checkers) for case in cases))


def _parse_expected(text: str, symbols: dict[str, sp.Symbol]) -> sp.Expr:
    if text.strip().lower() == "true":
        return sp.true
    if text.strip().lower() == "false":
        return sp.false
    parsed = parse_quantified_formula(text, symbols=symbols)
    return parsed.matrix_expr


__all__ = ["CaseValidationResult", "ValidationRunReport", "validate_case", "run_validation_cases"]
