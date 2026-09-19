from __future__ import annotations

from dataclasses import asdict, dataclass

from ..benchmarks import BenchmarkCaseResult, BenchmarkSuiteResult


@dataclass(frozen=True)
class CaseDashboardRow:
    case: str
    passed: bool
    taxonomy: tuple[str, ...]
    run_count: int
    checker_mismatch_count: int
    diff_mismatch_count: int
    total_cells: int
    total_visited_cells: int


def case_to_row(case: BenchmarkCaseResult) -> CaseDashboardRow:
    total_cells = sum(run.metrics.cells_total for run in case.runs.values())
    total_visited = sum(run.metrics.cells_visited for run in case.runs.values())
    return CaseDashboardRow(
        case=case.case.name,
        passed=case.passed,
        taxonomy=tuple(case.taxonomy),
        run_count=len(case.runs),
        checker_mismatch_count=len(case.checker_mismatches),
        diff_mismatch_count=len(case.differential_mismatches),
        total_cells=total_cells,
        total_visited_cells=total_visited,
    )


def suite_dashboard(result: BenchmarkSuiteResult) -> dict[str, object]:
    rows = [case_to_row(case) for case in result.cases]
    return {
        "suite": result.name,
        "passed": result.passed,
        "summary": result.summary(),
        "rows": [asdict(row) for row in rows],
        "totals": {
            "cases": len(rows),
            "runs": sum(row.run_count for row in rows),
            "cells": sum(row.total_cells for row in rows),
            "visited_cells": sum(row.total_visited_cells for row in rows),
            "checker_mismatches": sum(row.checker_mismatch_count for row in rows),
            "differential_mismatches": sum(row.diff_mismatch_count for row in rows),
        },
    }


def format_dashboard(result: BenchmarkSuiteResult) -> str:
    dash = suite_dashboard(result)
    lines = [
        f"suite={dash['suite']} passed={dash['passed']} cases={dash['totals']['cases']} runs={dash['totals']['runs']}"
    ]
    for row in dash["rows"]:
        lines.append(
            f"- {row['case']}: passed={row['passed']} taxonomy={','.join(row['taxonomy']) if row['taxonomy'] else '-'} "
            f"runs={row['run_count']} cells={row['total_cells']} visited={row['total_visited_cells']} "
            f"checker_mismatches={row['checker_mismatch_count']} differential_mismatches={row['diff_mismatch_count']}"
        )
    return "\n".join(lines)
