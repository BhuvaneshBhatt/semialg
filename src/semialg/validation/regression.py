from __future__ import annotations

import json
from collections.abc import Iterable, Mapping
from dataclasses import dataclass, field
from pathlib import Path

from .corpus import ValidationCase, read_jsonl_cases, write_jsonl_cases
from .runner import CaseValidationResult


@dataclass(frozen=True)
class RegressionRecord:
    """A replayable validation failure with enough context for triage."""

    case: ValidationCase
    reason: str
    solver_formula: str | None = None
    checker_status: tuple[str, ...] = ()
    witness: Mapping[str, str] | None = None
    metadata: Mapping[str, object] = field(default_factory=dict)

    def to_json_dict(self) -> dict[str, object]:
        return {
            "case": self.case.to_json_dict(),
            "reason": self.reason,
            "solver_formula": self.solver_formula,
            "checker_status": list(self.checker_status),
            "witness": None if self.witness is None else dict(self.witness),
            "metadata": dict(self.metadata),
        }

    @classmethod
    def from_json_dict(cls, data: Mapping[str, object]) -> RegressionRecord:
        return cls(
            case=ValidationCase.from_json_dict(data["case"]),  # type: ignore[arg-type]
            reason=str(data["reason"]),
            solver_formula=data.get("solver_formula"),  # type: ignore[arg-type]
            checker_status=tuple(str(s) for s in data.get("checker_status", ())),  # type: ignore[arg-type]
            witness=data.get("witness"),  # type: ignore[arg-type]
            metadata=dict(data.get("metadata", {})),  # type: ignore[arg-type]
        )


def record_from_result(result: CaseValidationResult) -> RegressionRecord | None:
    if result.passed:
        return None
    witness = None
    for check in result.equivalence_checks:
        if check.witness is not None:
            witness = {str(k): str(v) for k, v in check.witness.assignment.items()}
            break
    status = tuple(f"{o.checker_name}:{o.status}" for o in result.checker_results)
    reason = "; ".join(result.diagnostics) if result.diagnostics else result.solver_status
    return RegressionRecord(
        case=result.case,
        reason=reason,
        solver_formula=result.solver_formula,
        checker_status=status,
        witness=witness,
    )


def write_records(path: str | Path, records: Iterable[RegressionRecord]) -> None:
    out = Path(path)
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", encoding="utf-8") as handle:
        for record in records:
            handle.write(json.dumps(record.to_json_dict(), sort_keys=True) + "\n")


def read_records(path: str | Path) -> tuple[RegressionRecord, ...]:
    records: list[RegressionRecord] = []
    with Path(path).open("r", encoding="utf-8") as handle:
        for line in handle:
            text = line.strip()
            if text:
                records.append(RegressionRecord.from_json_dict(json.loads(text)))
    return tuple(records)


def write_failing_cases(
    path: str | Path, results: Iterable[CaseValidationResult]
) -> tuple[RegressionRecord, ...]:
    records = tuple(
        record for result in results if (record := record_from_result(result)) is not None
    )
    write_records(path, records)
    return records


def export_cases(path: str | Path, records: Iterable[RegressionRecord]) -> None:
    write_jsonl_cases(path, (record.case for record in records))


def import_cases(path: str | Path) -> tuple[ValidationCase, ...]:
    return read_jsonl_cases(path)


__all__ = [
    "RegressionRecord",
    "record_from_result",
    "write_records",
    "read_records",
    "write_failing_cases",
    "export_cases",
    "import_cases",
]
