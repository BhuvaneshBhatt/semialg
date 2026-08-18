from __future__ import annotations

import json
from collections.abc import Iterable, Mapping
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any

import sympy as sp


@dataclass(frozen=True)
class ValidationCase:
    """A serializable validation input for real-polynomial QE tests.

    The formula is stored as text rather than a Python object so cases can be
    minimized, written to JSONL, and replayed outside the original process.
    """

    name: str
    formula_text: str
    variables: tuple[str, ...]
    quantifiers: tuple[tuple[str, str], ...] = ()
    expected_text: str | None = None
    tags: tuple[str, ...] = ()
    metadata: Mapping[str, Any] = field(default_factory=dict)

    def symbols(self) -> dict[str, sp.Symbol]:
        return {name: sp.Symbol(name, real=True) for name in self.variables}

    def sympy_variables(self) -> tuple[sp.Symbol, ...]:
        syms = self.symbols()
        return tuple(syms[name] for name in self.variables)

    def sympy_quantifiers(self) -> tuple[tuple[str, sp.Symbol], ...]:
        syms = self.symbols()
        return tuple((q, syms[name]) for q, name in self.quantifiers)

    def to_json_dict(self) -> dict[str, Any]:
        data = asdict(self)
        data["variables"] = list(self.variables)
        data["quantifiers"] = [list(pair) for pair in self.quantifiers]
        data["tags"] = list(self.tags)
        data["metadata"] = dict(self.metadata)
        return data

    @classmethod
    def from_json_dict(cls, data: Mapping[str, Any]) -> ValidationCase:
        return cls(
            name=str(data["name"]),
            formula_text=str(data["formula_text"]),
            variables=tuple(str(v) for v in data["variables"]),
            quantifiers=tuple((str(q), str(v)) for q, v in data.get("quantifiers", ())),
            expected_text=data.get("expected_text"),
            tags=tuple(str(t) for t in data.get("tags", ())),
            metadata=dict(data.get("metadata", {})),
        )


def write_jsonl_cases(path: str | Path, cases: Iterable[ValidationCase]) -> None:
    out = Path(path)
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", encoding="utf-8") as handle:
        for case in cases:
            handle.write(json.dumps(case.to_json_dict(), sort_keys=True) + "\n")


def read_jsonl_cases(path: str | Path) -> tuple[ValidationCase, ...]:
    cases: list[ValidationCase] = []
    with Path(path).open("r", encoding="utf-8") as handle:
        for line in handle:
            stripped = line.strip()
            if stripped:
                cases.append(ValidationCase.from_json_dict(json.loads(stripped)))
    return tuple(cases)


def built_in_smoke_cases() -> tuple[ValidationCase, ...]:
    return (
        ValidationCase(
            name="exists_square_root_of_one",
            formula_text="x**2 - 1 = 0",
            variables=("x",),
            quantifiers=(("exists", "x"),),
            expected_text="True",
            tags=("sentence", "univariate", "exists"),
        ),
        ValidationCase(
            name="positive_square_plus_one",
            formula_text="x**2 + 1 > 0",
            variables=("x",),
            quantifiers=(("forall", "x"),),
            expected_text="True",
            tags=("sentence", "univariate", "forall"),
        ),
        ValidationCase(
            name="square_root_projection",
            formula_text="y**2 - x = 0",
            variables=("x", "y"),
            quantifiers=(("exists", "y"),),
            expected_text="x >= 0",
            tags=("projection", "bivariate", "exists"),
        ),
    )


__all__ = ["ValidationCase", "write_jsonl_cases", "read_jsonl_cases", "built_in_smoke_cases"]
