from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field
from typing import Literal

import sympy as sp

from ...algebraic.samples import Sample, sample_to_expr
from ...algebraic.signs import sign_at_sample
from ..polynomial_utils import polynomial_key as _poly_key

CellKind = Literal["sector", "section"]


@dataclass(frozen=True)
class CADCell:
    level: int
    index: tuple[int, ...]
    sample: tuple[Sample, ...]
    interval: tuple[Sample | None, Sample | None] | None = None
    section_polynomial: sp.Expr | None = None
    signs: Mapping[str, int] = field(default_factory=dict)
    kind: CellKind = "sector"
    parent_index: tuple[int, ...] | None = None
    stack_position: int = 0
    root_index: int | None = None
    defining_polynomial_key: str | None = None

    def __post_init__(self) -> None:
        if self.kind == "section" and self.section_polynomial is None:
            raise ValueError("section cells require a section polynomial")
        if self.kind == "sector" and self.root_index is not None:
            raise ValueError("sector cells must not have a root index")

    @property
    def is_section(self) -> bool:
        return self.kind == "section"

    @property
    def is_sector(self) -> bool:
        return self.kind == "sector"

    @property
    def sample_exprs(self) -> tuple[sp.Expr, ...]:
        return tuple(sample_to_expr(sample) for sample in self.sample)

    @property
    def lower_bound(self) -> Sample | None:
        return None if self.interval is None else self.interval[0]

    @property
    def upper_bound(self) -> Sample | None:
        return None if self.interval is None else self.interval[1]


def sign_table(polys: Sequence[sp.Poly], sample: Sequence[Sample]) -> dict[str, int]:
    return {_poly_key(poly): sign_at_sample(poly, sample) for poly in polys}
