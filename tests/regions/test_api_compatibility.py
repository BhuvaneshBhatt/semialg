from __future__ import annotations

import sympy as sp

from semialg import (
    BoxRegion,
    CADRegion,
    SemialgebraicRegion,
    as_cad_region,
    as_semialgebraic_region,
    region_intersection,
)
from semialg.reasoning import region_subset


def test_four_region_entry_styles_interoperate():
    x, y = sp.symbols("x y", real=True)
    formula = sp.And(x >= 0, x <= 1, y >= 0, y <= 1)
    standard = BoxRegion(((0, 1), (0, 1)))
    symbolic = SemialgebraicRegion(formula, (x, y))
    cad = as_cad_region(symbolic)

    assert isinstance(as_semialgebraic_region(standard, (x, y)), SemialgebraicRegion)
    assert isinstance(cad, CADRegion)
    assert region_subset(standard, symbolic, (x, y))

    mixed = region_intersection(symbolic, standard)
    assert isinstance(mixed, SemialgebraicRegion)
    assert mixed.equals_region(symbolic)
