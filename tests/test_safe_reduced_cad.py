from __future__ import annotations

import sympy as sp

from semialg.cad.projection.mccallum import build_mccallum_proj_set
from semialg.cad.reduced import decompose_reduced_safe
from semialg.formula import parse_formula
from semialg.tticad.safe import decompose_tticad_safe


def test_safe_reduced_01() -> None:
    x, y = sp.symbols("x y")
    tower = build_mccallum_proj_set(
        [y**2 - x, y - 1],
        [x, y],
        equational_constraints=[y**2 - x],
    )
    assert tower.requested_theory == "mccallum"
    assert tower.validity.complete_if_used is False
    assert tower.validity.fallback_backend == "collins-complete"
    assert tower.tower.variables == (x, y)


def test_safe_reduced_02() -> None:
    x, y = sp.symbols("x y")
    result = decompose_reduced_safe(
        [y**2 - x, y - 1],
        [x, y],
        backend="lazard",
        equational_constraints=[y**2 - x],
    )
    assert result.requested_backend == "lazard"
    assert result.complete is True
    assert result.certificate is not None
    assert result.side_conditions is not None
    assert result.effective_backend in {"collins-complete", "lazard-reduced-certified"}
    assert result.cad.complete is True
    assert result.cad.cells


def test_safe_reduced_03() -> None:
    x, y = sp.symbols("x y")
    formula = parse_formula(sp.Or(sp.Eq(y**2 - x, 0), sp.Eq(y - 1, 0)))
    result = decompose_tticad_safe(formula, [x, y])
    assert result.complete is True
    assert result.family_count == 2
    assert result.certificate is not None
    assert result.side_conditions is not None
    assert result.effective_backend in {"collins-complete", "tticad-reduced-certified"}
    assert result.cad.cells
