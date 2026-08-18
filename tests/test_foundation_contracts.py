from __future__ import annotations

import pytest
import sympy as sp

from semialg.algebraic import (
    AlgebraicRoot,
    RationalSample,
    compare_samples,
    isolate_real_roots,
    sign_at_sample,
)
from semialg.cad import build_collins_proj_set, decomp_collins_complete
from semialg.cad.reduced import decompose_reduced_safe
from semialg.formula import parse_quant_form_text
from semialg.qe.complete import qe_by_complete_cad
from semialg.validation import built_in_smoke_cases, run_validation_cases


def test_foundation_con_01() -> None:
    x, y = sp.symbols("x y")
    tower = build_collins_proj_set([(y**2 - x) * (y - 1), x + 2], [x, y])
    top = {sp.factor(poly.as_expr()) for poly in tower.level(2).polynomials}
    base = {sp.factor(poly.as_expr()) for poly in tower.level(1).polynomials}
    assert sp.sstr(x - y**2) in {sp.sstr(p) for p in top}
    assert y - 1 in top
    assert sp.sstr((x - y**2) * (y - 1)) in {sp.sstr(p) for p in top}
    assert x + 2 in base
    sources = {entry.source for entry in tower.level(1).entries}
    assert {"coefficient", "discriminant", "resultant"} & sources


def test_foundation_con_02() -> None:
    x = sp.symbols("x")
    cad = decomp_collins_complete([x**2 - 1], [x])
    assert cad.complete
    assert cad.verify_sign_invariance().ok
    assert all(isinstance(cell.sample[0], (RationalSample, AlgebraicRoot)) for cell in cad.cells)


def test_foundation_con_03() -> None:
    x = sp.symbols("x")
    root = isolate_real_roots((x - 2) ** 3, x)[0]
    assert root.multiplicity == 3
    with pytest.raises(TypeError):
        sign_at_sample(x - 2, [sp.Integer(2)])
    assert compare_samples(root, RationalSample(2)) == 0


def test_foundation_con_04() -> None:
    x, y = sp.symbols("x y", real=True)
    parsed = parse_quant_form_text("exists y. y^2 - x = 0", symbols={"x": x, "y": y})
    result = qe_by_complete_cad(parsed.vars, parsed.quantifiers, parsed.matrix, free_variables=[x])
    assert result.status == "complete"
    assert result.cell_union is not None
    assert bool(result.formula.subs(x, 4)) is True
    assert bool(result.formula.subs(x, -1)) is False


def test_foundation_con_05() -> None:
    report = run_validation_cases(built_in_smoke_cases())
    assert report.passed, report.to_json()
    x, y = sp.symbols("x y")
    reduced = decompose_reduced_safe(
        [y**2 - x, y - 1], [x, y], backend="mccallum", equational_constraints=[y**2 - x]
    )
    assert reduced.complete
    assert reduced.cad.complete
    assert reduced.certificate is not None
