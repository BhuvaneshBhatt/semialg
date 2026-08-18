import sympy as sp

from semialg.cad import (
    ProjectionPolynomial,
    build_collins_proj_set,
    decomp_collins_complete,
    verify_cad_sign_inv,
)


def test_auditable_coll_01():
    x, y = sp.symbols("x y")
    tower = build_collins_proj_set([y**2 - x, y - 1], [x, y])
    level_one = tower.level(1)
    assert level_one.entries
    assert all(isinstance(entry, ProjectionPolynomial) for entry in level_one.entries)
    sources = {entry.source for entry in level_one.entries}
    assert any("discriminant" in source or "resultant" in source for source in sources)
    assert tower.metadata["provenance"] is True
    assert tower.poly_count_by_level()[2] == 2


def test_auditable_coll_02():
    x = sp.symbols("x")
    cad = decomp_collins_complete([x**2 - 1], [x])
    assert [cell.kind for cell in cad.cells] == ["sector", "section", "sector", "section", "sector"]
    assert [cell.stack_position for cell in cad.cells] == [0, 1, 2, 3, 4]
    assert [cell.index for cell in cad.cells] == [(0,), (1,), (2,), (3,), (4,)]
    assert cad.cells[1].root_index == 0
    assert cad.cells[3].root_index == 1
    assert cad.cells[1].defining_polynomial_key in {
        sp.sstr(x**2 - 1),
        sp.sstr(x + 1),
        sp.sstr(x - 1),
    }


def test_auditable_coll_03():
    x = sp.symbols("x")
    cad = decomp_collins_complete([x**2 - 1], [x])
    assert cad.cell_count_by_level() == {1: 5}
    assert cad.proj_poly_count_by_level()[1] >= 1
    assert cad.max_stack_size() == 5
    assert cad.diagnostics is not None
    assert cad.diagnostics.cell_count_by_level == {1: 5}
    assert cad.diagnostics.proj_poly_count_by_level[1] >= 1
    assert cad.diagnostics.max_stack_size == 5
    assert "projection" in cad.diagnostics.timing_by_stage
    assert "lifting" in cad.diagnostics.timing_by_stage
    assert cad.diagnostics.invariant_failures == ()


def test_auditable_coll_04():
    x = sp.symbols("x")
    cad = decomp_collins_complete([x**2 - 1], [x])
    check = verify_cad_sign_inv(cad.cells_by_level, cad.tower)
    assert check.ok, check.failures
    assert check.checked_cells == 5
    assert check.checked_polynomials >= 5
