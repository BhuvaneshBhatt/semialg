"""Mathematical contracts for reusing CAD invariance proofs."""

import pytest
import sympy as sp

from semialg import as_cad_region
from semialg.cad_algorithms.decomposition import decomp_collins_complete, decomp_from_proj_tower
from semialg.cad_algorithms.projection.mccallum import build_mccallum_proj_set
from semialg.cad_algorithms.reduced import decomp_form_reduced_safe, decompose_reduced_safe
from semialg.cad_region import cad_extend
from semialg.errors import ExactComputationUnavailable
from semialg.formula import parse_formula
from semialg.qe.complete import qe_from_cad
from semialg.tticad.safe import decompose_tticad_safe


@pytest.mark.parametrize("quantified", [False, True])
@pytest.mark.parametrize("boundary", [sp.Integer(1), sp.sqrt(2)])
def test_qe_refines_missing_boundary(quantified, boundary):
    x = sp.Symbol("x", real=True)
    source = decomp_collins_complete((x,), (x,))
    matrix = parse_formula(sp.And(x > 0, x < boundary))
    prefix = (("exists", x),) if quantified else ()
    result = qe_from_cad(source, (x,), prefix, matrix, backend="reused", return_result=True)
    if quantified:
        assert result.formula is sp.true
    else:
        for point, expected in [
            (0, False),
            (boundary / 2, True),
            (boundary, False),
            (2 * boundary, False),
        ]:
            assert bool(result.formula.subs(x, point)) is expected
    assert result.cad is not source
    assert result.backend != "reused"


@pytest.mark.parametrize("backend", ["mccallum", "lazard", "tticad"])
def test_truth_scope_cannot_prove_signs(backend):
    x, y = sp.symbols("x y", real=True)
    polynomial = y**2 + y + x
    matrix = parse_formula(sp.And(sp.Eq(y, 0), polynomial >= 0))
    safe = (
        decompose_tticad_safe(matrix, (x, y))
        if backend == "tticad"
        else decomp_form_reduced_safe(matrix, (x, y), backend=backend)
    )
    assert not safe.used_fallback
    result = qe_from_cad(
        safe.cad, (x, y), (("exists", y),), parse_formula(polynomial < 0), return_result=True
    )
    for point in (-1, 0, sp.Rational(1, 8), sp.Rational(1, 4), 1):
        assert bool(result.formula.subs(x, point)) is bool(point < sp.Rational(1, 4))
    assert result.cad is not safe.cad


@pytest.mark.parametrize("backend", ["mccallum", "lazard", "tticad"])
@pytest.mark.parametrize("complement", [False, True])
def test_truth_scope_reuses_its_formula(backend, complement):
    x, y = sp.symbols("x y", real=True)
    expression = sp.And(sp.Eq(y, 0), y**2 + y + x >= 0)
    matrix = parse_formula(expression)
    safe = (
        decompose_tticad_safe(matrix, (x, y))
        if backend == "tticad"
        else decomp_form_reduced_safe(matrix, (x, y), backend=backend)
    )
    query = parse_formula(sp.Not(expression)) if complement else matrix
    result = qe_from_cad(safe.cad, (x, y), (("exists", y),), query, return_result=True)
    assert result.cad is safe.cad
    for point in (-1, 0, 1):
        assert bool(result.formula.subs(x, point)) is (complement or point >= 0)


def test_raw_reduced_lifting_has_no_certificate():
    x, y = sp.symbols("x y", real=True)
    tower = build_mccallum_proj_set(
        (y, y**2 + y + x), (x, y), equational_constraints=(y,), certify=True
    )
    raw = decomp_from_proj_tower(tower.tower)
    assert not raw.complete
    with pytest.raises(ExactComputationUnavailable):
        qe_from_cad(raw, (x, y), (("exists", y),), parse_formula(y**2 + y + x < 0))


@pytest.mark.parametrize(
    "expression", [lambda x: -7 * x**3 < 0, lambda x: sp.Ne(x**2, 0), lambda x: x >= 0]
)
def test_sign_scope_reuses_factors(expression):
    x = sp.Symbol("x", real=True)
    safe = decompose_reduced_safe((x,), (x,), backend="mccallum")
    result = qe_from_cad(safe.cad, (x,), (), parse_formula(expression(x)), return_result=True)
    assert result.cad is safe.cad
    for point in (-1, 0, 1):
        assert bool(result.formula.subs(x, point)) is bool(expression(sp.Integer(point)))


def test_region_extension_reuses_factored_boundaries():
    x = sp.Symbol("x", real=True)
    region = as_cad_region(sp.And(x >= 0, x <= 1), (x,))
    extended = cad_extend(region, -7 * x**3 < 0)
    assert extended.cad is region.cad
    for point, expected in [
        (-1, False),
        (0, False),
        (sp.Rational(1, 2), True),
        (1, True),
        (2, False),
    ]:
        assert extended.result.as_function()(point) is expected


def test_finite_extension_keeps_point_semantics():
    x = sp.Symbol("x", real=True)
    region = as_cad_region(sp.Eq(x * x, 1), (x,))
    extended = cad_extend(region, x > 0)
    assert extended.cad is region.cad
    for point in (-1, 0, 1, 2):
        assert extended.result.as_function()(point) is (point == 1)


@pytest.mark.parametrize("truth", [sp.true, sp.false])
def test_zero_coordinate_cad(truth):
    from semialg import cad

    result = cad(truth, (), return_result=True)
    function = result.as_function()
    assert function() is bool(truth)
    assert function.locate_cell(()).index == ()
    location = as_cad_region(result).locate_point(())
    assert location.selected is bool(truth)
    assert location.dimension == 0
    assert location.sign_vector == ()
    check = result.cad.verify_sign_invariance()
    assert check.ok
    assert check.checked_cells == 1
    assert check.checked_polynomials == 0


def test_zero_coordinate_relifting():
    source = decomp_collins_complete((1,), ())
    result = decomp_from_proj_tower(source.tower)
    assert tuple(cell.index for cell in result.cells) == ((),)
    assert result._point_locator.locate(()).sample == ()
    assert result.complete
    assert result.verify_sign_invariance().ok


@pytest.mark.parametrize("coordinate", [-1, 1])
def test_fully_restricted_cad(coordinate):
    x = sp.Symbol("x", real=True)
    original = as_cad_region(x > 0, (x,)).result.as_function()
    result = original.restrict({x: coordinate})
    assert result.variables == result.cad.tower.variables == ()
    assert result() is (coordinate > 0)
    assert result.locate_cell(()).index == ()
    assert result.cad.verify_sign_invariance().ok


def test_universal_query_refines_missing_boundary():
    x = sp.Symbol("x", real=True)
    source = decomp_collins_complete((1,), (x,))
    result = qe_from_cad(source, (x,), (("forall", x),), parse_formula(x < 1))
    assert result is sp.false


def test_full_projection_keeps_sign_scope():
    x = sp.Symbol("x", real=True)
    source = decomp_form_reduced_safe(parse_formula(sp.Or(x < 0, x > 1)), (x,), backend="mccallum")
    result = qe_from_cad(source.cad, (x,), (), parse_formula(sp.Eq(x, 0)), return_result=True)
    assert result.cad is source.cad
    for point in (-1, 0, 1):
        assert bool(result.formula.subs(x, point)) is (point == 0)


@pytest.mark.parametrize("truth", [sp.true, sp.false])
def test_ground_qe_preserves_zero_coordinates(truth):
    from semialg.qe.complete import qe_by_complete_cad

    result = qe_by_complete_cad((), (), parse_formula(truth), return_result=True)
    assert result.variables == result.cad.tower.variables == ()
    assert result.formula is truth
    assert result.cell_truth == {(): bool(truth)}
    assert result.satisfying_cell_indices == (((),) if truth is sp.true else ())
    assert result.witness_samples == (({},) if truth is sp.true else ())
    assert result.diagnostics.full_cell_count == 1
    assert result.diagnostics.projected_cell_count == int(bool(truth))
    reused = qe_from_cad(result.cad, (), (), parse_formula(truth), return_result=True)
    assert reused.cad is result.cad
    assert reused.cell_truth == result.cell_truth
    assert reused.witness_samples == result.witness_samples


@pytest.mark.parametrize("truth", [sp.true, sp.false])
@pytest.mark.parametrize("quantified", [False, True])
def test_constant_qe_records_selected_cells(truth, quantified):
    from semialg.qe.complete import qe_by_complete_cad

    x, y = sp.symbols("x y", real=True)
    prefix = (("exists", y),) if quantified else ()
    result = qe_by_complete_cad(
        (x, y), prefix, parse_formula(truth), use_presolve=False, return_result=True
    )
    assert result.truth_value is None
    assert not result.is_sentence
    assert result.cell_union.cell_indices == result.satisfying_cell_indices
    assert result.diagnostics.projected_cell_count == (1 if truth is sp.true else 0)
    assert result.cell_truth == {(0,) if quantified else (0, 0): bool(truth)}
    assert result.witness_samples == (({x: 0, y: 0},) if truth is sp.true else ())
