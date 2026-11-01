"""Degree and nullification conditions for intrinsic EC certification."""

import pytest
import sympy as sp

from semialg.cad_algorithms.reduced import decomp_form_reduced_safe
from semialg.formula import parse_formula
from semialg.qe.complete import qe_from_cad


@pytest.mark.parametrize("backend", ["mccallum", "lazard"])
@pytest.mark.parametrize("shift", [-2, 0, 1])
@pytest.mark.parametrize("quantifier", ["exists", "forall"])
def test_hidden_nullification(backend, shift, quantifier):
    x, y, z = sp.symbols("x y z", real=True)
    matrix = parse_formula(sp.Eq(y * z + x - shift, 0))
    result = decomp_form_reduced_safe(matrix, (x, y, z), backend=backend)
    projected = qe_from_cad(result.cad, (x, y, z), ((quantifier, z),), matrix)
    for x_value in (shift - 1, shift, shift + 1):
        for y_value in (-2, 0, 2):
            expected = (
                y_value != 0 or x_value == shift
                if quantifier == "exists"
                else y_value == 0 and x_value == shift
            )
            assert bool(projected.subs({x: x_value, y: y_value})) is expected
    assert result.used_fallback or result.certificate.checked_refinement_cells > 0


@pytest.mark.parametrize("backend", ["mccallum", "lazard"])
def test_positive_dimensional_nullification(backend):
    x, y, w, z = sp.symbols("x y w z", real=True)
    variables = (x, y, w, z)
    matrix = parse_formula(sp.Eq(w * z + x - 1, 0))
    result = decomp_form_reduced_safe(matrix, variables, backend=backend)
    projected = qe_from_cad(result.cad, variables, (("exists", z),), matrix)
    for y_value in (-3, 0, 4):
        assert bool(projected.subs({x: 1, y: y_value, w: 0}))
        assert not bool(projected.subs({x: 2, y: y_value, w: 0}))


@pytest.mark.parametrize("backend", ["mccallum", "lazard"])
@pytest.mark.parametrize("leading", [1, "positive"])
def test_degree_invariant_intrinsic_path(backend, leading):
    x, y, z = sp.symbols("x y z", real=True)
    coefficient = 1 if leading == 1 else y**2 + 1
    matrix = parse_formula(sp.Eq(coefficient * z + x - 1, 0))
    result = decomp_form_reduced_safe(matrix, (x, y, z), backend=backend)
    assert not result.used_fallback
    assert result.certificate.checked_refinement_cells == 0
    assert qe_from_cad(result.cad, (x, y, z), (("exists", z),), matrix) is sp.true


@pytest.mark.parametrize("backend", ["mccallum", "lazard"])
def test_recursive_reduction_preserves_order(backend):
    x, y, z = sp.symbols("x y z", real=True)
    matrix = parse_formula(sp.And(sp.Eq(y - x, 0), sp.Eq(z**2 - y, 0)))
    result = decomp_form_reduced_safe(matrix, (x, y, z), backend=backend)
    projected = qe_from_cad(result.cad, (x, y, z), (("exists", y), ("exists", z)), matrix)
    for value in (-2, -1, 0, 1, 2):
        assert bool(projected.subs(x, value)) is (value >= 0)
    assert not result.used_fallback
    assert result.certificate.checked_refinement_cells == 0
    assert any(plan.operator == "ec-order" for plan in result.cad.tower.metadata["ec_levels"])
