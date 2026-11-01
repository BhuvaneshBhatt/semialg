"""Regression contracts for reduced-CAD proofs and quantified queries."""

from dataclasses import replace

import pytest
import sympy as sp

from semialg import cad, find_instance
from semialg.cad_algorithms.decomposition import decomp_collins_complete
from semialg.cad_algorithms.reduced import (
    _certify_sign_invariance,
    certify_truth_inv,
    decomp_form_reduced_safe,
)
from semialg.errors import ExactComputationUnavailable
from semialg.exact_arithmetic import exact_sign
from semialg.formula import parse_formula
from semialg.qe.complete import qe_by_complete_cad, qe_from_cad
from semialg.simplify.cell_union import cell_union_to_formula
from semialg.solve.find_instance import find_instance_text
from semialg.status import SolverStatus
from semialg.tti import extract_formula_families
from semialg.tticad.safe import certify_truth_table_inv


def test_lazard_algebraic_nullification():
    x, y, z = sp.symbols("x y z", real=True)
    polynomial = x**5 - x - 1
    matrix = parse_formula(sp.And(sp.Eq(polynomial * z + y, 0), y < 5))
    reduced = decomp_form_reduced_safe(matrix, (x, y, z), backend="lazard")
    result = qe_from_cad(reduced.cad, (x, y, z), (("exists", z),), matrix)
    alpha = sp.CRootOf(polynomial, 0)
    assert bool(result.subs({x: alpha, y: 0}))
    assert not bool(result.subs({x: alpha, y: 1}))
    assert reduced.used_fallback or reduced.certificate.checked_refinement_cells > 0


@pytest.mark.slow
def test_lazard_quadratic_degree_drop():
    x, a, b, y = sp.symbols("x a b y", real=True)
    matrix = parse_formula(sp.And(sp.Eq(x * y**2 + a * y + b, 0), b < 5))
    reduced = decomp_form_reduced_safe(matrix, (x, a, b, y), backend="lazard")
    result = qe_from_cad(reduced.cad, (x, a, b, y), (("exists", y),), matrix)
    for values, expected in [
        ((0, 0, 0), True),
        ((0, 0, 1), False),
        ((0, 1, 1), True),
        ((1, 0, 1), False),
        ((-1, 0, 1), True),
    ]:
        assert bool(result.subs(dict(zip((x, a, b), values, strict=True)))) is expected


@pytest.mark.parametrize("invariant", ["sign", "truth", "truth-table"])
@pytest.mark.parametrize("reference", ["empty", "missing-boundary"])
def test_refinement_requires_valid_reference(invariant, reference):
    x = sp.Symbol("x", real=True)
    reduced = decomp_collins_complete((1,), (x,))
    refiner = (
        cad(sp.Eq(x**2 + 1, 0, evaluate=False), (x,), return_result=True).cad
        if reference == "empty"
        else reduced
    )
    matrix = parse_formula(x > 0)
    if invariant == "sign":
        result = _certify_sign_invariance(
            reduced=reduced, refiner=refiner, input_polys=(x,), variables=(x,)
        )
    elif invariant == "truth":
        result = certify_truth_inv(reduced=reduced, refiner=refiner, formula=matrix, variables=(x,))
    else:
        result = certify_truth_table_inv(
            reduced=reduced,
            refiner=refiner,
            families=extract_formula_families(matrix),
            variables=(x,),
        )
    assert not result.valid


@pytest.mark.parametrize("mutation", ["incomplete", "variable-order", "projection-proof"])
def test_refinement_rejects_invalid_provenance(mutation):
    x, y = sp.symbols("x y", real=True)
    source = decomp_collins_complete((x, y), (x, y))
    variables = (x, y)
    if mutation == "incomplete":
        refiner = replace(source, complete=False)
    elif mutation == "variable-order":
        refiner = source
        variables = (y, x)
    else:
        refiner = replace(source, tower=replace(source.tower, metadata={"complete": False}))
    certificate = certify_truth_inv(
        reduced=source, refiner=refiner, formula=parse_formula(x > 0), variables=variables
    )
    assert not certificate.valid
    assert certificate.checked_refinement_cells == 0


def test_refinement_accepts_factor_boundaries():
    x, y = sp.symbols("x y", real=True)
    source = decomp_collins_complete((y**2 - x,), (x, y))
    certificate = _certify_sign_invariance(
        reduced=source, refiner=source, input_polys=(-7 * x**3,), variables=(x, y)
    )
    assert certificate.valid
    assert certificate.checked_refinement_cells == len(source.cells)
    assert {values[0] for values in certificate.grouped_values.values()} == {(-1,), (0,), (1,)}


def test_qe_rejects_incomplete_ambient_source():
    x = sp.Symbol("x", real=True)
    source = replace(decomp_collins_complete((x,), (x,)), complete=False)
    with pytest.raises(ExactComputationUnavailable):
        qe_from_cad(source, (x,), (("exists", x),), parse_formula(x > 0))


@pytest.mark.parametrize("quantifier", ["exists", "forall"])
def test_finite_source_qe_uses_ambient_semantics(quantifier):
    x, y = sp.symbols("x y", real=True)
    expression = sp.And(sp.Eq(y**2, 1), sp.Eq(x, y))
    source = cad(expression, (x, y), return_result=True)
    result = qe_from_cad(source.cad, (x, y), ((quantifier, y),), parse_formula(expression))
    for point in (-2, -1, 0, 1, 2):
        expected = quantifier == "exists" and point in (-1, 1)
        assert bool(result.subs(x, point)) is expected


def test_finite_source_alternating_qe():
    x, y = sp.symbols("x y", real=True)
    expression = sp.And(sp.Eq(y**2, 1), sp.Eq(x, y))
    source = cad(expression, (x, y), return_result=True)
    result = qe_from_cad(
        source.cad, (x, y), (("forall", x), ("exists", y)), parse_formula(expression)
    )
    assert result is sp.false


def test_filtered_variety_cannot_discard_new_solutions():
    x, y = sp.symbols("x y", real=True)
    expression = sp.And(sp.Eq(y**2, 1), sp.Eq(x, y))
    source = cad(sp.And(expression, y > 0), (x, y), return_result=True)
    result = qe_from_cad(source.cad, (x, y), (("exists", y),), parse_formula(expression))
    assert bool(result.subs(x, -1))
    assert bool(result.subs(x, 1))
    assert not bool(result.subs(x, 0))


def test_finite_union_is_not_universal():
    x, y = sp.symbols("x y", real=True)
    matrix = parse_formula(sp.And(sp.Eq(y**2, 1), sp.Eq(x, y)))
    result = qe_by_complete_cad(
        (x, y), (("exists", y),), matrix, use_presolve=False, return_result=True
    )
    formula = cell_union_to_formula(result.cell_union)
    for point in (-2, -1, 0, 1, 2):
        assert bool(formula.subs(x, point)) is (point in (-1, 1))


def test_finite_union_keeps_coordinate_relations():
    x, y = sp.symbols("x y", real=True)
    matrix = parse_formula(sp.And(sp.Eq(y**2, 1), sp.Eq(x, y)))
    result = qe_by_complete_cad((x, y), (), matrix, use_presolve=False, return_result=True)
    formula = cell_union_to_formula(result.cell_union)
    for point, expected in [((1, 1), True), ((-1, -1), True), ((1, -1), False), ((0, 0), False)]:
        assert bool(formula.subs(dict(zip((x, y), point, strict=True)))) is expected


@pytest.mark.parametrize("problem", ["contradiction", "parameter"])
def test_quantified_instances_respect_assumptions(problem):
    x, y = sp.symbols("x y", real=True)
    text, assumptions = (
        ("exists x. x < 0", x > 0) if problem == "contradiction" else ("exists y. y**2 = x", x < 0)
    )
    result = find_instance_text(
        text, symbols={"x": x, "y": y}, assumptions=assumptions, return_result=True
    )
    assert result.status == SolverStatus.UNSAT
    assert result.result is None


def test_lazy_instance_respects_parameter_bounds():
    x, y = sp.symbols("x y", real=True)
    result = find_instance_text(
        "exists y. y**3 = x", symbols={"x": x, "y": y}, assumptions=x > 2, return_result=True
    )
    assert result.status == SolverStatus.SAT
    assert result.result[x] > 2
    assert exact_sign(result.result[y] ** 3 - result.result[x]) == 0


def test_assumptions_supply_missing_coordinates():
    x, y = sp.symbols("x y", real=True)
    result = find_instance_text(
        "exists x. x > 0", symbols={"x": x, "y": y}, assumptions=y > 3, return_result=True
    )
    assert set(result.result) == {x, y}
    assert result.result[x] > 0 and result.result[y] > 3


@pytest.mark.parametrize("interface", ["expression", "text"])
@pytest.mark.parametrize("count", [0, -1])
def test_zero_count_is_not_infeasibility(interface, count):
    x = sp.Symbol("x", real=True)
    result = (
        find_instance(x > 0, (x,), count=count, return_result=True)
        if interface == "expression"
        else find_instance_text(
            "exists x. x > 0", symbols={"x": x}, count=count, return_result=True
        )
    )
    assert result.status == SolverStatus.UNKNOWN
    assert result.instances == ()
