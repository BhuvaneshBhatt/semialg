"""Exact contracts for CAD root ordering, reuse, pruning and witnesses."""

import pytest
import sympy as sp

from semialg import as_cad_region
from semialg.algebraic.comparison import compare_samples
from semialg.algebraic.intervals import RationalInterval
from semialg.algebraic.samples import AlgebraicRoot, RationalSample
from semialg.cad_algorithms.decomposition import (
    _sample_for_fiber_arithmetic,
    decomp_collins_complete,
)
from semialg.cad_region import cad_combine
from semialg.exact_arithmetic import compare_exact_reals
from semialg.formula import Atom, parse_formula
from semialg.partial.pruning import evaluate_pruning_status
from semialg.partial.qe import lazy_find_inst_form, lazy_resolve_formula
from semialg.qe.complete import qe_by_complete_cad
from semialg.qe.virtual_substitution import try_quadratic_vs_witness
from semialg.qe.witnesses import existential_wit_assign, witness_from_qe_result
from semialg.reconstruct.root_functions import root_of
from semialg.status import SolverStatus


def test_repeated_fiber_root_keeps_its_rank():
    x = sp.Symbol("x", real=True)
    root = AlgebraicRoot(sp.Poly((x**2 - 2) ** 2, x), RationalInterval(1, 2), 1, 2)
    assert compare_exact_reals(_sample_for_fiber_arithmetic(root), sp.sqrt(2)) == 0


@pytest.mark.parametrize("earlier", [-1, sp.sqrt(2)])
def test_repeated_family_endpoint_equality(earlier):
    x = sp.Symbol("x", real=True)
    root = AlgebraicRoot(
        sp.Poly((x - earlier) ** 2 * (x - 2), x, extension=True),
        RationalInterval(sp.Rational(3, 2), 2),
        1,
    )
    assert compare_samples(root, RationalSample(2)) == 0
    assert compare_samples(RationalSample(2), root) == 0


@pytest.mark.parametrize("operator", ["=", "!="])
def test_opaque_root_does_not_prove_pruning(operator):
    x = sp.Symbol("x", real=True)
    root = root_of(x**3 + x - sp.sqrt(2), x, 0)
    decision = evaluate_pruning_status(Atom(root**3 + root - sp.sqrt(2), operator), {})
    assert decision.status is None
    assert not decision.should_prune


@pytest.mark.parametrize("truth", [sp.true, sp.false])
def test_constant_lazy_sentence(truth):
    result = lazy_resolve_formula((), (), parse_formula(truth))
    assert result.truth_value is bool(truth)
    assert result.status == "complete"
    assert result.stats.lifted_stacks == 0


@pytest.mark.parametrize("truth", [sp.true, sp.false])
def test_constant_lazy_instance(truth):
    result = lazy_find_inst_form((), parse_formula(truth))
    assert result.found is bool(truth)
    assert result.instance == ({} if truth is sp.true else None)


def test_short_circuit_instance_assigns_all_coordinates():
    x, y = sp.symbols("x y", real=True)
    formula = sp.Or(x < 0, y < 0)
    result = lazy_find_inst_form((x, y), parse_formula(formula))
    assert set(result.instance) == {x, y}
    assert bool(formula.subs(result.instance))


def test_failed_witness_search_is_unknown():
    x, y = sp.symbols("x y", real=True)
    result = try_quadratic_vs_witness(
        (x, y), (("exists", x),), sp.And(x**2 < 1, y > 5), lambda formula, variables: None
    )
    assert result.status == SolverStatus.UNKNOWN
    assert result.instance is None
    assert bool(result.reduced_formula.subs(y, 6))


def test_ground_false_vs_result_proves_unsat():
    x = sp.Symbol("x", real=True)
    result = try_quadratic_vs_witness(
        (x,), (("exists", x),), x**2 < 0, lambda formula, variables: None
    )
    assert result.status == SolverStatus.UNSAT
    assert result.reduced_formula is sp.false


def test_candidate_witness_preserves_fixed_coordinates():
    x, y = sp.symbols("x y", real=True)
    formula = parse_formula(y > x)
    decomposition = decomp_collins_complete((y - x,), (x, y))
    witness = existential_wit_assign(formula, (x, y), {x: 10}, decomposition.cells)
    assert witness is None


def test_projected_witness_completes_partial_point():
    x, y = sp.symbols("x y", real=True)
    result = qe_by_complete_cad((x, y), (), parse_formula(y > x), return_result=True)
    witness = witness_from_qe_result(result, {x: 10})
    assert witness.assignment[x] == 10
    assert witness.assignment[y] > 10
    assert (
        witness.source_cell.index
        == result.cad._point_locator.locate(
            tuple(witness.assignment[v] for v in result.free_variables)
        ).index
    )


def test_complete_point_has_matching_source_cell():
    x, y = sp.symbols("x y", real=True)
    result = qe_by_complete_cad((x, y), (), parse_formula(y > x), return_result=True)
    witness = witness_from_qe_result(result, {x: 10, y: 11})
    located = result.cad._point_locator.locate(
        tuple(witness.assignment[var] for var in result.free_variables)
    )
    assert witness.assignment == {x: 10, y: 11}
    assert witness.source_cell.index == located.index


def test_finite_cad_complement_has_ambient_cells():
    x = sp.Symbol("x", real=True)
    region = as_cad_region(sp.Eq(x, 0), (x,))
    complement = cad_combine(region, op=sp.Not)
    assert complement.signature().region_dimension == 1
    assert complement.result.as_function()(1)
    assert not complement.result.as_function()(0)


def test_callable_combination_refines_new_boundaries():
    x = sp.Symbol("x", real=True)
    region = as_cad_region(x > 0, (x,))
    combined = cad_combine(region, op=lambda formula: sp.Or(formula, x < -1))
    assert combined.result.as_function()(-2)
    assert not combined.result.as_function()(-1)
    assert not combined.result.as_function()(0)
    assert combined.result.as_function()(1)
    assert not combined.result.diagnostics["cad_reused"]


def test_xor_builds_common_refinement():
    x = sp.Symbol("x", real=True)
    left = as_cad_region(x > 0, (x,))
    right = as_cad_region(x > 1, (x,))
    combined = cad_combine(left, right, op="xor")
    for point, expected in [(0, False), (sp.Rational(1, 2), True), (1, True), (2, False)]:
        assert combined.result.as_function()(point) is expected


def test_endpoint_root_is_not_a_different_section():
    x = sp.Symbol("x", real=True)
    root = AlgebraicRoot(sp.Poly((x - 1) ** 2 * (x - 2), x), RationalInterval(1, 2), 1)
    assert compare_samples(root, RationalSample(1)) == 1
    assert compare_samples(RationalSample(1), root) == -1


def test_incomplete_qe_has_no_certified_witness():
    from dataclasses import replace

    from semialg.errors import ExactComputationUnavailable

    x = sp.Symbol("x", real=True)
    result = qe_by_complete_cad((x,), (), parse_formula(x > 0), return_result=True)
    with pytest.raises(ExactComputationUnavailable, match="incomplete"):
        witness_from_qe_result(replace(result, status="partial"))


def test_empty_point_search_does_not_prove_unsat():
    x, y = sp.symbols("x y", real=True)
    result = try_quadratic_vs_witness(
        (x, y), (("exists", x), ("exists", y)), x**2 + y**2 < 1, lambda formula, variables: None
    )
    assert result.status == SolverStatus.SAT
    assert bool((x**2 + y**2 < 1).subs(result.instance))


def test_point_outside_projected_solution_has_no_witness():
    x = sp.Symbol("x", real=True)
    result = qe_by_complete_cad((x,), (), parse_formula(x > 0), return_result=True)
    assert witness_from_qe_result(result, {x: -1}) is None


def test_finite_cad_boolean_union_reuses_cells():
    x = sp.Symbol("x", real=True)
    region = as_cad_region(sp.Eq(x, 0), (x,))
    combined = cad_combine(region, region, op="or")
    assert combined.result.cad is region.result.cad
    assert combined.signature().region_dimension == 0
    assert combined.result.as_function()(0)
    assert not combined.result.as_function()(1)


@pytest.mark.parametrize("multiplicities", [(1, 2, 3), (3, 1, 2), (2, 3, 1)])
@pytest.mark.parametrize("scale", [-1, 1])
def test_rational_endpoint_rank_family(multiplicities, scale):
    x = sp.Symbol("x", real=True)
    roots = (-2, 0, 3)
    polynomial = sp.Poly(
        scale
        * sp.prod((x - root) ** power for root, power in zip(roots, multiplicities, strict=True)),
        x,
    )
    for index, (root, power) in enumerate(zip(roots, multiplicities, strict=True)):
        sample = AlgebraicRoot(
            polynomial, RationalInterval(root - sp.Rational(1, 4), root), index, power
        )
        assert compare_samples(sample, RationalSample(root)) == 0
        assert compare_samples(sample, RationalSample(root - sp.Rational(1, 8))) == 1
        assert compare_samples(sample, RationalSample(root + sp.Rational(1, 8))) == -1
