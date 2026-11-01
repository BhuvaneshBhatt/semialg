"""Mathematical regressions for CAD fibers, root ranks and QE prefixes."""

import pytest
import sympy as sp

from semialg import cad, quantifier_eliminate
from semialg.cad_algorithms.cells import extract_cylindrical_solution
from semialg.cad_algorithms.structured_cells import extract_structured_cad_cells
from semialg.errors import ResourceLimitError
from semialg.formula import parse_formula
from semialg.qe.complete import evaluate_formula_on_cell
from semialg.qe.virtual_substitution import try_quadratic_vs_qe


@pytest.mark.parametrize(
    "formula, points",
    [
        (lambda x, y: y > x, [(-10, -9), (10, 9), (10, 10), (-10, -10)]),
        (lambda x, y: y <= x**2, [(3, 8), (3, 9), (-3, 10), (sp.sqrt(2), 2)]),
        (
            lambda x, y: x**2 + y**2 <= 1,
            [(sp.Rational(3, 4), sp.Rational(3, 4)), (sp.Rational(3, 5), sp.Rational(4, 5))],
        ),
        (lambda x, y: x * y > 1, [(2, 1), (2, sp.Rational(1, 2)), (-2, -1), (-2, 1)]),
    ],
    ids=["affine", "parabola", "disk", "hyperbola"],
)
def test_point_location_uses_query_fiber(formula, points):
    x, y = sp.symbols("x y", real=True)
    condition = formula(x, y)
    function = cad(condition, (x, y), output="function", strategy="full")
    for point in points:
        expected = bool(condition.subs(dict(zip((x, y), point, strict=True))))
        assert function(*point) is expected
        leaf = function.locate_cell(point)
        assert evaluate_formula_on_cell(parse_formula(condition), leaf, (x, y)) is expected


@pytest.mark.parametrize("full", [False, True])
@pytest.mark.parametrize("kinds", [("exists", "forall"), ("forall", "exists")])
def test_virtual_substitution_preserves_alternation(full, kinds):
    x, y = sp.symbols("x y", real=True)
    prefix = tuple(zip(kinds, (x, y), strict=True))
    assert try_quadratic_vs_qe((x, y), prefix, sp.Eq(x, y), full=full) is None


def test_repeated_roots_have_distinct_ranks():
    x, y = sp.symbols("x y", real=True)
    result = cad(
        sp.And(sp.Eq(x, 0), sp.Eq(y**3 - 3 * y - 2 + x, 0)),
        (x, y),
        strategy="full",
        return_result=True,
    )
    cells = extract_structured_cad_cells(result).cells
    assert {cell.sample[y] for cell in cells} == {-1, 2}
    for cell in cells:
        bound = cell.levels[-1].lower_bound
        assert bound.evaluate({x: 0}).as_expr() == cell.sample[y]
        assert bound.as_expr().subs(x, 0) == cell.sample[y]


@pytest.mark.parametrize("source", ["formula", "result", "cell_set"])
def test_full_extraction_includes_false_cells(source):
    x, y = sp.symbols("x y", real=True)
    condition = sp.And(sp.Eq(x, 1), sp.Eq(y, 2))
    result = cad(condition, (x, y), return_result=True)
    input_ = {"formula": condition, "result": result, "cell_set": result.cell_set}[source]
    solution = extract_cylindrical_solution(input_, (x, y), selected_only=False)
    assert any(cell.selected for cell in solution.cells)
    assert any(not cell.selected for cell in solution.cells)
    for point in [(0, 0), (1, 2), (3, -1)]:
        subs = dict(zip((x, y), point, strict=True))
        containing = [cell for cell in solution.cells if bool(cell.as_formula().subs(subs))]
        assert len(containing) == 1
        assert containing[0].selected is bool(condition.subs(subs))


def test_full_extraction_bypasses_selected_shortcut():
    x, y = sp.symbols("x y", real=True)
    condition = sp.And(x > 1, y**3 + y - x < 0)
    solution = extract_cylindrical_solution(condition, (x, y), selected_only=False)
    assert any(not cell.selected for cell in solution.cells)
    for point in [(0, 0), (2, 5), (2, -5)]:
        subs = dict(zip((x, y), point, strict=True))
        containing = [cell for cell in solution.cells if bool(cell.as_formula().subs(subs))]
        assert len(containing) == 1
        assert containing[0].selected is bool(condition.subs(subs))


def test_unknown_cad_cannot_certify_empty_solution():
    x = sp.Symbol("x", real=True)
    result = cad(x > 0, (x,), max_cells=1, return_result=True)
    assert result.status == "unknown"
    with pytest.raises(ResourceLimitError):
        extract_cylindrical_solution(result)


def test_root_selector_preserves_binding():
    from semialg.reconstruct.root_functions import root_of

    t, a = sp.symbols("t a", real=True)
    root = root_of(t**3 + a * t + 1, t, 0)
    assert root.free_symbols == {a}
    assert root.subs(t, 0) == root
    substituted = root.subs(a, t)
    assert substituted.free_symbols == {t}
    assert substituted.subs(t, 0) == -1
    assert substituted.subs(t, -1) == sp.real_roots(t**3 - t + 1)[0]


def test_refinement_checks_factorized_signs():
    from semialg.cad_algorithms.decomposition import decomp_collins_complete
    from semialg.cad_algorithms.reduced import _certify_sign_invariance

    x = sp.Symbol("x", real=True)
    polynomial = x * (x - 1)
    refiner = decomp_collins_complete((polynomial,), (x,))
    coarsened = decomp_collins_complete((1,), (x,))
    certificate = _certify_sign_invariance(
        reduced=coarsened, refiner=refiner, input_polys=(polynomial,), variables=(x,)
    )
    assert not certificate.valid
    assert set(certificate.grouped_values[(0,)]) == {(-1,), (0,), (1,)}


def test_refinement_locates_moving_sections():
    from semialg.cad_algorithms.decomposition import decomp_collins_complete
    from semialg.cad_algorithms.reduced import certify_truth_inv

    x, y = sp.symbols("x y", real=True)
    reduced = decomp_collins_complete((y - x,), (x, y))
    refiner = decomp_collins_complete((y - x, x - 10), (x, y))
    certificate = certify_truth_inv(
        reduced=reduced, refiner=refiner, formula=parse_formula(y > x), variables=(x, y)
    )
    assert certificate.valid
    assert certificate.checked_refinement_cells == len(refiner.cells)
    assert all(len(values) == 1 for values in certificate.grouped_values.values())


def test_decomposition_certificate_checks_coverage():
    from dataclasses import replace

    from semialg.cad_algorithms.cells import cylindrical_solution_from_structured

    x = sp.Symbol("x", real=True)
    decomposition = extract_structured_cad_cells(x**2 <= 1, (x,))
    complete = cylindrical_solution_from_structured(decomposition)
    assert complete.decomposition_cert.certified
    omitted = cylindrical_solution_from_structured(
        replace(decomposition, cells=decomposition.cells[:-1])
    )
    assert not omitted.decomposition_cert.coverage_verified
    duplicated = cylindrical_solution_from_structured(
        replace(decomposition, cells=decomposition.cells + decomposition.cells[:1])
    )
    assert not duplicated.decomposition_cert.pairwise_disjoint
    source_free = cylindrical_solution_from_structured(replace(decomposition, cad_result=None))
    assert not source_free.decomposition_cert.coverage_verified


@pytest.mark.parametrize("strategy", ["auto", "full"])
def test_root_ranks_over_algebraic_base(strategy):
    x, y = sp.symbols("x y", real=True)
    result = cad(
        sp.And(sp.Eq(x**2, 2), sp.Eq(y**2, x)), (x, y), strategy=strategy, return_result=True
    )
    solution = extract_cylindrical_solution(result)
    assert {cell.levels[-1].typed_lower.root_index for cell in solution.cells} == {0, 1}
    for cell in solution.cells:
        bound = cell.levels[-1].typed_lower
        expected = cell.sample[y]
        assert not expected.free_symbols
        actual = bound.evaluate({x: sp.sqrt(2)}).as_expr()
        assert sp.simplify(actual - expected) == 0
        assert sp.simplify(bound.as_expr().subs(x, sp.sqrt(2)) - expected) == 0


@pytest.mark.parametrize("strategy", ["auto", "collins"])
@pytest.mark.parametrize("prefix, expected", [("exists", sp.false), ("forall", sp.true)])
def test_alternating_prefix_truth(strategy, prefix, expected):
    x, y = sp.symbols("x y", real=True)
    inner = "forall" if prefix == "exists" else "exists"
    assert (
        quantifier_eliminate(sp.Eq(x, y), ((prefix, x), (inner, y)), strategy=strategy) is expected
    )


def test_nested_fiber_point_location():
    x, y, z = sp.symbols("x y z", real=True)
    formula = sp.And(y >= x, z > x + y)
    function = cad(formula, (x, y, z), output="function", strategy="full")
    for point in [(10, 11, 22), (10, 11, 20), (-10, -9, -18), (10, 9, 30)]:
        subs = dict(zip((x, y, z), point, strict=True))
        assert function(*point) is bool(formula.subs(subs))


def test_exterior_retains_internal_sections():
    x, y = sp.symbols("x y", real=True)
    condition = sp.And(y > x, y < -x)
    exterior = cad(condition, (x, y), operation="exterior", return_result=True)
    expected = sp.Or(y < x, y > -x)
    for xv in range(-2, 3):
        for yv in range(-2, 3):
            point = {x: xv, y: yv}
            assert bool(exterior.formula.subs(point)) is bool(expected.subs(point))
    assert exterior.as_function()(0, 1)
    assert not exterior.as_function()(0, 0)


def test_closure_requires_boundary_refinement():
    from semialg.errors import ExactComputationUnavailable

    x, y, z = sp.symbols("x y z", real=True)
    condition = sp.And(x > 0, sp.Eq(y, x * z))
    # The closure contains the full z-axis at x=y=0. Substituting those
    # coordinates into the path z=y/x cannot establish its incidence.
    with pytest.raises(ExactComputationUnavailable, match="refinement"):
        cad(condition, (x, y, z), operation="closure")


@pytest.mark.parametrize("method", ["project", "project_cell_set"])
def test_empty_projection_is_not_false(method):
    x = sp.Symbol("x", real=True)
    function = cad(x > 0, (x,), output="function")
    with pytest.raises(NotImplementedError, match="zero-coordinate"):
        getattr(function, method)(())
