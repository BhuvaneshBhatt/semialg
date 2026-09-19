"""Small exact case matrices shared by related public operations."""

from __future__ import annotations

import pytest
import sympy as sp

from semialg import (
    argmax_set,
    argmin_set,
    classify_real_roots,
    connected_component_count,
    connected_component_samples,
    function_range,
    integrate_over_region,
    is_bijective,
    is_empty,
    is_function_continuous,
    is_function_smooth,
    is_injective,
    is_satisfiable,
    is_surjective,
    is_tautology,
    matrix_pd_on,
    matrix_psd_on,
    matrix_rank_on,
    region_boundary,
    region_closure,
    region_dimension,
    region_interior,
    semialgebraic_maximize,
    semialgebraic_minimize,
)
from semialg.roadmaps import roadmap
from semialg.topology.semialgebraic import component_decomposition

X = sp.Symbol("x", real=True)
Y = sp.Symbol("y", real=True)


@pytest.mark.parametrize("structured", [False, True], ids=["boolean", "result"])
@pytest.mark.parametrize(
    "formula, expected",
    [
        pytest.param(sp.true, True, id="universe"),
        pytest.param(sp.false, False, id="empty"),
        pytest.param(X > 0, True, id="strict-ray"),
        pytest.param((X > 0) & (X <= 0), False, id="strict-conflict"),
        pytest.param(sp.Eq(X**2, 2), True, id="irrational-points"),
        pytest.param(sp.Eq(X**2, -1), False, id="no-real-roots"),
        pytest.param(sp.Ne(X, 0), True, id="punctured-line"),
    ],
)
def test_decision_matrix(formula, expected, structured):
    result = is_satisfiable(formula, (X,), return_result=structured)
    if structured:
        assert result.satisfiable is expected
        assert result.variables == (X,)
        if result.witness is not None:
            assert expected
            assert sp.simplify(formula.subs(result.witness)) is sp.true
    else:
        assert result is expected
    assert is_empty(formula, (X,)) is (not expected)
    assert is_tautology(sp.Not(formula), (X,)) is (not expected)


@pytest.mark.parametrize(
    "source, closed, interior, boundary, dimension, count",
    [
        pytest.param(sp.S.EmptySet, sp.S.EmptySet, sp.S.EmptySet, sp.S.EmptySet, -1, 0, id="empty"),
        pytest.param(sp.S.Reals, sp.S.Reals, sp.S.Reals, sp.S.EmptySet, 1, 1, id="universe"),
        pytest.param(
            sp.FiniteSet(0), sp.FiniteSet(0), sp.S.EmptySet, sp.FiniteSet(0), 0, 1, id="singleton"
        ),
        pytest.param(
            sp.Interval.open(0, 1),
            sp.Interval(0, 1),
            sp.Interval.open(0, 1),
            sp.FiniteSet(0, 1),
            1,
            1,
            id="open",
        ),
        pytest.param(
            sp.Interval(0, 1),
            sp.Interval(0, 1),
            sp.Interval.open(0, 1),
            sp.FiniteSet(0, 1),
            1,
            1,
            id="closed",
        ),
        pytest.param(
            sp.Interval.Ropen(0, 1),
            sp.Interval(0, 1),
            sp.Interval.open(0, 1),
            sp.FiniteSet(0, 1),
            1,
            1,
            id="half-open",
        ),
        pytest.param(
            sp.S.Reals - sp.FiniteSet(0),
            sp.S.Reals,
            sp.S.Reals - sp.FiniteSet(0),
            sp.FiniteSet(0),
            1,
            2,
            id="punctured-line",
        ),
    ],
)
def test_topology_matrix(source, closed, interior, boundary, dimension, count):
    formula = source.as_relational(X)
    assert region_closure(formula, (X,)).as_set() & sp.S.Reals == closed
    assert region_interior(formula, (X,)).as_set() & sp.S.Reals == interior
    assert region_boundary(formula, (X,)).as_set() & sp.S.Reals == boundary
    assert region_dimension(formula, (X,)) == dimension
    assert connected_component_count(formula, (X,)) == count


@pytest.mark.parametrize("left_open", [False, True], ids=["left-closed", "left-open"])
@pytest.mark.parametrize("right_open", [False, True], ids=["right-closed", "right-open"])
@pytest.mark.parametrize("slope", [-2, 0, 3], ids=["decreasing", "constant", "increasing"])
def test_optimization_endpoint_matrix(left_open, right_open, slope):
    domain = sp.Interval(-1, 2, left_open, right_open)
    formula = domain.as_relational(X)
    objective = slope * X + 1
    low, high = sorted((1 - slope, 1 + 2 * slope))
    min_closed = not (left_open if slope > 0 else right_open) if slope else True
    max_closed = not (right_open if slope > 0 else left_open) if slope else True
    minimum = semialgebraic_minimize(objective, formula, (X,), return_result=True)
    maximum = semialgebraic_maximize(objective, formula, (X,), return_result=True)
    image = function_range(objective, formula, (X,), return_result=True)
    for result, value, attained in ((minimum, low, min_closed), (maximum, high, max_closed)):
        assert result.value == value
        assert result.attained is attained
        assert result.certified
        for point in result.points:
            assert sp.simplify(formula.subs(point)) is sp.true
            assert objective.subs(point) == value
        if not attained:
            assert result.points == ()
    assert (image.infimum, image.supremum) == (low, high)
    assert image.minimum_attained is min_closed
    assert image.maximum_attained is max_closed
    min_set = domain if not slope else sp.FiniteSet(-1 if slope > 0 else 2) & domain
    max_set = domain if not slope else sp.FiniteSet(2 if slope > 0 else -1) & domain
    assert argmin_set(objective, formula, (X,)).as_set() == min_set
    assert argmax_set(objective, formula, (X,)).as_set() == max_set


@pytest.mark.parametrize(
    "polynomial, count, pattern",
    [
        pytest.param(sp.Integer(0), sp.oo, (), id="zero"),
        pytest.param(sp.Integer(5), 0, (), id="nonzero-constant"),
        pytest.param(X - 2, 1, (1,), id="linear"),
        pytest.param((X - 2) ** 3, 1, (3,), id="triple-root"),
        pytest.param((X + 1) ** 2 * (X - 2), 2, (1, 2), id="mixed-multiplicity"),
        pytest.param(X**2 - 2, 2, (1, 1), id="irrational"),
        pytest.param(X**2 + 1, 0, (), id="complex-only"),
    ],
)
def test_root_classification_matrix(polynomial, count, pattern):
    result = classify_real_roots(polynomial, X)
    assert result.parameters == ()
    assert result.generic_root_count == count
    assert result.generic_multiplicity_pattern == pattern
    assert len(result.cells) == 1
    assert result.cells[0].condition is sp.true
    assert result.cells[0].root_count == count


@pytest.mark.parametrize("leading", [-1, 0, 1])
@pytest.mark.parametrize("constant", [-1, 0, 1])
def test_root_degree_drop_matrix(leading, constant):
    a, b = sp.symbols("a b", real=True)
    family = classify_real_roots(a * X + b, X, parameters=(a, b))
    assignment = {a: leading, b: constant}
    selected = [cell for cell in family.cells if cell.condition.subs(assignment) is sp.true]
    expected = 1 if leading else (0 if constant else sp.oo)
    assert len(selected) == 1
    assert selected[0].root_count == expected
    specialized = classify_real_roots(leading * X + constant, X)
    assert specialized.generic_root_count == expected


@pytest.mark.parametrize(
    "formula, variables, count",
    [
        pytest.param(sp.Ne(X, 0), (X,), 2, id="punctured-line"),
        pytest.param(sp.Ne(X * (X - 1), 0), (X,), 3, id="two-punctures"),
        pytest.param(sp.Ne(X, 0), (X, Y), 2, id="separated-half-planes"),
        pytest.param(sp.Ne(X * Y, 0), (X, Y), 4, id="open-quadrants"),
        pytest.param(sp.Eq(X * Y, 0), (X, Y), 1, id="axes-with-bridge"),
    ],
)
def test_component_boundary_matrix(formula, variables, count):
    assert connected_component_count(formula, variables) == count
    decomposition = component_decomposition(formula, variables)
    assert decomposition.count == count
    samples = connected_component_samples(formula, variables)
    assert len(samples) == count
    assert len({tuple(sample[var] for var in variables) for sample in samples}) == count
    assert all(sp.simplify(formula.subs(sample)) is sp.true for sample in samples)
    if len(variables) == 1:
        assert roadmap(formula, variables).component_count == count
        expected = formula.as_set()
        pieces = [component.as_formula().as_set() for component in decomposition.components]
        assert sp.Union(*pieces) == expected
        assert all(piece.is_Interval for piece in pieces)


@pytest.mark.parametrize("structured", [False, True], ids=["value", "result"])
@pytest.mark.parametrize("degree", [0, 1, 2, 3], ids=["mass", "odd", "quadratic", "cubic"])
@pytest.mark.parametrize(
    "domain, intervals",
    [
        pytest.param(sp.Interval(-1, 1), ((-1, 1),), id="closed"),
        pytest.param(sp.Interval.open(-1, 1), ((-1, 1),), id="open"),
        pytest.param(sp.FiniteSet(0), (), id="ambient-singleton"),
        pytest.param(sp.S.EmptySet, (), id="empty"),
        pytest.param(
            sp.Union(sp.Interval(-2, -1), sp.Interval(1, 3)), ((-2, -1), (1, 3)), id="disconnected"
        ),
    ],
)
def test_integral_matrix(domain, intervals, degree, structured):
    # The elementary antiderivative is independent of the region integration backend.
    expected = sum(
        (sp.Rational(b) ** (degree + 1) - sp.Rational(a) ** (degree + 1)) / (degree + 1)
        for a, b in intervals
    )
    result = integrate_over_region(
        X**degree, domain.as_relational(X), (X,), return_result=structured
    )
    assert (result.value if structured else result) == expected


@pytest.mark.parametrize(
    "expression, injective, surjective, smooth",
    [
        pytest.param(X, True, True, True, id="identity"),
        pytest.param(-2 * X + 3, True, True, True, id="affine"),
        pytest.param(X**2, False, False, True, id="square"),
        pytest.param(sp.Integer(2), False, False, True, id="constant"),
        pytest.param(sp.Abs(X), False, False, False, id="absolute-value"),
    ],
)
def test_mapping_matrix(expression, injective, surjective, smooth):
    assert is_injective(expression, X) is injective
    assert is_surjective(expression, X) is surjective
    assert is_bijective(expression, X) is (injective and surjective)
    assert is_function_continuous(expression, X) is True
    assert is_function_smooth(expression, X) is smooth


@pytest.mark.parametrize(
    "first", [-1, 0, 2], ids=["negative-first", "zero-first", "positive-first"]
)
@pytest.mark.parametrize(
    "second", [-2, 0, 1], ids=["negative-second", "zero-second", "positive-second"]
)
def test_diagonal_matrix_contracts(first, second):
    matrix = sp.diag(first, second)
    assert matrix_pd_on(matrix) is (first > 0 and second > 0)
    assert matrix_psd_on(matrix) is (first >= 0 and second >= 0)
    assert matrix_rank_on(matrix) == int(first != 0) + int(second != 0)
