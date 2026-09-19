from __future__ import annotations

import sympy as sp

from semialg import (
    BallRegion,
    BoxRegion,
    ParametricRegion,
    SemialgebraicRegion,
    as_semialgebraic_region,
    region_closure,
    region_difference,
    region_dimension,
    region_interior,
    region_intersection,
    region_union,
)
from semialg.reasoning import region_disjoint, region_subset
from semialg.symbolic_regions import (
    RDisjoint,
    RegionElement,
    RegionNotElement,
    REqual,
    RSubset,
    region_element_conditions,
    region_relation_conditions,
)


def test_formula_region_carries_lazy_context_and_reusable_cad():
    x = sp.Symbol("x", real=True)
    region = SemialgebraicRegion(sp.And(x >= 0, x <= 1), (x,))

    assert region.formula == sp.And(x >= 0, x <= 1)
    assert region.variables == (x,)
    assert region.cad is None
    assert region.context is region.context

    first = region.ensure_cad()
    second = region.ensure_cad()
    assert first is second
    assert region.cad is first
    assert region.contains((sp.Rational(1, 2),)) is True
    assert region.contains((2,)) is False


def test_explicit_regions_lower_to_same_symbolic_region_layer():
    x, y = sp.symbols("x y", real=True)
    box = BoxRegion(((0, 1), (-1, 2)))
    ball = BallRegion((0, 0), 1)

    symbolic_box = box.as_semialgebraic_region((x, y))
    symbolic_ball = as_semialgebraic_region(ball, (x, y))

    assert isinstance(symbolic_box, SemialgebraicRegion)
    assert symbolic_box.formula == sp.And(x >= 0, x <= 1, y >= -1, y <= 2)
    assert symbolic_ball.formula == (x**2 + y**2 <= 1)


def test_membership_predicates_lower_inside_boolean_expressions():
    x, y = sp.symbols("x y", real=True)
    region = as_semialgebraic_region(BoxRegion(((0, 1), (0, 1))), (x, y))
    point = (sp.Rational(1, 2), y)

    element = RegionElement(point, region)
    not_element = RegionNotElement((2, 2), region)
    condition = region_element_conditions(sp.And(element, not_element))

    assert sp.simplify(condition) == sp.And(y >= 0, y <= 1)
    assert RegionElement((sp.Rational(1, 2), sp.Rational(1, 2)), region).evaluate() is True
    assert not_element.evaluate() is True


def test_region_relation_nodes_lower_to_quantified_formulas_and_evaluate():
    x = sp.Symbol("x", real=True)
    inner = SemialgebraicRegion(sp.And(x >= 0, x <= 1), (x,))
    outer = SemialgebraicRegion(sp.And(x >= -1, x <= 2), (x,))
    disjoint = SemialgebraicRegion(sp.And(x >= 3, x <= 4), (x,))

    subset = RSubset(inner, outer)
    assert isinstance(subset.as_formula(), sp.Basic)
    assert subset.evaluate() is True
    assert RDisjoint(inner, disjoint).evaluate() is True
    assert REqual(inner, inner).evaluate() is True

    lowered = region_relation_conditions(subset)
    assert "ForAll" in str(lowered)
    assert region_relation_conditions(subset, eliminate=True) is sp.true


def test_region_operations_accept_unified_and_explicit_regions():
    x = sp.Symbol("x", real=True)
    left = SemialgebraicRegion(sp.And(x >= 0, x <= 1), (x,))
    right = BoxRegion(((1, 2),))

    union = region_union(left, right)
    intersection = region_intersection(left, right)
    difference = region_difference(left, right)

    assert isinstance(union, SemialgebraicRegion)
    assert isinstance(intersection, SemialgebraicRegion)
    assert isinstance(difference, SemialgebraicRegion)
    assert union.contains((sp.Rational(3, 2),)) is True
    assert intersection.contains((1,)) is True
    assert difference.contains((sp.Rational(1, 2),)) is True
    assert difference.contains((1,)) is False

    assert region_subset(left, union) is True
    assert region_disjoint(difference, as_semialgebraic_region(BoxRegion(((1, 2),)), (x,))) is True


def test_parametric_region_keeps_existential_membership_until_needed():
    t = sp.Symbol("t", real=True)
    x, y = sp.symbols("x y", real=True)
    explicit = ParametricRegion((t,), ((t, 0, 1),), (t, t**2))
    region = as_semialgebraic_region(explicit, (x, y))

    assert "Exists" in str(region.formula)
    symbolic_membership = region.membership_formula((x, y), eliminate=False)
    assert "Exists" in str(symbolic_membership)
    assert region.contains((0, 0)) is True
    assert region.contains((1, 1)) is True


def test_topology_helpers_preserve_unified_region_objects():
    x = sp.Symbol("x", real=True)
    open_interval = SemialgebraicRegion(sp.And(x > 0, x < 1), (x,))

    closure = region_closure(open_interval, strategy="syntactic")
    interior = region_interior(closure, strategy="syntactic")

    assert isinstance(closure, SemialgebraicRegion)
    assert isinstance(interior, SemialgebraicRegion)
    assert closure.formula == sp.And(x >= 0, x <= 1)
    assert interior.formula == sp.And(x > 0, x < 1)
    assert region_dimension(closure) == 1
    assert closure.is_bounded() is True
    assert closure.is_closed(strategy="syntactic") is True
    assert closure.is_compact(strategy="syntactic") is True


def test_nested_quantified_polygon_branches_are_eliminated_for_membership():
    from semialg import PolygonRegion

    x, y = sp.symbols("x y", real=True)
    polygon = PolygonRegion(((0, 0), (3, 0), (3, 3), (2, 3), (2, 1), (1, 1), (1, 3), (0, 3)))
    region = as_semialgebraic_region(polygon, (x, y))

    qf = region.quantifier_free_formula()
    assert "Exists" not in str(qf)
    assert region.contains((sp.Rational(1, 2), 2)) is True
    assert region.contains((sp.Rational(3, 2), 2)) is False


def test_parametric_membership_survives_elimination():
    t = sp.Symbol("t", real=True)
    x, y = sp.symbols("x y", real=True)
    region = as_semialgebraic_region(ParametricRegion((t,), ((t, 0, 1),), (t, t**2)), (x, y))
    atom = RegionElement((x, y), region)

    lowered = region_element_conditions(atom, eliminate=False)
    assert lowered == atom
    assert region_element_conditions(atom, eliminate=True) != atom


def test_region_condition_real_parameters_adds_explicit_membership():
    x = sp.Symbol("x", real=True)
    a = sp.Symbol("a")
    region = SemialgebraicRegion(sp.And(x >= 0, x <= a), (x,))
    atom = RegionElement((x,), region)

    lowered = region_element_conditions(atom, eliminate=False, real_parameters=True)
    assert sp.Contains(a, sp.S.Reals, evaluate=False) in lowered.args

    relation = region_relation_conditions(RSubset(region, region), real_parameters=True)
    assert relation.has(sp.Contains(a, sp.S.Reals, evaluate=False))
