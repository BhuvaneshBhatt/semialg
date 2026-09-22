import sympy as sp

import semialg
from semialg import SemialgebraicRegion
from semialg._api_policy import EXPERT_EXPORTS, PRIMARY_EXPORTS


def test_region_relation_api_has_one_primary_decision_vocabulary():
    assert {"is_subset", "is_equal", "is_disjoint", "is_interior_disjoint", "intersects"} <= set(
        PRIMARY_EXPORTS
    )
    assert {"region_subset", "region_equal", "region_disjoint"} <= set(EXPERT_EXPORTS)
    assert {"RegionSubset", "RegionEqual", "RegionDisjoint"} <= set(EXPERT_EXPORTS)


def test_region_symmetric_difference_formula_and_region_object_agree():
    x = sp.Symbol("x", real=True)
    left = sp.And(x >= 0, x <= 2)
    right = sp.And(x >= 1, x <= 3)
    formula = semialg.region_symmetric_difference(left, right)
    assert semialg.is_equal(formula, sp.Or(sp.And(x >= 0, x < 1), sp.And(x > 2, x <= 3)), (x,))

    region = semialg.region_symmetric_difference(
        SemialgebraicRegion(left, (x,)), SemialgebraicRegion(right, (x,))
    )
    assert isinstance(region, SemialgebraicRegion)
    assert semialg.is_equal(region, formula, (x,))


def test_is_interior_disjoint_allows_boundary_contact():
    x = sp.Symbol("x", real=True)
    left = sp.And(x >= 0, x <= 1)
    touching = sp.And(x >= 1, x <= 2)
    overlapping = sp.And(x >= sp.Rational(1, 2), x <= 2)

    assert not semialg.is_disjoint(left, touching, (x,))
    assert semialg.is_interior_disjoint(left, touching, (x,))
    assert not semialg.is_interior_disjoint(left, overlapping, (x,))


def test_is_interior_disjoint_handles_lower_dimensional_regions():
    x, y = sp.symbols("x y", real=True)
    line_a = sp.And(sp.Eq(y, 0), x >= 0, x <= 1)
    line_b = sp.And(sp.Eq(y, 0), x >= sp.Rational(1, 2), x <= 2)
    assert semialg.is_interior_disjoint(line_a, line_b, (x, y))


def test_connected_components_is_the_only_primary_component_api():
    assert "connected_components" in PRIMARY_EXPORTS
    assert "region_components" not in PRIMARY_EXPORTS
    assert "explicit_region_components" in EXPERT_EXPORTS
    assert not hasattr(semialg, "region_components")


def test_semialgebraic_region_components_use_certified_connectivity():
    x, y = sp.symbols("x y", real=True)
    formula = sp.Or(
        sp.And((x + 2) ** 2 + y**2 <= 1),
        sp.And((x - 2) ** 2 + y**2 <= 1),
    )
    region = SemialgebraicRegion(formula, (x, y))
    components = region.components()
    assert len(components) == 2
    assert all(isinstance(component, SemialgebraicRegion) for component in components)


def test_normalized_moment_root_api_has_no_region_prefixed_duplicates():
    assert "centroid" in PRIMARY_EXPORTS
    assert "covariance_matrix" in PRIMARY_EXPORTS
    assert "region_centroid" not in PRIMARY_EXPORTS
    assert "region_covariance" not in PRIMARY_EXPORTS
    assert not hasattr(semialg, "region_centroid")
    assert not hasattr(semialg, "region_covariance")


def test_completed_root_naming_audit_has_one_image_vocabulary():
    assert {"affine_image", "affine_preimage", "region_image", "region_preimage"} <= set(
        PRIMARY_EXPORTS
    )
    assert {"semialgebraic_image", "semialgebraic_preimage"} <= set(EXPERT_EXPORTS)
    assert "affine_transform" not in PRIMARY_EXPORTS


def test_boolean_region_structured_constructors_are_class_methods():
    from semialg import BooleanRegion, Interval

    left = Interval(0, 1)
    right = Interval(1, 2)
    assert isinstance(BooleanRegion.union(left, right), BooleanRegion)
    assert isinstance(BooleanRegion.intersection(left, right), BooleanRegion)
    assert isinstance(BooleanRegion.difference(left, right), BooleanRegion)
    assert isinstance(BooleanRegion.symmetric_difference(left, right), BooleanRegion)
    for old_name in (
        "RegionUnion",
        "RegionIntersection",
        "RegionDifference",
        "RegionSymmetricDifference",
    ):
        assert old_name not in PRIMARY_EXPORTS


def test_regularization_composition_names_are_explicit():
    assert {"closure_of_interior", "interior_of_closure"} <= set(PRIMARY_EXPORTS)
    assert "region_interior_closure" not in PRIMARY_EXPORTS
    assert "region_closure_interior" not in PRIMARY_EXPORTS
