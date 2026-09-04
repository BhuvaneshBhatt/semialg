import sympy as sp

from semialg import (
    convexity_certificate,
    function_convexity,
    function_mapping_properties,
    function_monotonicity,
    function_sign_partition,
)


def test_one_dimensional_convexity_detects_removed_point_gap():
    x = sp.symbols("x", real=True)
    result = convexity_certificate(sp.Ne(x, 0), (x,))
    assert result.outcome is False
    assert result.method == "one-dimensional-ordered-cad-interval"
    assert result.details["gap_count"] == 1


def test_reciprocal_is_not_globally_decreasing_across_disconnected_domain():
    x = sp.symbols("x", real=True)
    result = function_monotonicity(1 / x, x, return_result=True)
    assert result.classification == "nonmonotonic"
    assert result.method == "partition_cross_component_counterexample"
    assert result.counterexample is not None


def test_rational_sign_partition_uses_single_direct_cad():
    x = sp.symbols("x", real=True)
    result = function_sign_partition((x - 1) / (x + 2), x, return_result=True)
    assert result.method == "univariate_rational_sign_cad"
    assert [kind for kind, _ in result.pieces] == ["positive", "negative", "zero", "positive"]


def test_all_univariate_quasi_properties_share_one_monotonicity_partition():
    x = sp.symbols("x", real=True)
    result = function_convexity(x**3, (x,), properties="all", return_result=True)
    assert result.quasiconvex is True
    assert result.pseudoconvex is False
    assert result.details["analysis_cache"]["monotonic_partitions"] == 1


def test_convex_implications_avoid_duplicate_monotonicity_partitions():
    x = sp.symbols("x", real=True)
    result = function_convexity(x**4, (x,), properties="all", return_result=True)
    assert result.strictly_convex is True
    assert result.quasiconvex is True
    assert result.strictly_quasiconvex is True
    assert result.pseudoconvex is True
    assert result.details["analysis_cache"]["monotonic_partitions"] <= 1


def test_strong_convexity_uses_univariate_second_derivative_range():
    x = sp.symbols("x", real=True)
    result = function_convexity(x**4 + x**2, (x,), return_result=True)
    assert result.classification == "strongly_convex"
    assert result.strong_convexity_modulus == 2


def test_scalar_strict_monotonicity_shortcuts_injectivity():
    x = sp.symbols("x", real=True)
    result = function_mapping_properties(x**3, x)
    assert result.injective is True
    assert "strict_monotonicity_injectivity" in result.method
    assert "pairwise_collision_qe" not in result.method
