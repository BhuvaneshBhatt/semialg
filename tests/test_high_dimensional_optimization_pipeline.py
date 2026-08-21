import pytest
import sympy as sp

from semialg import (
    OptimizationCertificationPolicy,
    polynomial_locus_dimension,
    semialgebraic_minimize,
)
from semialg import optimization as optimization_module


@pytest.mark.slow
def test_cost_controlled_range_certifies_three_variable_open_ball():
    x, y, z = sp.symbols("x y z", real=True)
    result = semialgebraic_minimize(
        x,
        [x**2 + y**2 + z**2 < 1],
        [x, y, z],
        certification="auto",
        range_cost_limit=2500,
    )
    assert result.value == -1
    assert result.attained is False
    assert result.certified is True
    assert "cad_range_certificate" in result.method


def test_cost_policy_can_decline_expensive_range_fallback():
    x, y, z = sp.symbols("x y z", real=True)
    with pytest.raises(NotImplementedError):
        semialgebraic_minimize(
            x,
            [x**2 + y**2 + z**2 < 1],
            [x, y, z],
            certification="auto",
            range_cost_limit=1,
        )
    with pytest.raises(NotImplementedError):
        semialgebraic_minimize(
            x,
            [x**2 + y**2 + z**2 < 1],
            [x, y, z],
            certification="candidate",
        )


@pytest.mark.slow
def test_complete_policy_ignores_auto_cost_limit():
    x, y, z = sp.symbols("x y z", real=True)
    result = semialgebraic_minimize(
        x,
        [x**2 + y**2 + z**2 < 1],
        [x, y, z],
        certification="complete",
        range_cost_limit=0,
    )
    assert result.value == -1
    assert result.certified


def test_optimization_certification_policy_validates_inputs():
    with pytest.raises(ValueError):
        OptimizationCertificationPolicy("unknown")
    with pytest.raises(ValueError):
        OptimizationCertificationPolicy(range_cost_limit=-1)
    with pytest.raises(ValueError):
        OptimizationCertificationPolicy(recursion_limit=-1)


def test_polynomial_locus_dimension_detects_positive_zero_and_empty_loci():
    x, y, z = sp.symbols("x y z", real=True)
    assert polynomial_locus_dimension([], [x, y, z]) == 3
    assert polynomial_locus_dimension([x], [x, y, z]) == 2
    assert polynomial_locus_dimension([x * y], [x, y]) == 1
    assert polynomial_locus_dimension([x, y, z], [x, y, z]) == 0
    assert polynomial_locus_dimension([x, x - 1], [x]) == -1


@pytest.mark.slow
def test_positive_dimensional_kkt_locus_is_optimized_exactly():
    x, y = sp.symbols("x y", real=True)
    # The KKT locus contains the entire circle (lambda = 1), so finite-point
    # KKT solving alone is insufficient.  The positive-dimensional locus path
    # certifies the constant objective on that locus.
    result = semialgebraic_minimize(
        x**2 + y**2,
        [sp.Eq(x**2 + y**2, 1)],
        [x, y],
    )
    assert result.value == 1
    assert result.attained
    assert result.certified


def test_active_set_pruning_removes_scaled_duplicates_and_inconsistent_sets():
    x, y = sp.symbols("x y", real=True)
    subsets = optimization_module._pruned_active_subsets(
        (),
        (x, 2 * x, x - 1, y),
        (x, y),
    )
    # x and 2*x describe the same boundary, so no subset contains both.
    assert all(sum(sp.factor(item) in {x, 2 * x} for item in subset) <= 1 for subset in subsets)
    # x=0 and x-1=0 are algebraically inconsistent and must be pruned.
    assert not any(x in subset and x - 1 in subset for subset in subsets)


def test_equality_reduction_precedes_kkt_and_lifts_optimizer():
    x, y, z = sp.symbols("x y z", real=True)
    result = semialgebraic_minimize(
        x**2 + y**2 + z**2,
        [sp.Eq(x + y + z, 1)],
        [x, y, z],
    )
    assert result.value == sp.Rational(1, 3)
    assert result.point == {x: sp.Rational(1, 3), y: sp.Rational(1, 3), z: sp.Rational(1, 3)}
    assert result.method.startswith("equality_reduction+")
    assert result.certified


def test_equality_reduction_does_not_divide_by_variable_coefficient():
    x, y = sp.symbols("x y", real=True)
    objective, condition, variables, substitutions = optimization_module._reduce_linear_equalities(
        x**2 + y**2,
        sp.Eq(y * x + 1, 0),
        (x, y),
    )
    assert variables == (x, y)
    assert substitutions == {}
    assert objective == x**2 + y**2
    assert condition == sp.Eq(x * y + 1, 0)
