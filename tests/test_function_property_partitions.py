import sympy as sp

from semialg import function_convex_partition, function_monotonic_partition
from semialg.conditional import ParameterStratifiedResult
from semialg.function_analysis import FunctionPropertyPartitionResult


def _classes(partition):
    return tuple(classification for classification, _ in partition)


def test_monotonic_partition_merges_isolated_derivative_zero():
    x = sp.symbols("x", real=True)
    partition = function_monotonic_partition(x**3, x)
    assert _classes(partition) == ("strictly_increasing",)


def test_monotonic_partition_finds_turning_points():
    x = sp.symbols("x", real=True)
    partition = function_monotonic_partition(x**3 - 3 * x, x)
    assert _classes(partition) == (
        "strictly_increasing",
        "strictly_decreasing",
        "strictly_increasing",
    )


def test_monotonic_partition_preserves_disconnected_domain():
    x = sp.symbols("x", real=True)
    partition = function_monotonic_partition(1 / x, x)
    assert _classes(partition) == ("strictly_decreasing", "strictly_decreasing")
    assert sp.simplify_logic(sp.Xor(partition[0][1], partition[1][1])) != sp.false


def test_monotonic_partition_handles_abs_graph_derivative():
    x = sp.symbols("x", real=True)
    partition = function_monotonic_partition(sp.Abs(x), x)
    assert _classes(partition) == ("strictly_decreasing", "strictly_increasing")


def test_convex_partition_exact_algebraic_inflection_points():
    x = sp.symbols("x", real=True)
    partition = function_convex_partition(x**4 - x**2, x)
    assert _classes(partition) == ("convex", "concave", "convex")
    boundary = sp.sqrt(6) / 6
    assert any(region.has(boundary) for _, region in partition)


def test_convex_partition_disconnected_rational_domain():
    x = sp.symbols("x", real=True)
    partition = function_convex_partition(1 / x, x)
    assert _classes(partition) == ("concave", "convex")


def test_partition_result_mode():
    x = sp.symbols("x", real=True)
    result = function_monotonic_partition(x**2, x, return_result=True)
    assert isinstance(result, FunctionPropertyPartitionResult)
    assert result.certified
    assert result.property_name == "monotonicity"


def test_parameterized_partitions_use_parameter_first_cad():
    x, a = sp.symbols("x a", real=True)
    monotonic = function_monotonic_partition(x**3 + a * x, x)
    assert isinstance(monotonic, ParameterStratifiedResult)
    negative = monotonic.select({a: -3})
    assert tuple(value for value, _ in negative) == (
        "strictly_increasing",
        "constant",
        "strictly_decreasing",
        "constant",
        "strictly_increasing",
    )
    convex = function_convex_partition(a * x**2, x)
    assert convex.select({a: 1})[0][0] == "convex"
    assert convex.select({a: 0})[0][0] == "affine"
    assert convex.select({a: -1})[0][0] == "concave"


def test_function_sign_partition_polynomial_components():
    from semialg import function_sign_partition

    x = sp.symbols("x", real=True)
    partition = function_sign_partition(x * (x - 1), x)
    assert _classes(partition) == ("positive", "zero", "negative", "zero", "positive")


def test_function_sign_partition_supported_nonsmooth_graph():
    from semialg import function_sign_partition

    x = sp.symbols("x", real=True)
    partition = function_sign_partition(sp.Abs(x), x)
    assert _classes(partition) == ("positive", "zero", "positive")


def test_function_sign_partition_preserves_disconnected_natural_domain():
    from semialg import function_sign_partition

    x = sp.symbols("x", real=True)
    partition = function_sign_partition(1 / x, x)
    assert _classes(partition) == ("negative", "positive")


def test_function_sign_partition_parameter_conditions():
    from semialg import function_sign_partition

    x, a = sp.symbols("x a", real=True)
    result = function_sign_partition(a * x, x)
    assert isinstance(result, ParameterStratifiedResult)
    assert tuple(value for value, _ in result.select({a: 1})) == ("negative", "zero", "positive")
    assert tuple(value for value, _ in result.select({a: -1})) == ("positive", "zero", "negative")
    assert all(value == "zero" for value, _ in result.select({a: 0}))
