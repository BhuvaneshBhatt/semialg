import sympy as sp

from semialg.instances import (
    coordinate_bounds,
    eliminate_linear_equations,
    expand_with_subs,
    fast_factor_list,
    find_real_witnesses,
    is_bounded_solution_set,
    normalize_relations,
    sample_bounded_witnesses,
    satisfies_formula,
    try_fast_witness,
)
from semialg.solve import find_instance


def test_zero_rhs_normalization_uses_descriptive_api():
    x, y = sp.symbols("x y")
    normalized = normalize_relations(sp.And(x + 1 <= y, x > 0))
    assert sp.Le(x - y + 1, 0) in normalized.args
    assert sp.Lt(-x, 0) in normalized.args


def test_fast_real_heuristics_find_nonzero_polynomial_witness():
    x, y = sp.symbols("x y")
    result = try_fast_witness(sp.Ne(x**2 + y**2, 0), (x, y))
    assert result.found
    assert satisfies_formula(sp.Ne(x**2 + y**2, 0), result.instances[0])


def test_bounded_sampling_finds_disk_interior_point():
    x, y = sp.symbols("x y")
    result = sample_bounded_witnesses(x**2 + y**2 < 1, (x, y), seed=11)
    assert result.found
    assert satisfies_formula(x**2 + y**2 < 1, result.instances[0])


def test_real_fallback_pipeline_handles_decomposed_univariate_root():
    x = sp.symbols("x")
    result = find_real_witnesses(sp.Eq((x**2 - 2) ** 2, 0), (x,))
    assert result.found
    assert satisfies_formula(sp.Eq((x**2 - 2) ** 2, 0), result.instances[0])


def test_find_instance_auto_uses_fallback_diagnostics():
    x, y = sp.symbols("x y")
    result = find_instance(x**2 + y**2 < 1, [x, y], strategy="fallback", return_result=True)
    assert result.instances
    assert result.diagnostics["fallback_status"] in {"satisfied", "unknown"}
    assert satisfies_formula(x**2 + y**2 < 1, result.instances[0])


def test_equality_substitution_expands_conjunctions():
    x, y = sp.symbols("x y")
    simplified = expand_with_subs(sp.And(sp.Eq(x, y + 1), x > 3))
    assert sp.Gt(y + 1, 3) in simplified.args


def test_coordinate_bounds_and_boundedness_helpers():
    x, y = sp.symbols("x y")
    bounds = coordinate_bounds(sp.And(x >= -1, x <= 2, y > 0, y < 5), (x, y))
    assert not bounds.inconsistent
    assert is_bounded_solution_set(sp.And(x >= -1, x <= 2, y >= 0, y <= 5), (x, y)) is True


def test_linear_equation_elimination_helper():
    x, y = sp.symbols("x y")
    result = eliminate_linear_equations((x + y - 1, y - 2), (x, y), {x: 0, y: 2})
    assert len(result.replacements) >= 1
    assert result.variables != (x, y)


def test_fast_factor_list_avoids_expanding_large_power():
    x, y = sp.symbols("x y")
    factors = fast_factor_list((x + y) ** 200)
    assert ((x + y), 200) in factors


def test_coordinate_bounds_marks_unsupported_relations_incomplete():
    x, y = sp.symbols("x y", real=True)
    bounds = coordinate_bounds(x + y <= 1, (x, y))

    assert not bounds.complete
