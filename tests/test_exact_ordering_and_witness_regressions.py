from __future__ import annotations

from functools import cmp_to_key

import sympy as sp

from semialg._range_special_cases import _compare_bounds
from semialg._region_integrate_intrinsic import _zero_dimensional_points
from semialg.decision.sampling_helpers import _dedupe_samples
from semialg.simplify.boolean import simplify_boolean
from semialg.solve.domains import SolveDomain
from semialg.solve.find_instance import _make_result
from semialg.structural_keys import point_key


def test_exact_bound_comparison_distinguishes_close_rationals() -> None:
    lower = sp.Integer(1)
    upper = sp.Integer(1) + sp.Rational(1, 10**20)
    ordered = sorted((upper, lower), key=cmp_to_key(_compare_bounds))
    assert ordered == [lower, upper]


def test_dependent_witness_is_grounded_before_approximation() -> None:
    x, y = sp.symbols("x y", real=True)
    result = _make_result(
        instances=[{x: y, y: sp.Integer(1)}],
        variables=(x, y),
        domain=SolveDomain.REALS,
        method="regression",
    )
    assert result.instances[0] == {x: y, y: 1}
    assert result.approximate[0] == {x: 1.0, y: 1.0}


def test_structural_point_keys_keep_same_name_symbols_distinct() -> None:
    xr = sp.Symbol("x", real=True)
    xi = sp.Symbol("x", integer=True)
    assert xr != xi
    assert point_key({xr: 0}) != point_key({xi: 0})
    samples = _dedupe_samples(({xr: 0}, {xi: 0}))
    assert len(samples) == 2


def test_structural_point_keys_distinguish_swapped_same_name_coordinates() -> None:
    xr = sp.Symbol("x", real=True)
    xi = sp.Symbol("x", integer=True)
    assert point_key({xr: 0, xi: 1}) != point_key({xr: 1, xi: 0})


def test_boolean_simplifier_treats_relations_as_propositional_atoms() -> None:
    x = sp.Symbol("x", real=True)
    relation = sp.Eq(sp.sqrt(2) * x**2 + x - 1, 0)
    formula = sp.Or(relation, sp.And(relation, x > 0))
    assert simplify_boolean(formula) == relation


def test_zero_dimensional_reality_filter_rejects_complex_solutions_exactly() -> None:
    x, y = sp.symbols("x y", real=True)
    condition = sp.And(sp.Eq(x**2 + 1, 0, evaluate=False), sp.Eq(y, 0), evaluate=False)
    assert _zero_dimensional_points(condition, (x, y)) == ()


def test_nonsingular_affine_image_uses_exact_inverse_description() -> None:
    from semialg import affine_transform, is_compact, region_dimension

    x, y = sp.symbols("x y", real=True)
    box = sp.And(x >= -1, x <= 1, y >= -2, y <= 2)
    image = affine_transform(box, [[1, 1], [0, 2]], [3, -1], [x, y])
    expected = sp.And(y >= -5, y <= 3, 2 * x - y >= 5, 2 * x - y <= 9)

    assert sp.simplify_logic(sp.Equivalent(image, expected)) is sp.true
    assert region_dimension(image, [x, y]) == 2
    assert is_compact(image, [x, y])
