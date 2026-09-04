from __future__ import annotations

import pytest
import sympy as sp
from hypothesis import given, settings
from hypothesis import strategies as st

from semialg import is_equal, root_count_conditions

pytestmark = pytest.mark.slow

x, y, u = sp.symbols("x y u", real=True)


@given(
    st.integers(min_value=1, max_value=7),
    st.integers(min_value=-5, max_value=5),
)
@settings(max_examples=12, deadline=None)
def test_positive_scaling_preserves_polynomial_inequality(scale: int, offset: int) -> None:
    lhs = x**2 + offset * x - 2
    assert is_equal(lhs <= 0, sp.expand(scale * lhs) <= 0, [x])


@given(st.permutations((x >= -2, x <= 3, x**2 <= 4)))
@settings(max_examples=6, deadline=None)
def test_conjunction_permutation_is_semantically_invariant(atoms) -> None:
    original = sp.And(x >= -2, x <= 3, x**2 <= 4)
    assert is_equal(original, sp.And(*atoms), [x])


@given(st.integers(min_value=-4, max_value=4))
@settings(max_examples=9, deadline=None)
def test_parameter_root_count_is_variable_renaming_invariant(a_value: int) -> None:
    a = sp.Symbol("a", real=True)
    b = sp.Symbol("b", real=True)
    first = root_count_conditions(x**5 + a * x + 1, x, [a])
    second = root_count_conditions(u**5 + b * u + 1, u, [b])

    selected_first = [
        count
        for count, condition in first.items()
        if sp.simplify(condition.subs(a, a_value)) is sp.true
    ]
    selected_second = [
        count
        for count, condition in second.items()
        if sp.simplify(condition.subs(b, a_value)) is sp.true
    ]
    assert selected_first == selected_second
