"""Norms compared with independent determinant/resultant computations."""

import random

import pytest
import sympy as sp

from semialg._function_field_factorization import _det_poly
from semialg.algebraic_function_fields import (
    MonogenicFunctionField,
    RationalFunctionField,
    _norm_polynomial,
)


def test_polynomial_determinants_match_sympy_for_singular_and_pivoted_matrices():
    rng = random.Random(913)
    t, z = sp.symbols("t z")
    field = RationalFunctionField((z,))
    matrices = [
        [],
        [[[1, z]]],
        [[[], [1]], [[1], []]],
        [[[], [1], []], [[], [], [1]], [[1], [], []]],
    ]
    for n in (2, 3, 4):
        matrices.append(
            [
                [
                    [sp.Integer(rng.randrange(-2, 3)), sp.Integer(rng.randrange(-1, 2))]
                    for _ in range(n)
                ]
                for _ in range(n)
            ]
        )
    matrices.append([[[1], [2]], [[2], [4]]])
    for matrix in matrices:
        coeffs = _det_poly(matrix, field)
        actual = sum(c * t**i for i, c in enumerate(coeffs))
        expected = sp.Matrix(
            [[sum(c * t**i for i, c in enumerate(p)) for p in row] for row in matrix]
        ).det()
        assert sp.cancel(actual - expected) == 0
    with pytest.raises(ValueError, match="square"):
        _det_poly([[[1], [2]]], field)


@pytest.mark.parametrize("degree", [2, 3, 5, 8])
def test_shifted_generator_norm_matches_resultant(degree):
    a, z, t = sp.symbols("a z t")
    base = RationalFunctionField((z,))
    field = MonogenicFunctionField(
        base, a, (-z,) + (sp.Integer(0),) * (degree - 1) + (sp.Integer(1),)
    )
    coefficients = _norm_polynomial((-(field.alpha + 1), field.one), field)
    actual = sum(c * t**i for i, c in enumerate(coefficients))
    expected = sp.resultant(a**degree - z, t - a - 1, a)
    assert sp.expand(actual - expected) == 0
    assert _norm_polynomial((field.zero,), field) == []


def test_norm_over_a_monogenic_base_preserves_tower_elements():
    a, b, z = sp.symbols("a b z")
    lower = MonogenicFunctionField(RationalFunctionField((z,)), a, (-z, 0, 1))
    top = MonogenicFunctionField(lower, b, (-lower.alpha, lower.zero, lower.one))
    result = _norm_polynomial((-top.alpha, top.one), top)
    expected = (-lower.alpha, lower.zero, lower.one)
    assert len(result) == len(expected)
    assert all(lower.is_zero(c - d) for c, d in zip(result, expected, strict=True))
