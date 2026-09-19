from dataclasses import replace

import sympy as sp

from semialg import replay_certificate
from semialg.algebraic.gtz_zero_dim import (
    verify_zero_dimensional_primary_certificate,
    zero_dimensional_primary_decomposition,
)
from semialg.algebraic_function_fields import (
    MonogenicFunctionField,
    RationalFunctionField,
)


def test_gtz3_qq_repeated_components():
    x = sp.symbols("x")
    result = zero_dimensional_primary_decomposition(((x - 1) ** 2 * (x + 2),), (x,))
    assert result.complete
    assert result.irredundant
    assert sorted(component.degree for component in result.components) == [1, 2]
    assert verify_zero_dimensional_primary_certificate(result.certificate)
    assert replay_certificate(result).verified is True


def test_gtz3_quadratic_extension_splits_over_coefficient_field():
    x, a = sp.symbols("x a")
    qq = RationalFunctionField(())
    field = MonogenicFunctionField(qq, a, (-2, 0, 1))
    result = zero_dimensional_primary_decomposition(
        ((x**2 - 2) ** 2,), (x,), coefficient_field=field
    )
    assert result.complete
    assert len(result.components) == 2
    assert all(component.degree == 2 for component in result.components)
    radicals = [component.radical for component in result.components]
    assert any(any(sp.expand(g - (a - x)) == 0 for g in radical) for radical in radicals)
    assert any(any(sp.expand(g - (a + x)) == 0 for g in radical) for radical in radicals)
    assert verify_zero_dimensional_primary_certificate(result.certificate)


def test_gtz3_deep_algebraic_extension():
    x, a, b = sp.symbols("x a b")
    qq = RationalFunctionField(())
    k1 = MonogenicFunctionField(qq, a, (-2, 0, 1))
    k2 = MonogenicFunctionField(k1, b, (-k1.alpha, k1.zero, k1.one))
    result = zero_dimensional_primary_decomposition((x - b,), (x,), coefficient_field=k2)
    assert result.complete
    assert len(result.components) == 1
    assert result.components[0].degree == 1
    assert verify_zero_dimensional_primary_certificate(result.certificate)


def test_gtz3_rejects_reducible_coefficient_tower():
    x, a = sp.symbols("x a")
    qq = RationalFunctionField(())
    bad = MonogenicFunctionField(qq, a, (-1, 0, 1))
    try:
        zero_dimensional_primary_decomposition((x,), (x,), coefficient_field=bad)
    except ValueError as exc:
        assert "certified" in str(exc)
    else:
        raise AssertionError("reducible coefficient tower should be rejected")


def test_gtz3_tampered_separator_is_rejected():
    x, a = sp.symbols("x a")
    field = MonogenicFunctionField(RationalFunctionField(()), a, (-2, 0, 1))
    result = zero_dimensional_primary_decomposition(
        ((x**2 - 2) ** 2,), (x,), coefficient_field=field
    )
    certificate = result.certificate
    first = replace(certificate.components[0], separator=sp.Integer(1))
    forged = replace(certificate, components=(first, *certificate.components[1:]))
    assert not verify_zero_dimensional_primary_certificate(forged)


def test_gtz3_requires_zero_dimensional_source():
    x, y = sp.symbols("x y")
    try:
        zero_dimensional_primary_decomposition((x * y,), (x, y))
    except ValueError as exc:
        assert "zero-dimensional" in str(exc)
    else:
        raise AssertionError("positive-dimensional input should be rejected")
