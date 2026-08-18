from __future__ import annotations

import sympy as sp

from semialg import classify_real_roots


def test_classify_unparameterized_polynomial_roots_and_multiplicities() -> None:
    x = sp.symbols("x")
    result = classify_real_roots((x - 1) ** 2 * (x + 2), x)
    assert result.parameters == ()
    assert result.generic_root_count == 2
    assert result.generic_multiplicity_pattern == (1, 2)
    assert len(result.cells) == 1
    assert result.cells[0].condition is sp.true
    assert result.cells[0].root_count == 2


def test_classify_quadratic_parameter_family_by_discriminant() -> None:
    x, a, b = sp.symbols("x a b", real=True)
    result = classify_real_roots(x**2 + a * x + b, x, parameters=[a, b])
    conditions = {sp.sstr(cell.condition): cell.root_count for cell in result.cells}
    assert conditions["a**2 - 4*b > 0"] == 2
    assert conditions["Eq(a**2 - 4*b, 0)"] == 1
    assert conditions["a**2 - 4*b < 0"] == 0


def test_classify_quadratic_family_records_sample_parameters() -> None:
    x, a = sp.symbols("x a", real=True)
    result = classify_real_roots(x**2 + a, x, parameters=[a])
    by_count = {cell.root_count: cell for cell in result.cells}
    assert by_count[2].sample is not None
    assert bool((a < 0).subs(by_count[2].sample)) is True
    assert by_count[1].sample == {a: 0}
    assert bool((a > 0).subs(by_count[0].sample)) is True


def test_classify_linear_parameter_family_includes_degenerate_cases() -> None:
    x, a, b = sp.symbols("x a b", real=True)
    result = classify_real_roots(a * x + b, x, parameters=[a, b])
    conditions = {sp.sstr(cell.condition): cell.root_count for cell in result.cells}
    assert conditions["Ne(a, 0)"] == 1
    assert conditions["Eq(a, 0) & Eq(b, 0)"] == sp.oo
    assert conditions["Eq(a, 0) & Ne(b, 0)"] == 0
