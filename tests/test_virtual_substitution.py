import pytest
import sympy as sp

from semialg.algebraic.rational_univariate import evaluate_boolean_formula_at_point
from semialg.qe.virtual_substitution import (
    QuadraticVirtualSubstitutionResult,
    VirtualSubstitutionError,
    eliminate_exists_quadratic_variable,
    eliminate_quadratic_variable,
)


def _truth(formula, **assignment):
    return evaluate_boolean_formula_at_point(formula, assignment)


def test_eliminates_quadratic_disk_projection_without_the_eliminated_variable():
    x, y = sp.symbols("x y")

    result = eliminate_quadratic_variable(x**2 + y**2 <= 1, x, simplify=False)

    assert x not in result.free_symbols
    assert _truth(result, y=0)
    assert _truth(result, y=1)
    assert not _truth(result, y=2)


def test_eliminates_open_half_disk_projection():
    x, y = sp.symbols("x y")

    result = eliminate_quadratic_variable(sp.And(x > 0, x**2 + y**2 <= 1), x, simplify=False)

    assert x not in result.free_symbols
    assert _truth(result, y=0)
    assert not _truth(result, y=1)
    assert not _truth(result, y=-1)
    assert not _truth(result, y=2)


def test_linear_equality_degeneracy_is_preserved():
    x, a, b = sp.symbols("x a b")

    result = eliminate_quadratic_variable(sp.Eq(a * x + b, 0), x, simplify=False)

    assert x not in result.free_symbols
    assert _truth(result, a=2, b=3)
    assert _truth(result, a=0, b=0)
    assert not _truth(result, a=0, b=1)


def test_quadratic_equality_with_positive_root_side_condition():
    x, y = sp.symbols("x y")

    result = eliminate_quadratic_variable(sp.And(sp.Eq(x**2 - y, 0), x > 0), x, simplify=False)

    assert x not in result.free_symbols
    assert _truth(result, y=1)
    assert not _truth(result, y=0)
    assert not _truth(result, y=-1)


def test_interval_with_strict_bounds_uses_epsilon_candidates():
    x = sp.symbols("x")

    result = eliminate_quadratic_variable(sp.And(x > 1, x < 2), x, simplify=False)

    assert result == sp.true or _truth(result)


def test_structured_result_records_backend_and_eliminated_variable():
    x, y = sp.symbols("x y")

    result = eliminate_exists_quadratic_variable(x**2 + y <= 0, x, simplify=False)

    assert isinstance(result, QuadraticVirtualSubstitutionResult)
    assert result.eliminated_variable == x
    assert result.backend == "quadratic-virtual-substitution"
    assert x not in result.formula.free_symbols
    assert _truth(result.formula, y=-1)
    assert not _truth(result.formula, y=1)


def test_rejects_atoms_above_quadratic_degree():
    x = sp.symbols("x")

    with pytest.raises(VirtualSubstitutionError):
        eliminate_quadratic_variable(x**3 - 1 <= 0, x)
