import pytest
import sympy as sp

from semialg.applications import geometric_probability, polynomial_probability


def test_polynomial_probability_normalizes_density_exactly():
    x = sp.Symbol("x", real=True)
    result = polynomial_probability(
        x <= sp.Rational(1, 2),
        [x],
        support=sp.true,
        density=2 * x,
        bounds={x: (0, 1)},
    )

    assert result.certified
    assert result.normalizing_mass == 1
    assert result.event_mass == sp.Rational(1, 4)
    assert result.probability == sp.Rational(1, 4)


def test_geometric_probability_is_exact_uniform_measure_ratio():
    x = sp.Symbol("x")
    result = geometric_probability(
        x <= sp.Rational(1, 2),
        ["x"],
        support=sp.And(x >= 0, x <= 1),
    )
    assert result.variables == (x,)
    assert result.probability == sp.Rational(1, 2)


def test_probability_rejects_negative_density_and_zero_mass_support():
    x = sp.Symbol("x", real=True)
    with pytest.raises(ValueError, match="nonnegative"):
        polynomial_probability(x >= 0, [x], support=sp.And(x >= 0, x <= 1), density=-1)
    with pytest.raises(ValueError, match="positive normalizing mass"):
        polynomial_probability(sp.true, [x], support=sp.Eq(x, 0), density=1)


def test_probability_rejects_nonpolynomial_density_and_undeclared_parameters():
    x, a = sp.symbols("x a", real=True)
    with pytest.raises(ValueError, match="polynomial"):
        polynomial_probability(x >= 0, [x], support=sp.And(x >= 0, x <= 1), density=sp.exp(x))
    with pytest.raises(ValueError, match="undeclared"):
        polynomial_probability(x >= 0, [x], support=x <= a, density=1)
