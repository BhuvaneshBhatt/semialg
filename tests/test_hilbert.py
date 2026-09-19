import sympy as sp

from semialg.algebraic import (
    hilbert_function,
    hilbert_polynomial,
    ideal_degree,
    ideal_hilbert_data,
)


def test_hilbert_hypersurface_degree_and_polynomial():
    x, y = sp.symbols("x y")
    data = ideal_hilbert_data((x**3 + y**3,), (x, y))
    assert data.dimension == 1
    assert data.degree == 3
    assert hilbert_polynomial((x**3 + y**3,), (x, y)) == 3


def test_hilbert_zero_dimensional_length():
    x, y = sp.symbols("x y")
    equations = (x**2, y**3)
    assert ideal_degree(equations, (x, y)) == 6
    assert [hilbert_function(equations, (x, y), d) for d in range(6)] == [1, 2, 2, 1, 0, 0]
    assert hilbert_polynomial(equations, (x, y)) == 0


def test_hilbert_affine_space():
    x, y, z = sp.symbols("x y z")
    data = ideal_hilbert_data((), (x, y, z))
    assert data.dimension == 3
    assert data.degree == 1
    d = sp.Symbol("d")
    assert sp.expand(hilbert_polynomial((), (x, y, z), symbol=d) - (d + 1) * (d + 2) / 2) == 0


def test_hilbert_degree_of_two_components():
    x, y = sp.symbols("x y")
    assert ideal_degree((x * y,), (x, y)) == 2
