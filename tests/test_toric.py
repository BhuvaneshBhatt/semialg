import sympy as sp

from semialg.algebraic import (
    binomial_ideal,
    integer_kernel,
    lattice_ideal,
    markov_basis,
    toric_ideal,
)


def test_integer_kernel_is_saturated_lattice_basis():
    matrix = sp.Matrix([[2, 0, 1], [0, 2, 1]])
    basis = integer_kernel(matrix)
    assert basis == ((1, 1, -2),)
    assert matrix * sp.Matrix(basis[0]) == sp.zeros(2, 1)


def test_toric_ideal_for_two_by_two_independence():
    p = sp.symbols("p0:4")
    matrix = [[1, 1, 1, 1], [1, 1, 0, 0], [1, 0, 1, 0]]
    assert toric_ideal(matrix, p) == (p[0] * p[3] - p[1] * p[2],)
    assert markov_basis(matrix) == ((1, -1, -1, 1),)


def test_lattice_ideal_saturates_basis_binomials():
    x = sp.symbols("x0:3")
    moves = ((1, 1, -2),)
    assert binomial_ideal(moves, x) == (x[0] * x[1] - x[2] ** 2,)
    assert lattice_ideal(moves, x) == (x[0] * x[1] - x[2] ** 2,)


def test_toric_ideal_handles_negative_exponents():
    x0, x1 = sp.symbols("x0 x1")
    assert toric_ideal([[1, -1]], (x0, x1)) == (x0 * x1 - 1,)


def test_kernel_lattice_and_toric_ideal_agree_for_homogeneous_matrix():
    x = sp.symbols("x0:3")
    matrix = ((1, 1, 1), (0, 1, 2))
    kernel = integer_kernel(matrix)
    assert lattice_ideal(kernel, x) == toric_ideal(matrix, x)
    moves = markov_basis(matrix)
    points = {(0, 2, 0), (1, 0, 1)}
    assert any(
        tuple(a + sign * delta for a, delta in zip(point, move, strict=True)) in points
        for point in points
        for move in moves
        for sign in (-1, 1)
    )
