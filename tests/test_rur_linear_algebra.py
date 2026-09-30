import sympy as sp

from semialg.algebraic.rational_univariate.linear_algebra import (
    krylov_minimal_polynomial,
    matrix_density,
    preferred_matrix_storage,
)


def test_krylov_polynomial_matches_cyclic_characteristic_polynomial():
    t = sp.Symbol("t")
    matrix = sp.Matrix([[0, 2], [1, 0]])
    assert krylov_minimal_polynomial(matrix, t).as_expr() == t**2 - 2


def test_matrix_storage_uses_density():
    dense = sp.eye(2)
    sparse = sp.diag(1, 0, 0, 0)
    assert matrix_density(dense) == sp.Rational(1, 2)
    assert preferred_matrix_storage(dense) == "dense"
    assert preferred_matrix_storage(sparse) == "sparse"
