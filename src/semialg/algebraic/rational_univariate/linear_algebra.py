"""Exact linear-algebra kernels used by rational univariate representations."""

from algroots.rational_univariate.linear_algebra import (
    krylov_minimal_polynomial,
    matrix_density,
    packed_krylov_coefficients,
    packed_krylov_minimal_polynomial,
    preferred_matrix_storage,
)

__all__ = [
    "krylov_minimal_polynomial",
    "packed_krylov_minimal_polynomial",
    "packed_krylov_coefficients",
    "matrix_density",
    "preferred_matrix_storage",
]
