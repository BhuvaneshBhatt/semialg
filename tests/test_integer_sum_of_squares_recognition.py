from __future__ import annotations

import sympy as sp

from semialg.solve.integer.diophantine import detect_sum_eqn, solve_sum_of_two_squares
from semialg.solve.integer.families import detect_sum_fam
from semialg.solve.integer.special_families import (
    detect_diag_fam,
    detect_pell_family,
    detect_sum_fam2,
    solve_pell_family,
)


def test_sum_of_two_squares_recognizes_exact_supported_shape():
    x, y = sp.symbols("x y", integer=True)
    equation = sp.Eq(x**2 + y**2, 25)

    assert detect_sum_eqn(equation, (x, y)) == 25
    tag = detect_sum_fam(equation, (x, y))
    assert tag is not None
    assert tag.metadata["target_norm"] == 25
    special_tag = detect_sum_fam2(equation, (x, y))
    assert special_tag is not None
    assert special_tag.metadata["target"] == 25

    result = solve_sum_of_two_squares(equation, (x, y))
    assert result is not None
    assert result.complete is True
    assert (3, 4) in result.solutions


def test_sum_of_two_squares_rejects_extra_polynomial_terms():
    x, y = sp.symbols("x y", integer=True)
    unsupported = (
        sp.Eq(x**2 + y**2 + x, 25),
        sp.Eq(x**2 + y**2 + x * y, 25),
        sp.Eq(x**2 + y**2 + x**3, 25),
    )

    for equation in unsupported:
        assert detect_sum_eqn(equation, (x, y)) is None
        assert detect_sum_fam(equation, (x, y)) is None
        assert detect_sum_fam2(equation, (x, y)) is None
        assert solve_sum_of_two_squares(equation, (x, y)) is None


def test_sum_of_two_squares_rejects_nonpolynomial_and_symbolic_targets():
    x, y, n = sp.symbols("x y n", integer=True)
    unsupported = (
        sp.Eq(x**2 + y**2 + sp.sin(x), 25),
        sp.Eq(x**2 + y**2, n),
        sp.Eq(x**2 + y**2 + n, 25),
    )

    for equation in unsupported:
        assert detect_sum_eqn(equation, (x, y)) is None
        assert detect_sum_fam(equation, (x, y)) is None
        assert detect_sum_fam2(equation, (x, y)) is None


def test_related_diagonal_detector_rejects_unsupported_inputs_conservatively():
    x, y, n = sp.symbols("x y n", integer=True)

    assert detect_diag_fam(sp.Eq(x**2 + sp.sin(y), 1), (x, y)) is None
    assert detect_diag_fam(sp.Eq(x**2 + y**2, n), (x, y)) is None
    assert detect_diag_fam(sp.Eq(x**2, 1), (x, y)) is None


def test_pell_detector_and_solver_use_concrete_integer_parameters():
    x, y, d = sp.symbols("x y d", integer=True)
    pell = sp.Eq(x**2 - 2 * y**2, 1)

    desc = detect_pell_family(pell, (x, y))
    assert desc is not None
    assert desc.metadata == {"D": 2, "N": 1}
    result = solve_pell_family(pell, (x, y), search_bound=3)
    assert result is not None
    assert (3, 2) in result.solutions
    assert detect_pell_family(sp.Eq(x**2 - d * y**2, 1), (x, y)) is None
    assert detect_pell_family(sp.Eq(x**2 - 2 * y**2 + sp.sin(x), 1), (x, y)) is None
