from __future__ import annotations

import sympy as sp

from semialg import sign_at, sign_vector
from semialg.algebraic import compute_rational_univariate_representation, solve_rur_points


def test_sign_at_exact_rootof_expression() -> None:
    x = sp.symbols("x")
    root = sp.RootOf(x**2 - 2, 1)
    assert sign_at(x**2 - 2, {x: root}) == 0
    assert sign_at(x - 1, {x: root}) == 1
    assert sign_at(x - 2, {x: root}) == -1


def test_sign_at_rur_point_reduces_expression_to_parameter() -> None:
    x, y = sp.symbols("x y")
    representation = compute_rational_univariate_representation([x**2 + y**2 - 1, x - y], [x, y])
    points = solve_rur_points(representation)
    positive = next(point for point in points if sign_at(x, point, variables=[x, y]) > 0)
    assert sign_at(x - y, positive, variables=[x, y]) == 0
    assert sign_at(x + y, positive, variables=[x, y]) == 1
    assert sign_at(x**2 + y**2 - 1, positive, variables=[x, y]) == 0


def test_sign_vector_accepts_rur_point_and_dict_output() -> None:
    x, y = sp.symbols("x y")
    representation = compute_rational_univariate_representation([x**2 + y**2 - 1, x - y], [x, y])
    point = next(
        point
        for point in solve_rur_points(representation)
        if sign_at(x, point, variables=[x, y]) > 0
    )
    signs = sign_vector([x - y, x + y, x**2 + y**2 - 1], point, variables=[x, y], as_dict=True)
    assert signs[x - y] == 0
    assert signs[x + y] == 1
    assert signs[x**2 + y**2 - 1] == 0
