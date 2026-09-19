from __future__ import annotations

import sympy as sp

from semialg.algebraic import (
    compute_rational_univariate_representation,
    sign_conditions_at_rur_points,
    solve_rur_points,
    thom_encoding,
    thom_encodings,
)
from semialg.algebraic.thom import ThomEncoding


def test_thom_encodings_distinguish_real_roots():
    x = sp.symbols("x", real=True)
    encodings = thom_encodings((x + 2) ** 2 * (x - 1), x)

    assert tuple(item.root for item in encodings) == (-2, 1)
    assert tuple(item.multiplicity for item in encodings) == (2, 1)
    assert encodings[0].derivative_signs == (0, -1, 1)
    assert encodings[1].derivative_signs == (1, 1, 1)


def test_thom_sign_determination_at_algebraic_root():
    x = sp.symbols("x", real=True)
    root = sp.sqrt(2)
    encoding = thom_encoding(x**2 - 2, root, x)

    assert isinstance(encoding, ThomEncoding)
    assert encoding.sign_of(x - 1) == 1
    assert encoding.sign_of(x**2 - 2) == 0
    assert encoding.sign_of(-x - 1) == -1


def test_rur_points_carry_thom_root_certificates():
    x, y = sp.symbols("x y", real=True)
    representation = compute_rational_univariate_representation([x**2 - 2, y - x], (x, y))
    points = solve_rur_points(representation)

    assert tuple(point.coordinates for point in points) == (
        (-sp.sqrt(2), -sp.sqrt(2)),
        (sp.sqrt(2), sp.sqrt(2)),
    )
    assert tuple(point.thom_encoding.derivative_signs[0] for point in points) == (-1, 1)
    assert tuple(point.sign_of(x + y) for point in points) == (-1, 1)


def test_rur_sign_conditions_share_univariate_certificate():
    x, y = sp.symbols("x y", real=True)
    representation = compute_rational_univariate_representation([x**2 - 2, y - x], (x, y))
    points = solve_rur_points(representation)

    signs = sign_conditions_at_rur_points((x, y, x + y, x**2 - 2), points)

    assert signs == ((-1, -1, -1, 0), (1, 1, 1, 0))
