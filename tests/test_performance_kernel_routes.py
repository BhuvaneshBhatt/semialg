import pytest
import sympy as sp
from sympy.polys.domains import GF, QQ
from sympy.polys.matrices import DomainMatrix

from semialg.algebraic.rational_univariate.linear_algebra import (
    packed_krylov_coefficients,
    packed_krylov_minimal_polynomial,
)
from semialg.heuristics import rank_variables_by_polynomials, suggest_cad_variable_order


def test_packed_krylov_matches_charpoly():
    t = sp.Symbol("t")
    matrix = sp.Matrix([[0, 0, -6], [1, 0, -11], [0, 1, -6]])
    assert (
        sp.expand(
            packed_krylov_minimal_polynomial(matrix, t).as_expr() - matrix.charpoly(t).as_expr()
        )
        == 0
    )


def test_packed_krylov_rejects_noncyclic_seed():
    t = sp.Symbol("t")
    with pytest.raises(ValueError):
        packed_krylov_minimal_polynomial(sp.eye(3), t)


def test_auto_order_is_brown_or_chordal_candidate():
    a, b, c, d = sp.symbols("a b c d")
    polys = (a * b - 1, b * c - 1, c * d - 1)
    score = suggest_cad_variable_order(polys, (a, b, c, d), strategy="auto")
    assert score.strategy in {"brown", "chordal"}
    assert score.order in {
        rank_variables_by_polynomials(polys, (a, b, c, d), strategy="brown"),
        rank_variables_by_polynomials(polys, (a, b, c, d), strategy="chordal"),
    }


def test_domain_packed_krylov_matches_charpoly():
    for domain in (QQ, GF(65521)):
        rows = [
            [domain.zero, domain.zero, -domain.one],
            [domain.one, domain.zero, -domain.one],
            [domain.zero, domain.one, -domain.one],
        ]
        packed = packed_krylov_coefficients(rows, domain)
        expected_descending = DomainMatrix.from_list(rows, domain).charpoly()
        assert tuple(reversed(packed)) == tuple(expected_descending)


def test_auto_order_handles_degree_growth_miss():
    a, x, y = sp.symbols("a x y")
    polys = (
        x - 1,
        2 - a,
        1 + x,
        9 - x,
        y + 6,
        6 - y,
        -(a**2 + 1) * (-y - a * x + 2 * a**2 + 2),
        (a**2 + 1) ** 2 - (-x + a * y) * (a**2 + 1),
        a**2 + 1,
        (a**2 + 1) * (a**2 + 1 - x + a * y),
        4 * (a**2 + 1) ** 2 - (a**2 + 1) * (y + a * x),
    )
    score = suggest_cad_variable_order(polys, (a, y, x), strategy="auto")
    chordal = rank_variables_by_polynomials(polys, (a, y, x), strategy="chordal")
    assert score.strategy == "chordal"
    assert score.order == chordal


def test_auto_order_is_invariant_under_polynomial_permutation():
    a, b, c, d = sp.symbols("a b c d")
    polys = (a**2 + b * c, b**3 + c * d, a * d + c**2, b * d - 1)
    expected = suggest_cad_variable_order(polys, (a, b, c, d), strategy="auto")
    reordered = suggest_cad_variable_order(tuple(reversed(polys)), (a, b, c, d), strategy="auto")
    assert reordered.order == expected.order
    assert reordered.strategy == expected.strategy


def test_auto_order_commutes_with_symbol_renaming():
    a, b, c, d = sp.symbols("a b c d")
    w, x, y, z = sp.symbols("w x y z")
    polys = (a**2 + b * c, b**3 + c * d, a * d + c**2, b * d - 1)
    renamed = {a: w, b: x, c: y, d: z}
    original = suggest_cad_variable_order(polys, (a, b, c, d), strategy="auto")
    transformed = suggest_cad_variable_order(
        tuple(poly.xreplace(renamed) for poly in polys), (w, x, y, z), strategy="auto"
    )
    assert transformed.order == tuple(renamed[var] for var in original.order)
    assert transformed.strategy == original.strategy


def test_auto_order_is_deterministic_on_adversarial_ties():
    a, b, c, d = sp.symbols("a b c d")
    polys = (a * b + c * d, a * c + b * d, a * d + b * c)
    first = suggest_cad_variable_order(polys, (a, b, c, d), strategy="auto")
    for _ in range(5):
        assert suggest_cad_variable_order(polys, (a, b, c, d), strategy="auto") == first
