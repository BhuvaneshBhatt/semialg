import sympy as sp

from semialg.heuristics import rank_variables_by_polynomials, suggest_cad_variable_order


def test_chordal_order_reverses_min_fill_elimination_order():
    a, b, c, d = sp.symbols("a b c d")
    polys = (a * b - 1, b * c - 1, c * d - 1)
    elimination = tuple(
        reversed(rank_variables_by_polynomials(polys, (a, b, c, d), strategy="chordal"))
    )
    suggested = suggest_cad_variable_order(polys, (a, b, c, d), strategy="chordal")
    assert suggested.order == tuple(reversed(elimination))
    assert set(suggested.order) == {a, b, c, d}
