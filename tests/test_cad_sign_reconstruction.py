"""Existing ambient CAD signs can certify a compact projected truth formula."""

from unittest.mock import patch

import sympy as sp

from semialg.cad_algorithms import decomp_collins_complete
from semialg.cad_algorithms.reduced import decomp_form_reduced_safe
from semialg.formula import parse_formula
from semialg.qe.complete import qe_from_cad
from semialg.simplify.cell_signs import formula_from_cad_signs


def test_quadratic_projection_reuses_signs():
    x, a, b, y = sp.symbols("x a b y", real=True)
    variables = (x, a, b, y)
    matrix = parse_formula(sp.And(sp.Eq(x * y * y + a * y + b, 0), b < 5))
    source = decomp_form_reduced_safe(matrix, variables, backend="lazard").cad
    with patch(
        "semialg.simplify.implication._is_unsatisfiable_cached",
        side_effect=AssertionError("new CAD implication"),
    ):
        result = qe_from_cad(source, variables, (("exists", y),), matrix)
    for xv in (-1, 0, 1):
        for av in (-2, 0, 2):
            for bv in (-1, 0, 1, 5):
                expected = bv < 5 and (av * av - 4 * xv * bv >= 0 if xv else av != 0 or bv == 0)
                assert bool(result.subs({x: xv, a: av, b: bv})) is expected


def test_equal_signs_with_different_truth_decline():
    x, y = sp.symbols("x y", real=True)
    source = decomp_collins_complete((x * x - 2, y), (x, y))
    truth = {cell.index: cell.sample[0].interval.left > 0 for cell in source.cells}
    assert formula_from_cad_signs(source, (x, y), truth) is None


def test_truth_only_source_declines():
    x, y, z = sp.symbols("x y z", real=True)
    matrix = parse_formula(sp.And(sp.Eq(y - x, 0), sp.Eq(z * z - y, 0)))
    source = decomp_form_reduced_safe(matrix, (x, y, z), backend="mccallum").cad
    assert formula_from_cad_signs(source, (x, y), {}) is None
