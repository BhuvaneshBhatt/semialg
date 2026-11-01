"""Mathematical and operation-count contracts for exact-work reuse."""

from unittest.mock import patch

import pytest
import sympy as sp

from semialg.algebraic.equality_ideal import EqualityIdealContext, _equality_groebner_basis
from semialg.cache_control import clear_caches
from semialg.cad_algorithms.decomposition import decomp_collins_complete
from semialg.formula import parse_formula
from semialg.instances.real_utils import satisfies_formula
from semialg.qe.complete import evaluate_formula_on_cell
from semialg.sampling import sample_points


def test_sampling_validates_each_candidate_once():
    x, y = sp.symbols("x y", real=True)
    formula = sp.And(sp.Eq(x, 1), sp.Eq(y, 1))
    with patch("semialg.sampling.satisfies_formula", wraps=satisfies_formula) as validate:
        result = sample_points(
            formula, (x, y), strategy="grid", bounds=((0, 1), (0, 1)), grid_resolution=7
        )
    assert result == ({x: 1, y: 1},)
    assert validate.call_count == 49


def test_equality_bases_preserve_fields_and_teardown():
    x, y = sp.symbols("x y", real=True)
    clear_caches(collect=False)
    rational = EqualityIdealContext.build((x**2 - 2, y - x), (x, y))
    repeated = EqualityIdealContext.build((x**2 - 2, y - x), (x, y))
    assert rational.groebner_basis is repeated.groebner_basis
    algebraic = EqualityIdealContext.build((x - sp.sqrt(2), y - x), (x, y))
    assert rational.quotient_dimension == 2
    assert algebraic.quotient_dimension == 1
    assert algebraic.normal_form(x - sp.sqrt(2)) == 0
    assert rational.normal_form(y) == y
    assert algebraic.normal_form(x) == sp.sqrt(2)
    clear_caches(collect=False)
    assert _equality_groebner_basis.cache_info().currsize == 0


@pytest.mark.parametrize("op", ["eq", "ne", "lt", "le", "gt", "ge"])
def test_qe_uses_recorded_exact_signs(op):
    x = sp.Symbol("x", real=True)
    source = decomp_collins_complete((x,), (x,))
    relation = {
        "eq": sp.Eq(x, 0),
        "ne": sp.Ne(x, 0),
        "lt": x < 0,
        "le": x <= 0,
        "gt": x > 0,
        "ge": x >= 0,
    }[op]
    with patch(
        "semialg.qe.complete.sign_at_sample",
        side_effect=AssertionError("redundant sign evaluation"),
    ):
        for cell in source.cells:
            assert evaluate_formula_on_cell(parse_formula(relation), cell, (x,)) is bool(
                relation.subs(x, cell.sample_exprs[0])
            )
    # Scaled and repeated factors need their own signs; a normalized factor's
    # recorded sign cannot be reused without its scalar and multiplicity.
    for expression in (-7 * x, x**2, -(x**2)):
        matrix = parse_formula(expression >= 0)
        for cell in source.cells:
            assert evaluate_formula_on_cell(matrix, cell, (x,)) is bool(
                (expression >= 0).subs(x, cell.sample_exprs[0])
            )


def test_point_location_reuses_shared_prefixes():
    x, y = sp.symbols("x y", real=True)
    source = decomp_collins_complete((x, y - x), (x, y))
    locator = source._point_locator
    points = ((0, 0), (0, 1), (0, 2), (1, 2), (0, -1))
    with patch.object(locator, "_locate_in_stack", wraps=locator._locate_in_stack) as locate:
        found = [locator.locate(point) for point in points]
    assert [cell.index for cell in found] == [(1, 1), (1, 2), (1, 2), (2, 2), (1, 0)]
    assert locate.call_count == 8
