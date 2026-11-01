"""Free-cell coverage and alternating-prefix differential checks."""

import pickle
from itertools import product
from unittest.mock import patch

import pytest
import sympy as sp

from semialg import quantifier_eliminate, replay_certificate
from semialg.cache_utils import cache_disabled
from semialg.cad_algorithms import decomp_collins_complete
from semialg.errors import ProjectionNullification, ResourceLimitError
from semialg.formula import parse_formula
from semialg.partial.qe import lazy_quantifier_eliminate
from semialg.qe.complete import qe_from_cad


@pytest.mark.parametrize("kinds", list(product(("exists", "forall"), repeat=2)))
@pytest.mark.parametrize("form", ["equality", "open", "birth"])
def test_free_prefix_differential(kinds, form):
    a, x, y = sp.symbols("a x y", real=True)
    matrix = {
        "equality": sp.Eq(y * y, x * x + a),
        "open": sp.Or(y * y < a, x * x > a),
        "birth": sp.And(a >= 0, sp.Eq(y * y, a * x * x)),
    }[form]
    prefix = tuple(zip(kinds, (x, y), strict=True))
    partial = lazy_quantifier_eliminate((a, x, y), prefix, matrix, return_result=True)
    full = decomp_collins_complete(
        (
            matrix_poly.lhs - matrix_poly.rhs
            for matrix_poly in matrix.atoms(sp.core.relational.Relational)
        ),
        (a, x, y),
    )
    reference = qe_from_cad(full, (a, x, y), prefix, parse_formula(matrix), return_result=True)
    assert partial.truth_by_index == reference.cell_truth
    assert set(partial.truth_by_index) == {cell.index for cell in partial.cells_by_level[1]}
    assert sum(partial.stats.visited_cells_by_level.values()) <= sum(
        map(len, full.cells_by_level.values())
    )
    for value in (-2, -1, 0, 1, 2):
        assert bool(partial.formula.subs(a, value)) == bool(reference.formula.subs(a, value))


def test_two_free_coordinates():
    a, b, x = sp.symbols("a b x", real=True)
    result = lazy_quantifier_eliminate(
        (a, b, x), (("exists", x),), sp.Eq(x * x, a * b), return_result=True
    )
    assert len(result.truth_by_index) == len(result.cells_by_level[2])
    for av, bv in product((-1, 0, 1), repeat=2):
        assert bool(result.formula.subs({a: av, b: bv})) == (av * bv >= 0)


def test_free_dispatch_and_serialization():
    a, x = sp.symbols("a x", real=True)
    result = quantifier_eliminate(
        sp.Eq(x * x, a), (("exists", x),), strategy="partial", return_result=True
    )
    restored = pickle.loads(pickle.dumps(result.backend_result))
    assert restored.truth_by_index == {(0,): False, (1,): True, (2,): True}
    assert replay_certificate(restored).verified is True
    with pytest.raises(TypeError):
        restored.truth_by_index[(0,)] = True


def test_free_budget_and_restart():
    a, x = sp.symbols("a x", real=True)
    with pytest.raises(ResourceLimitError):
        lazy_quantifier_eliminate((a, x), (("exists", x),), sp.Eq(x * x, a), max_cells=2)
    from semialg.partial.qe import _stack_for_parent

    def fail_once(state, parent, level, **kwargs):
        if state.tower.metadata.get("projection") != "lazard-complete":
            raise ProjectionNullification("forced ordinary failure")
        return _stack_for_parent(state, parent, level, **kwargs)

    with cache_disabled(), patch("semialg.partial.qe._stack_for_parent", side_effect=fail_once):
        result = lazy_quantifier_eliminate(
            (a, x), (("exists", x),), sp.Eq(x * x, a), return_result=True
        )
    assert result.backend == "partial-lazard-free"
    assert result.truth_by_index == {(0,): False, (1,): True, (2,): True}


@pytest.mark.parametrize(
    "kinds", [kinds for size in (1, 2, 3) for kinds in product(("exists", "forall"), repeat=size)]
)
def test_short_free_prefixes(kinds):
    a = sp.Symbol("a", real=True)
    coordinates = sp.symbols("x:z", real=True)[: len(kinds)]
    matrix = sp.Or(*(sp.Eq(variable**2, a) for variable in coordinates))
    result = lazy_quantifier_eliminate(
        (a, *coordinates), tuple(zip(kinds, coordinates, strict=True)), matrix, return_result=True
    )
    for value in (-1, 0, 1):
        assert bool(result.formula.subs(a, value)) == (value >= 0 and "exists" in kinds)
    assert len(result.truth_by_index) == len(result.cells_by_level[1])


def test_temporary_cells_are_budgeted():
    a, x = sp.symbols("a x", real=True)
    matrix = sp.Eq(x * x, a)
    with cache_disabled():
        with pytest.raises(ResourceLimitError):
            lazy_quantifier_eliminate((a, x), (("exists", x),), matrix, max_cells=11)
        result = lazy_quantifier_eliminate(
            (a, x), (("exists", x),), matrix, max_cells=12, return_result=True
        )
    assert sum(result.stats.visited_cells_by_level.values()) == 6
    assert result.stats.ec_skipped_sector_cells == 6
    assert result.formula == (a >= 0)


def test_forced_complete_backend_failure():
    a, x = sp.symbols("a x", real=True)
    with patch(
        "semialg.partial.qe._stack_for_parent",
        side_effect=ProjectionNullification("forced both backends"),
    ):
        with pytest.raises(ProjectionNullification, match="forced both backends"):
            quantifier_eliminate(sp.Eq(x * x, a), (("exists", x),), strategy="partial")


def test_partial_coordinate_roles():
    a, x, y = sp.symbols("a x y", real=True)
    prefix = (("forall", y), ("exists", x))
    result = lazy_quantifier_eliminate(
        (x, a, y), prefix, sp.Eq(x * x, y * y + a), return_result=True
    )
    assert result.variables == (a, y, x)
    assert result.free_variables == (a,)
    assert result.quantifiers == prefix
    assert result.formula == (a >= 0)


def test_positive_dimensional_free_nullification():
    x, y, t, z, w = sp.symbols("x y t z w", real=True)
    result = lazy_quantifier_eliminate(
        (x, y, t, z, w),
        (("forall", z), ("exists", w)),
        sp.Eq(w * w + y * z - x * x, 0),
        return_result=True,
    )
    assert result.backend == "partial-lazard-free"
    assert result.stats.ec_section_lifts == 0
    for xv, yv, tv in product((-1, 0, 1), repeat=3):
        assert bool(result.formula.subs({x: xv, y: yv, t: tv})) == (yv == 0)
