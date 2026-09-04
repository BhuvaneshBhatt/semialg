from __future__ import annotations

import pytest
import sympy as sp

from semialg._linear_relations import certified_sign, safe_linear_solution
from semialg.cad_algorithms import clear_cad_caches
from semialg.cad_algorithms.decomposition import decomp_collins_complete
from semialg.cad_algorithms.projection.collins import ProjectionTower, build_collins_proj_set
from semialg.decision import _metadata
from semialg.instances import coordinate_bounds, is_bounded_solution_set
from semialg.simplify.equality import simplify_equalities


def test_simplify_equalities_preserves_parameter_zero_stratum():
    a, x, y = sp.symbols("a x y", real=True)
    formula = sp.And(sp.Eq(a * x, 0), y > 0)

    simplified = simplify_equalities(formula)

    assert bool(formula.subs({a: 2, x: 0, y: 1}))
    assert bool(simplified.subs({a: 2, x: 0, y: 1}))
    assert sp.Eq(a * x, 0) in sp.And.make_args(simplified)


def test_simplify_equalities_may_divide_by_assumption_certified_coefficient():
    a = sp.Symbol("a", positive=True)
    x, y = sp.symbols("x y", real=True)

    simplified = simplify_equalities(sp.And(sp.Eq(a * x, 0), y > 0))

    assert sp.Eq(x, 0) in sp.And.make_args(simplified)


def test_safe_linear_solution_rejects_unresolved_parameter_coefficient():
    a, x = sp.symbols("a x", real=True)
    assert safe_linear_solution(a * x, x) is None
    assert certified_sign(a) is None


def test_coordinate_bounds_reject_unresolved_parameter_sign():
    a, x = sp.symbols("a x", real=True)
    bounds = coordinate_bounds(a * x <= 1, (x,))

    assert bounds.bounds == ((x, -sp.oo, sp.oo),)
    assert bounds.complete is False


def test_coordinate_bounds_accept_assumption_certified_parameter_sign():
    a = sp.Symbol("a", positive=True)
    x = sp.Symbol("x", real=True)
    bounds = coordinate_bounds(a * x <= 1, (x,))

    assert bounds.bounds == ((x, -sp.oo, 1 / a),)
    assert bounds.complete is True


def test_boundedness_does_not_drop_parameter_zero_stratum():
    a, x = sp.symbols("a x", real=True)
    formula = sp.And(a * x >= -1, a * x <= 1)

    assert is_bounded_solution_set(formula, (x,)) is None
    assert bool(formula.subs({a: 0, x: 10**6}))


def test_complete_cad_cache_result_is_transitively_read_only():
    x = sp.Symbol("x", real=True)
    clear_cad_caches()
    cad = decomp_collins_complete((x**2 - 1,), (x,))

    with pytest.raises(TypeError):
        cad.cells_by_level[1] = ()  # type: ignore[index]
    with pytest.raises(TypeError):
        cad.cells[0].signs["corrupt"] = 1  # type: ignore[index]

    again = decomp_collins_complete((x**2 - 1,), (x,))
    assert again is cad
    assert again.cells


def test_projection_tower_metadata_is_recursively_read_only():
    x = sp.Symbol("x", real=True)
    tower = ProjectionTower(
        variables=(x,),
        levels=(),
        original_polynomials=(),
        metadata={"nested": {"value": 1}},
    )

    with pytest.raises(TypeError):
        tower.metadata["new"] = 2  # type: ignore[index]
    nested = tower.metadata["nested"]
    assert isinstance(nested, dict) is False
    with pytest.raises(TypeError):
        nested["value"] = 2  # type: ignore[index]


def test_cached_projection_metadata_cannot_be_poisoned():
    x = sp.Symbol("x", real=True)
    clear_cad_caches()
    tower = build_collins_proj_set((x**2 - 1,), (x,))
    with pytest.raises(TypeError):
        tower.metadata["complete"] = False  # type: ignore[index]
    assert build_collins_proj_set((x**2 - 1,), (x,)) is tower
    assert tower.metadata["complete"] is True


def test_affine_dimension_uses_rank_not_equality_count():
    x, y = sp.symbols("x y", real=True)
    formula = sp.And(sp.Eq(x, 0), sp.Eq(2 * x, 0))

    assert _metadata._affine_equality_dimension(formula, (x, y)) == 1


def test_dimension_fallback_is_unknown_when_inequalities_can_change_dimension():
    x = sp.Symbol("x", real=True)
    formula = sp.And(x >= 0, x <= 0)

    metadata = {"dimension": None}
    _metadata._infer_dimension(metadata, formula, (x,))
    assert metadata["dimension"] is None


def test_affine_dimension_rejects_parameter_dependent_strata():
    a, x = sp.symbols("a x", real=True)
    assert _metadata._affine_equality_dimension(sp.Eq(a * x, 0), (x,)) is None
