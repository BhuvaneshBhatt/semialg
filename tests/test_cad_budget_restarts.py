"""Failed lifts and final verification remain inside a cooperative CAD budget."""

from unittest.mock import patch

import pytest
import sympy as sp

from semialg.cache_utils import cache_disabled
from semialg.cad_algorithms import decomp_collins_complete, decomp_lazard_complete
from semialg.cad_algorithms.decomposition import decomp_from_proj_tower
from semialg.errors import ProjectionNullification, ResourceLimitError


def test_restart_counts_failed_lift_work():
    x = sp.Symbol("x", real=True)

    def fail_first(tower, **kwargs):
        if kwargs["backend"] == "collins-complete":
            kwargs["budget"].check(3)
            raise ProjectionNullification("forced restart after three cells")
        return decomp_from_proj_tower(tower, **kwargs)

    with (
        cache_disabled(),
        patch(
            "semialg.cad_algorithms.decomposition.decomp_from_proj_tower", side_effect=fail_first
        ),
    ):
        with pytest.raises(ResourceLimitError, match="cell limit"):
            decomp_collins_complete((x,), (x,), max_cells=5)
        result = decomp_collins_complete((x,), (x,), max_cells=6)
    assert result.backend == "lazard-complete"
    assert len(result.cells) == 3


def test_lazard_checks_deadline_after_verification():
    x = sp.Symbol("x", real=True)
    source = decomp_collins_complete((x,), (x,))
    with (
        cache_disabled(),
        patch("semialg.cad_algorithms.decomposition.decomp_from_proj_tower", return_value=source),
        patch("semialg.cad_algorithms.resource_limits.perf_counter", side_effect=[0, 0, 2]),
    ):
        with pytest.raises(ResourceLimitError, match="deadline"):
            decomp_lazard_complete((x,), (x,), timeout=1)


def test_direct_lift_inherits_budget():
    from semialg.cad_algorithms.projection.collins import build_collins_proj_set
    from semialg.cad_algorithms.resource_limits import CADBudget, budget_scope

    x = sp.Symbol("x", real=True)
    tower = build_collins_proj_set((x,), (x,))
    with budget_scope(CADBudget(2, None)):
        with pytest.raises(ResourceLimitError):
            decomp_from_proj_tower(tower)
    with budget_scope(CADBudget(3, None)):
        result = decomp_from_proj_tower(tower)
    assert len(result.cells) == 3


@pytest.mark.parametrize("strategy", ["mccallum", "lazard", "tticad"])
def test_reduced_and_reference_work(strategy):
    from semialg import quantifier_eliminate
    from semialg.cad_algorithms.decomposition import _lift_projection_tower

    x, a = sp.symbols("x a", real=True)
    constructed = []

    def observe(*args, **kwargs):
        result = _lift_projection_tower(*args, **kwargs)
        constructed.append(sum(map(len, result.cells_by_level.values())))
        return result

    with (
        cache_disabled(),
        patch("semialg.cad_algorithms.decomposition._lift_projection_tower", side_effect=observe),
    ):
        reference = quantifier_eliminate(x**3 + a < 0, (("exists", x),), strategy=strategy)
    assert reference is sp.true
    assert len(constructed) >= 2
    total = sum(constructed)
    with cache_disabled():
        with pytest.raises(ResourceLimitError):
            quantifier_eliminate(
                x**3 + a < 0, (("exists", x),), strategy=strategy, max_cells=total - 1
            )
        assert (
            quantifier_eliminate(x**3 + a < 0, (("exists", x),), strategy=strategy, max_cells=total)
            is sp.true
        )


def test_failed_direct_lift_is_charged():
    from semialg.cad_algorithms.projection.collins import build_collins_proj_set
    from semialg.cad_algorithms.resource_limits import CADBudget, budget_scope

    x = sp.Symbol("x", real=True)
    tower = build_collins_proj_set((x,), (x,))

    def fail(*args, **kwargs):
        kwargs["budget"].check(3)
        raise ProjectionNullification("forced direct lift failure")

    with budget_scope(CADBudget(5, None)):
        with patch("semialg.cad_algorithms.decomposition._lift_projection_tower", side_effect=fail):
            with pytest.raises(ProjectionNullification):
                decomp_from_proj_tower(tower)
        with pytest.raises(ResourceLimitError):
            decomp_collins_complete((x,), (x,))


def test_instance_searches_share_work():
    from semialg.cad_algorithms.resource_limits import CADBudget, budget_scope
    from semialg.formula import parse_formula
    from semialg.partial.qe import lazy_find_inst_form

    x = sp.Symbol("x", real=True)
    with budget_scope(CADBudget(5, None)):
        first = lazy_find_inst_form((x,), parse_formula(x > 0))
        assert first.instance[x] == 1
        with pytest.raises(ResourceLimitError):
            lazy_find_inst_form((x,), parse_formula(x > 0))
