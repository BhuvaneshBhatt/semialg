from __future__ import annotations

import pytest

from .known_regions import corpus


def _params():
    return tuple(
        pytest.param(case, id=case.name, marks=pytest.mark.slow if case.slow else ())
        for case in corpus()
    )


@pytest.mark.parametrize("case", _params())
def test_known_region_basic_invariants(case):
    region = case.region
    assert region.dimension() == case.dimension
    assert len(region.components()) == case.components
    assert region.is_bounded() is case.bounded
    assert region.is_closed() is case.closed
    if case.euler_compact_support is not None:
        assert region.euler_characteristic() == case.euler_compact_support


@pytest.mark.slow
def test_native_region_analysis_reused_across_invariants():
    case = next(item for item in corpus() if item.name == "solid_torus")
    region = case.region
    analysis = region.ensure_native_analysis()
    assert region.dimension() == 3
    assert len(region.components()) == 1
    assert region.is_bounded() is True
    assert region.is_closed() is True
    assert region.euler_characteristic() == 0
    assert region.ensure_native_analysis() is analysis
