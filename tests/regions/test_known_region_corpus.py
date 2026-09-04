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
