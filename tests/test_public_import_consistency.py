import pickle

import sympy as sp

import semialg
from semialg.optimization import function_range
from semialg.region_integrate import (
    ReducedRegionIntegral,
    RegionIntegralPiece,
    RegionIntegralResult,
)


def test_public_objects_report_their_defining_modules():
    assert function_range.__module__ == "semialg._optimization_range"
    assert RegionIntegralResult.__module__ == "semialg.region_integral_results"
    assert RegionIntegralPiece.__module__ == "semialg.region_integral_results"
    assert ReducedRegionIntegral.__module__ == "semialg.region_integral_results"


def test_region_integral_result_is_pickleable():
    x = sp.Symbol("x", real=True)
    result = RegionIntegralResult(sp.Integer(1), x, sp.true, (x,), "test")
    restored = pickle.loads(pickle.dumps(result))
    assert restored == result
    assert type(restored) is RegionIntegralResult


def test_every_declared_root_export_resolves_from_lazy_registry():
    for name in semialg.__all__:
        assert getattr(semialg, name) is not None
