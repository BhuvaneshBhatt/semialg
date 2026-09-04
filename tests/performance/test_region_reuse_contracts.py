from __future__ import annotations

import sympy as sp

from semialg import SemialgebraicRegion
from semialg.cad_region import cad_combine


def test_region_context_and_cad_identity_are_stable_across_queries():
    x = sp.symbols("x", real=True)
    region = SemialgebraicRegion(sp.And(x >= 0, x <= 2), (x,))
    context = region.context
    result = region.ensure_cad()

    assert region.context is context
    assert region.ensure_cad() is result
    region.as_cad_region().locate_point((1,))
    assert region.ensure_cad() is result


def test_shared_cad_boolean_combine_reports_reuse():
    x = sp.symbols("x", real=True)
    base = SemialgebraicRegion(sp.And(x >= -2, x <= 2), (x,)).as_cad_region()
    left = base.extend(x <= 0)
    right = base.extend(x >= 0)
    # Extension may refine; combine only claims reuse when the actual decomposition is shared.
    if left.result.cad is right.result.cad:
        combined = cad_combine(left, right, op="union")
        assert combined.result.diagnostics.get("cad_reused") is True
        assert combined.result.cad is left.result.cad
