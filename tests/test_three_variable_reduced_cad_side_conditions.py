from __future__ import annotations

import pytest
import sympy as sp

from semialg.cad_algorithms.decomposition import decomp_from_proj_tower
from semialg.cad_algorithms.reduced import (
    _scan_reduced_conditions,
    build_reduced_proj,
    decompose_reduced_safe,
    delin_repair_polys,
)


def _reduced_attempt(polys, variables, ec, backend):
    projection = build_reduced_proj(
        polys,
        variables,
        backend=backend,
        equational_constraints=(ec,),
        certify_reduced=True,
    )
    cad = decomp_from_proj_tower(projection.tower)
    return projection, cad, _scan_reduced_conditions(cad, backend=backend)


def test_mccallum_detects_nested_three_variable_nullification():
    x, y, z = sp.symbols("x y z", real=True)
    ec = x * z + y

    _, _, side = _reduced_attempt((ec, z - 1), (x, y, z), ec, "mccallum")

    assert side.valid is False
    levels = {event.level for event in side.nullification_events}
    assert {2, 3} <= levels
    assert any(event.polynomial == x * z + y for event in side.nullification_events)


def test_three_variable_nullification_produces_lower_level_repair_factors():
    x, y, z = sp.symbols("x y z", real=True)
    ec = x * z + y
    _, _, side = _reduced_attempt((ec, z - 1), (x, y, z), ec, "mccallum")

    repairs = {
        event.level: tuple(poly.as_expr() for poly in delin_repair_polys(event, (x, y, z)))
        for event in side.nullification_events
    }

    assert x in repairs[2]
    assert x in repairs[3]
    assert y in repairs[3]


def test_lazard_valuation_avoids_false_mccallum_nullification_for_nested_case():
    x, y, z = sp.symbols("x y z", real=True)
    ec = x * z + y

    _, _, mccallum = _reduced_attempt((ec, z - 1), (x, y, z), ec, "mccallum")
    _, _, lazard = _reduced_attempt((ec, z - 1), (x, y, z), ec, "lazard")

    assert mccallum.nullification_events
    assert lazard.valid
    assert not lazard.nullification_events


def test_full_safe_mccallum_never_accepts_uncertified_nested_nullification():
    x, y, z = sp.symbols("x y z", real=True)
    ec = x * z + y

    result = decompose_reduced_safe(
        (ec, z - 1),
        (x, y, z),
        backend="mccallum",
        equational_constraints=(ec,),
    )

    assert result.complete
    if result.used_fallback:
        assert result.fallback_cad is not None
    else:
        assert result.certificate is not None and result.certificate.valid
        assert result.side_conditions is not None and result.side_conditions.valid


@pytest.mark.parametrize(
    "ec",
    [
        lambda x, y, z: x * z + y,
        lambda x, y, z: (x**2 + y**2) * z + x * y,
        lambda x, y, z: x * y * z + x + y,
        lambda x, y, z: x * z**2 + y * z + x,
    ],
)
def test_generated_three_variable_side_reports_are_self_consistent(ec):
    x, y, z = sp.symbols("x y z", real=True)
    equation = ec(x, y, z)

    for backend in ("mccallum", "lazard"):
        _, _, side = _reduced_attempt((equation, z - 1), (x, y, z), equation, backend)
        assert side.valid is (not side.nullification_events)
        if side.nullification_events:
            assert side.failed_conditions
            for event in side.nullification_events:
                assert event.level in (2, 3)
                assert event.variable in (y, z)
