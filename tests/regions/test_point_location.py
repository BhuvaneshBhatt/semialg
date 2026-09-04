from __future__ import annotations

import sympy as sp

from semialg import SemialgebraicRegion


def test_point_location_signs_agree_with_exact_substitution():
    x = sp.symbols("x", real=True)
    cad_region = SemialgebraicRegion(sp.And(x >= 0, x <= 2), (x,)).as_cad_region()
    location = cad_region.locate_point((sp.Rational(1, 2),))

    assert location.selected
    for entry in location.signs:
        value = sp.simplify(entry.polynomial.subs({x: sp.Rational(1, 2)}))
        expected = 0 if value == 0 else (1 if value > 0 else -1)
        assert entry.sign == expected


def test_section_point_has_zero_in_sign_vector():
    x = sp.symbols("x", real=True)
    cad_region = SemialgebraicRegion(sp.And(x >= 0, x <= 2), (x,)).as_cad_region()
    location = cad_region.locate_point((0,))
    assert location.selected
    assert 0 in location.sign_vector


def test_repeated_point_queries_reuse_identical_cad_object():
    x = sp.symbols("x", real=True)
    region = SemialgebraicRegion(sp.And(x >= 0, x <= 2), (x,))
    result = region.ensure_cad()
    cad_id = id(result.cad)

    wrapper = region.as_cad_region()
    wrapper.locate_point((sp.Rational(1, 3),))
    wrapper.locate_point((sp.sqrt(2),))

    assert region.ensure_cad() is result
    assert id(wrapper.result.cad) == cad_id
