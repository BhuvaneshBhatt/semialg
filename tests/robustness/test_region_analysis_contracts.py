import pytest
import sympy as sp

from semialg.region_analysis import _point_subs, region_singular_locus
from semialg.structural_keys import ordered_symbols, symbol_identity_key


def test_ordered_symbols_preserves_full_same_name_symbol_identity():
    xr = sp.Symbol("x", real=True)
    xp = sp.Symbol("x", positive=True)

    assert ordered_symbols({xp, xr}) == tuple(sorted((xp, xr), key=symbol_identity_key))


def test_region_singular_locus_respects_explicit_empty_ambient_variables():
    x = sp.Symbol("x", real=True)

    # Here x is a parameter because the caller explicitly requested a
    # zero-dimensional ambient space.  It must not be silently promoted to an
    # ambient coordinate merely because it occurs free in the formula.
    assert region_singular_locus(x**2 <= 0, variables=()) is sp.false


def test_region_analysis_rejects_missing_coordinates():
    x, y = sp.symbols("x y", real=True)

    with pytest.raises(ValueError, match="missing coordinate.*y"):
        _point_subs({x: 1}, (x, y))


def test_region_nonsmooth_locus_detects_square_corners_but_not_smooth_circle():
    from semialg.region_analysis import region_nonsmooth_locus

    x, y = sp.symbols("x y", real=True)
    square = sp.And(x >= 0, x <= 1, y >= 0, y <= 1)
    corners = region_nonsmooth_locus(square, (x, y))
    for point in ({x: 0, y: 0}, {x: 0, y: 1}, {x: 1, y: 0}, {x: 1, y: 1}):
        assert sp.simplify(corners.subs(point)) is sp.true
    assert sp.simplify(corners.subs({x: sp.Rational(1, 2), y: 0})) is sp.false
    assert region_nonsmooth_locus(x**2 + y**2 <= 1, (x, y)) is sp.false


def test_region_nonsmooth_locus_does_not_call_independent_equalities_corners():
    from semialg.region_analysis import region_nonsmooth_locus

    x, y = sp.symbols("x y", real=True)
    assert region_nonsmooth_locus(sp.And(sp.Eq(x, 0), sp.Eq(y, 0)), (x, y)) is sp.false


def test_active_boundary_strata_separate_square_edges_and_corners():
    from semialg.region_analysis import region_active_boundary_strata

    x, y = sp.symbols("x y", real=True)
    square = sp.And(x >= 0, x <= 1, y >= 0, y <= 1)
    strata = region_active_boundary_strata(square, (x, y))

    assert sum(stratum.active_count == 1 for stratum in strata) == 4
    assert sum(stratum.active_count == 2 for stratum in strata) == 4
    corner = {x: 0, y: 0}
    matching = [
        stratum for stratum in strata if sp.simplify(stratum.formula.subs(corner)) is sp.true
    ]
    assert len(matching) == 1
    assert matching[0].active_count == 2
