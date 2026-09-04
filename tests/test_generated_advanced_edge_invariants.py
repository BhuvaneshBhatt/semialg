from __future__ import annotations

import pytest
import sympy as sp

hypothesis = pytest.importorskip("hypothesis")
from hypothesis import given, settings
from hypothesis import strategies as st

from semialg import connected_components, root_count_conditions
from semialg.algebraic.rational_univariate import compute_rational_univariate_representation
from semialg.cad_algorithms.decomposition import decomp_from_proj_tower
from semialg.cad_algorithms.reduced import _scan_reduced_conditions, build_reduced_proj

x = sp.Symbol("x", real=True)
a, b, c, d, e = sp.symbols("a b c d e", real=True)
quartic = a * x**4 + b * x**3 + c * x**2 + d * x + e
quartic_conditions = root_count_conditions(quartic, x, (a, b, c, d, e))


def _selected_count(conditions, assignment):
    selected = [
        count
        for count, condition in conditions.items()
        if sp.simplify(condition.subs(assignment)) in (sp.true, True)
    ]
    assert len(selected) == 1
    return selected[0]


@settings(max_examples=40, deadline=None)
@given(
    aa=st.integers(-3, 3),
    bb=st.integers(-3, 3),
    cc=st.integers(-3, 3),
    dd=st.integers(-3, 3),
    ee=st.integers(-3, 3),
)
def test_generated_quartic_certified_counts_match_exact_roots(aa, bb, cc, dd, ee):
    assignment = {a: aa, b: bb, c: cc, d: dd, e: ee}
    specialized = sp.expand(quartic.subs(assignment))
    selected = _selected_count(quartic_conditions, assignment)

    if selected == -1:
        return
    if specialized == 0:
        assert selected is sp.oo
        return
    roots = tuple(dict.fromkeys(sp.real_roots(sp.Poly(specialized, x))))
    assert selected == len(roots)


@settings(max_examples=12, deadline=None)
@given(
    px=st.integers(-2, 2),
    py=st.integers(-2, 2),
    qx=st.integers(-2, 2),
    qy=st.integers(-2, 2),
)
def test_generated_three_variable_mccallum_side_scan_is_consistent(px, py, qx, qy):
    x0, y0, z0 = sp.symbols("x0 y0 z0", real=True)
    ec = (px * x0 + py * y0) * z0 + qx * x0 + qy * y0
    if ec == 0:
        return
    projection = build_reduced_proj(
        (ec, z0 - 1),
        (x0, y0, z0),
        backend="mccallum",
        equational_constraints=(ec,),
        certify_reduced=True,
    )
    cad = decomp_from_proj_tower(projection.tower)
    side = _scan_reduced_conditions(cad, backend="mccallum")

    assert side.valid is (not side.nullification_events)
    if side.nullification_events:
        assert side.failed_conditions


@settings(max_examples=12, deadline=None)
@given(
    multiplicity=st.integers(2, 4),
    field_index=st.integers(0, 2),
)
def test_generated_nonreduced_rur_dimension_is_stable_over_algebraic_fields(
    multiplicity, field_index
):
    x0, y0 = sp.symbols("x0 y0", real=True)
    alpha = (sp.Integer(1), sp.sqrt(2), sp.sqrt(3))[field_index]
    system = ((x0 - alpha) ** multiplicity, y0 - x0)

    representation = compute_rational_univariate_representation(system, (x0, y0))

    assert representation.quotient_dimension == multiplicity
    assert representation.geometric_solution_count == 1


@settings(max_examples=8, deadline=None)
@given(parameter=st.integers(-3, 3))
def test_generated_irreducible_singular_fiber_topology(parameter):
    x0, y0 = sp.symbols("x0 y0", real=True)
    poly = y0**2 - x0**2 * (x0 - parameter)
    factors = sp.factor_list(poly, x0, y0)[1]
    assert len(factors) == 1 and factors[0][1] == 1

    components = connected_components(sp.Eq(poly, 0), (x0, y0))

    expected = 2 if parameter > 0 else 1
    assert len(components) == expected
