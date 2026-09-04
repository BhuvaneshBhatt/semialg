from __future__ import annotations

from functools import lru_cache

import pytest
import sympy as sp

hypothesis = pytest.importorskip("hypothesis")
from hypothesis import given, settings
from hypothesis import strategies as st

from semialg import connected_components, root_count_conditions
from semialg.algebraic.rational_univariate import solve_zero_dimensional_system_with_rur
from semialg.cad_algorithms.reduced import decomp_form_reduced_safe
from semialg.formula import parse_formula

small_int = st.integers(min_value=-3, max_value=3)


@lru_cache(maxsize=1)
def _cubic_classifier():
    x, aa, bb, cc, dd = sp.symbols("x aa bb cc dd", real=True)
    family = aa * x**3 + bb * x**2 + cc * x + dd
    conditions = root_count_conditions(family, x, (aa, bb, cc, dd))
    return x, (aa, bb, cc, dd), family, conditions


def _selected_count(conditions, assignment):
    selected = [
        count
        for count, condition in conditions.items()
        if sp.simplify(condition.subs(assignment)) in (sp.true, True)
    ]
    assert len(selected) == 1
    return selected[0]


@settings(max_examples=30, deadline=None)
@given(a=small_int, b=small_int, c=small_int, d=small_int)
def test_generated_cubic_strata_match_exact_specialization(a, b, c, d):
    x, params, family, conditions = _cubic_classifier()
    aa, bb, cc, dd = params
    assignment = {aa: a, bb: b, cc: c, dd: d}
    specialized = sp.expand(family.subs(assignment))
    selected = _selected_count(conditions, assignment)

    if specialized == 0:
        assert selected is sp.oo
    else:
        roots = tuple(dict.fromkeys(sp.real_roots(sp.Poly(specialized, x))))
        assert selected == len(roots)


@settings(max_examples=20, deadline=None)
@given(
    px=st.integers(min_value=-2, max_value=2),
    constant=st.integers(min_value=-2, max_value=2),
)
def test_generated_reduced_cad_attempt_is_proof_carrying(px, constant):
    x, y = sp.symbols("x y", real=True)
    ec = y**2 + px * x * y + x + constant
    formula = parse_formula(sp.And(sp.Eq(ec, 0), x >= -2, x <= 2))

    for backend in ("mccallum", "lazard"):
        result = decomp_form_reduced_safe(formula, (x, y), backend=backend)
        assert result.complete
        if result.used_fallback:
            assert result.fallback_cad is not None
        else:
            assert result.certificate is not None and result.certificate.valid
            assert result.side_conditions is not None and result.side_conditions.valid


@settings(max_examples=15, deadline=None)
@given(index=st.integers(min_value=0, max_value=4))
def test_generated_rur_algebraic_field_transition_preserves_solution(index):
    x, y = sp.symbols("x y", real=True)
    scalars = (
        sp.Integer(1),
        sp.sqrt(2),
        sp.sqrt(3),
        1 + sp.sqrt(2),
        sp.sqrt(2) + sp.sqrt(3),
    )
    scalar = scalars[index]
    equations = (
        sp.expand(scalar * (x**2 - 2)),
        sp.expand(scalar * (y - x)),
    )

    solutions = solve_zero_dimensional_system_with_rur(equations, (x, y), real=True)

    assert solutions == ((-sp.sqrt(2), -sp.sqrt(2)), (sp.sqrt(2), sp.sqrt(2)))


@settings(max_examples=12, deadline=None)
@given(parameter=st.sampled_from((-2, -1, 0, 1, 2)))
def test_generated_singular_fiber_component_count(parameter):
    x, y = sp.symbols("x y", real=True)
    fiber = sp.Eq(y**2, x**2 * (x + parameter))
    components = connected_components(fiber, (x, y))

    expected = 2 if parameter < 0 else 1
    assert len(components) == expected
