from __future__ import annotations

from collections import OrderedDict

import sympy as sp
from hypothesis import given, settings
from hypothesis import strategies as st

from semialg._linear_relations import linear_relation_bound, safe_linear_solution
from semialg.instances import coordinate_bounds, is_bounded_solution_set
from semialg.instances.real_fallbacks import is_valid_numeric_value
from semialg.qe.complete import _canonicalize_same_kind_quantifier_blocks
from semialg.structural_keys import point_key, symbol_identity_key


@given(
    coeff=st.integers(min_value=-7, max_value=7).filter(lambda value: value != 0),
    offset=st.integers(min_value=-7, max_value=7),
    sample=st.integers(min_value=-10, max_value=10),
)
def test_certified_affine_bound_agrees_with_exact_specialization(coeff, offset, sample):
    x = sp.Symbol("x", real=True)
    relation = sp.Le(coeff * x + offset, 0)
    bound = linear_relation_bound(relation, x)

    assert bound is not None
    lhs_truth = bool(relation.subs(x, sample))
    if bound.side == "upper":
        rhs_truth = bool(sp.Le(sample, bound.value))
    else:
        rhs_truth = bool(sp.Ge(sample, bound.value))
    assert lhs_truth == rhs_truth


@given(
    coeff=st.integers(min_value=-7, max_value=7).filter(lambda value: value != 0),
    offset=st.integers(min_value=-7, max_value=7),
)
def test_safe_affine_equality_solution_replays_exactly(coeff, offset):
    x = sp.Symbol("x", real=True)
    solution = safe_linear_solution(coeff * x + offset, x)
    assert solution is not None
    assert sp.simplify((coeff * x + offset).subs(x, solution)) == 0


def test_parameter_dependent_affine_fast_paths_refuse_uncertified_strata():
    a, x = sp.symbols("a x", real=True)
    assert safe_linear_solution(a * x - 1, x) is None
    assert linear_relation_bound(a * x <= 1, x) is None
    assert coordinate_bounds(a * x <= 1, (x,)).complete is False


def test_coordinate_bounds_preserve_strictness_and_empty_semantics():
    x = sp.Symbol("x", real=True)
    contradictory = sp.And(x > 0, x < 0, evaluate=False)
    bounds = coordinate_bounds(contradictory, (x,))
    assert bounds.bounds == ((x, 0, 0),)
    assert bounds.strictness == ((x, True, True),)
    assert bounds.inconsistent is True
    assert is_bounded_solution_set(contradictory, (x,)) is True

    empty = coordinate_bounds(sp.false, (x,))
    assert empty.inconsistent is True
    assert empty.complete is True
    assert is_bounded_solution_set(sp.false, (x,)) is True


def test_real_numeric_validator_rejects_finite_complex_values():
    assert is_valid_numeric_value(sp.Integer(2))
    assert is_valid_numeric_value(sp.sqrt(2))
    assert not is_valid_numeric_value(sp.I)
    assert not is_valid_numeric_value(1 + sp.I)
    assert not is_valid_numeric_value(sp.oo)
    assert not is_valid_numeric_value(sp.nan)


def test_symbol_identity_key_distinguishes_assumptions_and_mapping_order():
    xr = sp.Symbol("x", real=True)
    xp = sp.Symbol("x", positive=True)
    assert symbol_identity_key(xr) != symbol_identity_key(xp)

    first = OrderedDict(((xr, sp.Integer(1)), (xp, sp.Integer(2))))
    second = OrderedDict(((xp, sp.Integer(2)), (xr, sp.Integer(1))))
    assert point_key(first) == point_key(second)


def test_quantifier_canonicalization_respects_assumptions():
    xr = sp.Symbol("x", real=True)
    xp = sp.Symbol("x", positive=True)
    q1 = (("exists", xr), ("exists", xp))
    q2 = (("exists", xp), ("exists", xr))
    assert _canonicalize_same_kind_quantifier_blocks(
        q1
    ) == _canonicalize_same_kind_quantifier_blocks(q2)


@settings(max_examples=40, deadline=None)
@given(
    kind=st.sampled_from(["and", "or", "not", "mixed"]),
    left=st.integers(-5, 0),
    right=st.integers(0, 5),
)
def test_coordinate_bounds_complete_tracks_supported_boolean_shape(kind, left, right):
    x = sp.Symbol("x", real=True)
    if kind == "and":
        formula = sp.And(x >= left, x <= right)
        expected = True
    elif kind == "or":
        formula = sp.Or(x <= left, x >= right, evaluate=False)
        expected = False
    elif kind == "not":
        formula = sp.Not(x <= right, evaluate=False)
        expected = False
    else:
        marker = sp.Eq(x**2, 1)
        formula = sp.And(x >= left, marker, evaluate=False)
        expected = False
    assert coordinate_bounds(formula, (x,)).complete is expected


@settings(max_examples=50, deadline=None)
@given(
    family=st.sampled_from(["a", "ab", "a2p1", "a2m1", "square"]),
    aval=st.integers(-4, 4),
    bval=st.integers(-3, 3),
    offset=st.integers(-4, 4),
    sample=st.integers(-6, 6),
)
def test_parameterized_affine_certificate_replays_on_exact_specializations(
    family, aval, bval, offset, sample
):
    a, b, x = sp.symbols("a b x", real=True)
    coeff = {
        "a": a,
        "ab": a * b,
        "a2p1": a**2 + 1,
        "a2m1": a**2 - 1,
        "square": (a - 1) ** 2,
    }[family]
    relation = coeff * x + offset <= 0
    bound = linear_relation_bound(relation, x)
    assignment = {a: aval, b: bval}
    specialized_coeff = sp.simplify(coeff.subs(assignment))
    specialized = relation.subs(assignment)
    if bound is None:
        # An unresolved global parameter sign is deliberately not a certificate.
        return
    specialized_bound = sp.simplify(bound.value.subs(assignment))
    if specialized_coeff == 0:
        # A global certificate must never have divided by a coefficient that
        # vanishes on an admissible specialization.
        raise AssertionError("certified affine bound crossed a zero-coefficient stratum")
    lhs = bool(specialized.subs(x, sample))
    rhs = (
        bool(sample <= specialized_bound)
        if bound.side == "upper"
        else bool(sample >= specialized_bound)
    )
    assert lhs == rhs


@settings(max_examples=30, deadline=None)
@given(reverse_mapping=st.booleans(), reverse_quantifiers=st.booleans())
def test_same_name_symbol_canonicalization_is_permutation_invariant(
    reverse_mapping, reverse_quantifiers
):
    xr = sp.Symbol("x", real=True)
    xp = sp.Symbol("x", positive=True)
    items = [(xr, sp.Integer(1)), (xp, sp.Integer(2))]
    if reverse_mapping:
        items.reverse()
    canonical = point_key(OrderedDict(items))
    assert canonical == point_key(OrderedDict(reversed(items)))

    quantifiers = [("exists", xr), ("exists", xp)]
    if reverse_quantifiers:
        quantifiers.reverse()
    assert _canonicalize_same_kind_quantifier_blocks(tuple(quantifiers)) == (
        ("exists", min((xr, xp), key=symbol_identity_key)),
        ("exists", max((xr, xp), key=symbol_identity_key)),
    )
