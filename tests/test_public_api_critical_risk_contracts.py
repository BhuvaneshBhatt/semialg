"""Independent behavioral contracts for thin critical root APIs."""

from __future__ import annotations

import sympy as sp

from semialg import (
    apply_quantifiers,
    find_negative_point,
    find_negative_witness_fast,
    polynomial_nonnegative,
    prove_nonzero,
    prove_zero,
    resolve_formula,
    solve_real_algebraic_set,
)
from semialg.formula import parse_quant_form_text


def test_quantifier_application_and_resolution_agree_on_closed_sentences() -> None:
    x = sp.Symbol("x", real=True)
    existential = apply_quantifiers(sp.Eq(x**2, 2), (("exists", x),))
    universal = apply_quantifiers(x**2 >= 0, (("forall", x),))

    assert existential.free_symbols == set()
    assert universal.free_symbols == set()
    assert resolve_formula(parse_quant_form_text("exists x. x^2 = 2")) is True
    assert resolve_formula(parse_quant_form_text("forall x. x^2 >= 0")) is True
    assert resolve_formula(parse_quant_form_text("exists x. x^2 = -1")) is False


def test_negative_point_apis_verify_witnesses() -> None:
    x, y = sp.symbols("x y", real=True)
    polynomial = x**2 + y**2 - 1

    exact = find_negative_point(polynomial, (x, y))
    assert exact is not None
    assert sp.simplify(polynomial.subs(exact)) < 0

    fast = find_negative_witness_fast(polynomial, (x, y), seed=19)
    assert fast.assignment is not None
    assert sp.simplify(polynomial.subs(fast.assignment)) < 0

    assert find_negative_point(x**2 + y**2 + 1, (x, y)) is None
    no_fast_witness = find_negative_witness_fast(x**2 + y**2 + 1, (x, y), seed=19)
    assert no_fast_witness.assignment is None


def test_polynomial_nonnegative_distinguishes_global_signs() -> None:
    x = sp.Symbol("x", real=True)

    assert polynomial_nonnegative((x - 3) ** 2, (x,)) is True
    assert polynomial_nonnegative(x**2 - 1, (x,)) is False


def test_zero_and_nonzero_proofs_respect_assumptions_and_structured_results() -> None:
    x = sp.Symbol("x", real=True)

    zero = prove_zero((x - 1) * (x + 1), (x,), assumptions=sp.Eq(x**2, 1), return_result=True)
    nonzero = prove_nonzero(x, (x,), assumptions=x > 0, return_result=True)

    assert zero.proven is True
    assert zero.relation == "zero"
    assert nonzero.proven is True
    assert nonzero.relation == "nonzero"
    assert prove_zero(x, (x,), assumptions=x > 0) is False
    assert prove_nonzero(x, (x,), assumptions=sp.Eq(x, 0)) is False


def test_real_algebraic_solver_returns_verified_witness_or_none() -> None:
    x, y = sp.symbols("x y", real=True)

    witness = solve_real_algebraic_set((x + y - 1, x - y), (x, y))
    assert witness is not None
    assert all(sp.simplify(poly.subs(witness)) == 0 for poly in (x + y - 1, x - y))
    assert solve_real_algebraic_set((x**2 + 1,), (x,)) is None
