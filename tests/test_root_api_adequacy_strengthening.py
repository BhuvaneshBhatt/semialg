"""Independent semantic contracts for root APIs found thin by the adequacy audit."""

from __future__ import annotations

import sympy as sp

from semialg import (
    apply_quantifiers,
    discretize_solution,
    find_negative_point,
    prove_negative,
    prove_nonzero,
    prove_zero,
    reduce_formula,
    resolve_formula,
    solve_real_algebraic_set,
    solve_semialgebraic,
)
from semialg.formula import parse_quant_form_text
from semialg.quantifiers import Exists, ForAll


def test_apply_quantifiers_preserves_prefix_order_and_free_parameters() -> None:
    x, y, a = sp.symbols("x y a", real=True)
    result = apply_quantifiers(sp.Eq(x + y, a), (("forall", x), ("exists", y)))

    assert result == ForAll(x, Exists(y, sp.Eq(x + y, a)))
    assert result.free_symbols == {a}


def test_find_negative_point_handles_boundary_minimum_and_certified_absence() -> None:
    x = sp.Symbol("x", real=True)

    witness = find_negative_point((x - 2) ** 2 - 1, (x,))
    assert witness is not None
    assert sp.simplify(((x - 2) ** 2 - 1).subs(witness)) < 0
    assert find_negative_point(x**2, (x,)) is None


def test_thin_sign_provers_cover_true_false_and_structured_contracts() -> None:
    x = sp.Symbol("x", real=True)

    negative = prove_negative(-(x**2 + 1), (x,), return_result=True)
    assert negative.proven is True
    assert negative.relation == "negative"
    assert prove_negative(-(x**2), (x,)) is False

    zero = prove_zero(x - 1, (x,), assumptions=sp.Eq(x, 1), return_result=True)
    assert zero.proven is True
    assert zero.relation == "zero"
    assert prove_zero(x - 1, (x,), assumptions=x > 1) is False

    nonzero = prove_nonzero(x, (x,), assumptions=x < 0, return_result=True)
    assert nonzero.proven is True
    assert nonzero.relation == "nonzero"
    assert prove_nonzero(x, (x,), assumptions=sp.Eq(x, 0)) is False


def test_reduce_and_resolve_formula_cover_false_universal_and_detailed_result() -> None:
    false_universal = parse_quant_form_text("forall x. x^2 > 0")
    true_existential = parse_quant_form_text("exists x. x^2 = 4")

    reduced = reduce_formula(false_universal, return_result=True)
    assert reduced.result is sp.false
    assert resolve_formula(false_universal) is False
    assert reduce_formula(true_existential) is sp.true
    assert resolve_formula(true_existential) is True


def test_real_algebraic_solver_handles_multiplicity_and_inconsistent_system() -> None:
    x, y = sp.symbols("x y", real=True)

    repeated = solve_real_algebraic_set(((x - 2) ** 2,), (x,))
    assert repeated is not None
    assert repeated[x] == 2

    inconsistent = solve_real_algebraic_set((x + y, x + y - 1), (x, y))
    assert inconsistent is None


def test_discretize_solution_root_function_covers_one_dimensional_components() -> None:
    x = sp.Symbol("x", real=True)
    solution = solve_semialgebraic((x >= -1, x <= 2), (x,), count=0)

    data = discretize_solution(solution, bounds=((-2, 3),))
    assert data.dimension == 1
    assert data.source == "interval-components"
    assert len(data.segments) == 1
