import sympy as sp

from semialg.solve.integer.diophantine import (
    IntEqnSolveResult,
    solve_int_pruning,
    solve_int_sys_via_factor,
)


def test_factorization_solver_preserves_non_equality_constraints():
    x = sp.Symbol("x", integer=True)
    result = solve_int_sys_via_factor(sp.And(sp.Eq((x - 1) * (x - 2), 0), x > 1), (x,))

    assert result is not None
    assert sp.simplify(result.formula.subs(x, 1)) is sp.false
    assert sp.simplify(result.formula.subs(x, 2)) is sp.true


def test_factorization_solver_does_not_claim_complete_with_unresolved_branch(monkeypatch):
    import semialg.solve.integer.diophantine as dio

    x = sp.Symbol("x", integer=True)

    def fake_solver(expr, variables):
        if sp.simplify(expr.subs(x, 1)) is sp.true:
            return IntEqnSolveResult(
                variables=tuple(variables),
                solutions=[(1,)],
                formula=sp.Eq(x, 1),
                method="fake_complete",
                complete=True,
            )
        return None

    monkeypatch.setattr(dio, "solve_int_methods", fake_solver)
    result = solve_int_sys_via_factor(sp.Eq((x - 1) * (x - 2), 0), (x,))

    assert result is not None
    assert not result.complete
    assert sp.simplify(result.formula.subs(x, 1)) is sp.true
    assert sp.simplify(result.formula.subs(x, 2)) is sp.true


def test_bounded_modular_pruning_does_not_claim_completeness_from_witnesses():
    x = sp.Symbol("x", integer=True)
    result = solve_int_pruning(sp.Eq(x, 1), (x,), search_radius=1, moduli=(2, 3))

    assert result is not None
    assert result.solutions
    assert not result.complete
    assert result.metadata["search_radius"] == 1


def test_integer_root_solver_failure_is_not_treated_as_unsat(monkeypatch):
    import semialg.solve.integer.diophantine as dio

    x, y = sp.symbols("x y", integer=True)

    def failed_roots(poly, variable):
        return [], False

    monkeypatch.setattr(dio, "_integer_roots_complete", failed_roots)
    result = dio.solve_int_recursion2(sp.And(sp.Eq(x**2 - 1, 0), sp.Eq(y, 0)), (x, y))

    assert result is None or not result.complete


def test_sum_of_squares_detector_rejects_unconstrained_zero_coefficient_variable():
    from semialg.solve.integer.special_families import detect_sum_fam2

    x, y = sp.symbols("x y", integer=True)
    assert detect_sum_fam2(sp.Eq(x**2, 1), (x, y)) is None


def test_mixed_sign_even_diagonal_family_does_not_claim_false_obstruction():
    from semialg.solve.integer.special_families import solve_diag_fam

    x, y = sp.symbols("x y", integer=True)
    result = solve_diag_fam(sp.Eq(x**2 - y**2, -1), (x, y), search_bound=2)

    assert result is not None
    assert result.formula is not sp.false
    assert not result.complete
    assert (0, 1) in result.solutions or (0, -1) in result.solutions
