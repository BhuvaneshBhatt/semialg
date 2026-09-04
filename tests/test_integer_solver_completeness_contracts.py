import sympy as sp

from semialg.solve.integer.diophantine import (
    IntEqnSolveResult,
    solve_int_sys_via_factor,
    solve_integer_with_modular_pruning,
)


def test_factorization_solver_preserves_non_equality_constraints():
    x = sp.Symbol("x", integer=True)
    result = solve_int_sys_via_factor(sp.And(sp.Eq((x - 1) * (x - 2), 0), x > 1), (x,))

    assert result is not None
    assert sp.simplify(result.formula.subs(x, 1)) is sp.false
    assert sp.simplify(result.formula.subs(x, 2)) is sp.true


def test_factorization_solver_does_not_claim_complete_with_unresolved_branch():

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

    result = solve_int_sys_via_factor(
        sp.Eq((x - 1) * (x - 2), 0),
        (x,),
        branch_solver=fake_solver,
    )

    assert result is not None
    assert not result.complete
    assert sp.simplify(result.formula.subs(x, 1)) is sp.true
    assert sp.simplify(result.formula.subs(x, 2)) is sp.true


def test_bounded_modular_pruning_does_not_claim_completeness_from_witnesses():
    x = sp.Symbol("x", integer=True)
    result = solve_integer_with_modular_pruning(sp.Eq(x, 1), (x,), search_radius=1, moduli=(2, 3))

    assert result is not None
    assert result.solutions
    assert not result.complete
    assert result.metadata["search_radius"] == 1


def test_integer_root_solver_failure_is_not_treated_as_unsat():
    import semialg.solve.integer.diophantine as dio

    x, y = sp.symbols("x y", integer=True)

    def failed_roots(poly, variable):
        return [], False

    result = dio.solve_recursive_diophantine(
        sp.And(sp.Eq(x**2 - 1, 0), sp.Eq(y, 0)),
        (x, y),
        root_solver=failed_roots,
    )

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


def test_recursive_groebner_preserves_symbol_assumptions_in_cache_keys():
    from semialg.solve.integer.groebner_recursion import compute_groebner_basis

    x, y = sp.symbols("x y", integer=True)
    expr = sp.And(sp.Eq(x**2 - 1, 0), sp.Eq(y - x, 0))

    basis = compute_groebner_basis(expr, (x, y))

    assert basis is not None
    assert sp.Integer(1) not in basis
    assert all(poly.free_symbols <= {x, y} for poly in basis)


def test_recursive_groebner_uses_acyclic_linear_elimination():
    from semialg.solve.integer.groebner_recursion import rec_reduce_sys

    x, y = sp.symbols("x y", integer=True)
    expr = sp.And(sp.Eq(x**2 - 1, 0), sp.Eq(y - x, 0))

    result = rec_reduce_sys(expr, (x, y))

    assert result is not None
    assert result.complete
    assert set(result.solutions) == {(-1, -1), (1, 1)}


def test_recursive_groebner_branches_without_coupled_relation():
    from semialg.solve.integer.groebner_recursion import rec_reduce_sys

    x, y = sp.symbols("x y", integer=True)
    expr = sp.And(sp.Eq(x**2 - 1, 0), sp.Eq(y**2 - 4, 0))

    result = rec_reduce_sys(expr, (x, y))

    assert result is not None
    for point in ((-1, -2), (-1, 2), (1, -2), (1, 2)):
        assignment = dict(zip((x, y), point, strict=True))
        assert sp.simplify(result.formula.subs(assignment)) is sp.true
