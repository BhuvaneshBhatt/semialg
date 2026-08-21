from __future__ import annotations

import pytest
import sympy as sp

from semialg.cad.bounds import (
    AlgebraicRootFunction,
    DelineabilityCertificate,
    ExplicitCADBound,
    verify_cad_cell_bounds,
)
from semialg.cad.cells import (
    CylindricalCoordinateConstraint,
    CylindricalDecompositionCertificate,
    CylindricalSolution,
    CylindricalSolutionCell,
    extract_cylindrical_solution,
    extract_explicit_cylindrical_solution,
)
from semialg.cad.integration import full_dimensional_solution_integrals, intrinsic_cell_integral


def test_negative_leading_quadratic_uses_correct_cad_root_branch():
    x, y = sp.symbols("x y", real=True)
    solution = extract_cylindrical_solution(sp.And(x > -2, x < -1, x * y**2 + y > 0), (x, y))
    cell = solution.full_dimensional_cells[0]
    y_level = cell.levels[1]
    assert sp.simplify(y_level.lower) == 0
    assert sp.simplify(y_level.upper + 1 / x) == 0
    assert verify_cad_cell_bounds(cell).verify()


def test_verifier_rejects_wrong_radical_presentation_even_with_root_certificate():
    x, y = sp.symbols("x y", real=True)
    # On x=-3/2, root 0 of x*y**2+y is 0.  Branch 1 is -1/x=2/3.
    cert = DelineabilityCertificate(
        polynomial=x * y**2 + y,
        fiber_variable=y,
        root_index=0,
        base_variables=(x,),
        sign_invariant=True,
        stack_order_verified=True,
        sample_root_verified=True,
        sample_root_value=sp.Integer(0),
        radical_branch_index=0,
        representation_verified=True,
        regular_section_verified=True,
    )
    wrong = AlgebraicRootFunction(
        x * y**2 + y,
        y,
        0,
        base_variables=(x,),
        certificate=cert,
        closed=True,
    )
    x_level = CylindricalCoordinateConstraint(
        x,
        1,
        "section",
        sp.Rational(-3, 2),
        sp.Rational(-3, 2),
        sp.Rational(-3, 2),
        (1,),
        lower_bound=ExplicitCADBound(sp.Rational(-3, 2), True),
        upper_bound=ExplicitCADBound(sp.Rational(-3, 2), True),
        lower_closed=True,
        upper_closed=True,
    )
    y_level = CylindricalCoordinateConstraint(
        y,
        2,
        "section",
        wrong.as_expr(),
        wrong.as_expr(),
        0,
        (1, 1),
        lower_bound=wrong,
        upper_bound=wrong,
        lower_closed=True,
        upper_closed=True,
        delineability=cert,
    )
    cell = CylindricalSolutionCell(
        (x, y), (x_level, y_level), {x: sp.Rational(-3, 2), y: 0}, (1, 1)
    )
    assert not verify_cad_cell_bounds(cell).verify()


def test_explicit_fast_path_declines_incomparable_symbolic_bounds():
    x, a, b = sp.symbols("x a b", real=True)
    result = extract_explicit_cylindrical_solution(sp.And(x >= a, x >= b, x <= 10), (x,))
    assert result is None


def test_explicit_or_declines_overlapping_pieces_to_avoid_double_counting():
    x = sp.symbols("x", real=True)
    formula = sp.Or(sp.And(x >= 0, x <= 2), sp.And(x >= 1, x <= 3))
    assert extract_explicit_cylindrical_solution(formula, (x,)) is None


def test_explicit_or_accepts_provably_disjoint_pieces_and_certifies_decomposition():
    x = sp.symbols("x", real=True)
    formula = sp.Or(sp.And(x >= 0, x <= 1), sp.And(x > 1, x <= 2))
    solution = extract_explicit_cylindrical_solution(formula, (x,))
    assert solution is not None
    assert solution.decomposition_cert is not None
    assert solution.decomposition_cert.verify()
    values = full_dimensional_solution_integrals(solution, 1, evaluate=True)
    assert sp.simplify(sum(item.integral for item in values) - 2) == 0


def test_explicit_fast_path_rejects_empty_and_open_singleton_cells():
    x = sp.symbols("x", real=True)
    assert extract_explicit_cylindrical_solution(sp.And(x >= 2, x <= 1), (x,)) is None
    open_empty = extract_explicit_cylindrical_solution(sp.And(x > 1, x <= 1), (x,))
    assert open_empty is not None and open_empty.cells == ()
    assert open_empty.decomposition_cert is not None and open_empty.decomposition_cert.verify()
    point = extract_explicit_cylindrical_solution(sp.And(x >= 1, x <= 1), (x,))
    assert point is not None
    assert point.cells[0].levels[0].is_section


def test_intrinsic_integration_rejects_uncertified_singular_algebraic_section():
    y = sp.symbols("y", real=True)
    cert = DelineabilityCertificate(
        polynomial=y**2,
        fiber_variable=y,
        root_index=0,
        sign_invariant=True,
        stack_order_verified=True,
        sample_root_verified=True,
        sample_root_value=sp.Integer(0),
        representation_verified=False,
        regular_section_verified=False,
    )
    root = AlgebraicRootFunction(y**2, y, 0, certificate=cert, closed=True)
    level = CylindricalCoordinateConstraint(
        y,
        1,
        "section",
        0,
        0,
        0,
        (1,),
        lower_bound=root,
        upper_bound=root,
        lower_closed=True,
        upper_closed=True,
        delineability=cert,
    )
    cell = CylindricalSolutionCell((y,), (level,), {y: 0}, (1,))
    assert verify_cad_cell_bounds(cell).verify()
    with pytest.raises(ValueError, match="regular algebraic sections"):
        intrinsic_cell_integral(cell, require_verified=True)


def test_solution_integrators_reject_failed_decomposition_cert():
    x = sp.symbols("x", real=True)
    base = extract_explicit_cylindrical_solution(sp.And(x >= 0, x <= 1), (x,))
    assert base is not None
    failed = CylindricalSolution(
        base.variables,
        base.cells,
        base.formula,
        decomposition_cert=CylindricalDecompositionCertificate(True, False, True, "test"),
    )
    with pytest.raises(ValueError, match="decomposition"):
        full_dimensional_solution_integrals(failed)


def test_programming_errors_in_root_presentation_are_not_swallowed(monkeypatch):
    import semialg.reconstruct.radicals as radicals

    x, y = sp.symbols("x y", real=True)
    cert = DelineabilityCertificate(
        polynomial=y**2 - x,
        fiber_variable=y,
        root_index=0,
        base_variables=(x,),
        sign_invariant=True,
        stack_order_verified=True,
        sample_root_verified=True,
        radical_branch_index=0,
        representation_verified=True,
    )
    root = AlgebraicRootFunction(y**2 - x, y, 0, (x,), certificate=cert)

    def broken(*args, **kwargs):
        raise NameError("programmer bug")

    monkeypatch.setattr(radicals, "fiber_root_candidates", broken)
    with pytest.raises(NameError, match="programmer bug"):
        root.as_expr()


def test_pytest_configuration_supports_src_layout(pytestconfig):
    configured = [str(path) for path in pytestconfig.getini("pythonpath")]
    assert any(path.endswith("src") for path in configured)
