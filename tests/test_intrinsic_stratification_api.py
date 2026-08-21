from __future__ import annotations

from dataclasses import replace

import pytest
import sympy as sp

from semialg import (
    DelineabilityCertificate,
    IntrinsicCellStratum,
    IntrinsicStratification,
    extract_explicit_cylindrical_solution,
    intrinsic_solution_integrals,
    stratify_intrinsic_solution,
)


def test_intrinsic_stratum_properties_and_stratification_properties():
    regular = IntrinsicCellStratum(object(), (1,), 1, True)
    singular = IntrinsicCellStratum(object(), (2,), 0, False, ("singular",))
    result = IntrinsicStratification((regular, singular), target_dimension=None)
    assert not regular.singular
    assert singular.singular
    assert result.regular_strata == (regular,)
    assert result.singular_strata == (singular,)
    assert len(result.regular_cells) == 1
    assert len(result.singular_cells) == 1


def test_regular_explicit_graph_is_regular_stratum_and_integrates():
    x, y = sp.symbols("x y", real=True)
    sol = extract_explicit_cylindrical_solution(sp.And(x >= 0, x <= 1, sp.Eq(y, x)), [x, y])
    assert sol is not None
    strat = stratify_intrinsic_solution(sol)
    assert len(strat.regular_strata) == 1
    assert not strat.singular_strata
    integrals = intrinsic_solution_integrals(sol, 1, dimension=1, evaluate=True)
    assert len(integrals) == 1
    assert sp.simplify(integrals[0].integral - sp.sqrt(2)) == 0


def test_uncertified_algebraic_section_is_exposed_as_singular_stratum():
    x, y = sp.symbols("x y", real=True)
    sol = extract_explicit_cylindrical_solution(sp.And(x >= 0, x <= 1, sp.Eq(y, x)), [x, y])
    assert sol is not None
    cell = sol.cells[0]
    cert = DelineabilityCertificate(
        polynomial=(y - x) ** 2,
        fiber_variable=y,
        root_index=0,
        base_variables=(x,),
        sign_invariant=True,
        stack_order_verified=True,
        sample_root_verified=True,
        sample_root_value=cell.sample[y],
        regular_section_verified=False,
    )
    ylevel = replace(cell.levels[1], delineability=cert)
    singular_cell = replace(cell, levels=(cell.levels[0], ylevel))
    singular_sol = replace(sol, cells=(singular_cell,))
    strat = stratify_intrinsic_solution(singular_sol)
    assert not strat.regular_strata
    assert len(strat.singular_strata) == 1
    assert "regularity certificate" in strat.singular_strata[0].reasons[0]
    with pytest.raises(ValueError, match="singular strata"):
        intrinsic_solution_integrals(singular_sol, dimension=1)


def test_dimension_filter_in_stratification():
    x, y = sp.symbols("x y", real=True)
    sol = extract_explicit_cylindrical_solution(sp.And(x >= 0, x <= 1, sp.Eq(y, x)), [x, y])
    assert sol is not None
    assert len(stratify_intrinsic_solution(sol, dimension=1).strata) == 1
    assert not stratify_intrinsic_solution(sol, dimension=0).strata


def test_public_intrinsic_stratification_exports_are_importable():
    import semialg

    for name in (
        "IntrinsicCellStratum",
        "IntrinsicStratification",
        "stratify_intrinsic_solution",
        "intrinsic_cell_integral",
        "intrinsic_solution_integrals",
    ):
        assert name in semialg.__all__
        assert getattr(semialg, name) is not None


def test_algebraic_section_without_certificate_is_never_classified_regular() -> None:
    from dataclasses import replace

    from semialg import AlgebraicRootFunction

    x, y = sp.symbols("x y", real=True)
    sol = extract_explicit_cylindrical_solution(sp.And(x >= 0, x <= 1, sp.Eq(y, x)), [x, y])
    assert sol is not None
    cell = sol.cells[0]
    root = AlgebraicRootFunction(y - x, y, 0, base_variables=(x,), certificate=None)
    ylevel = replace(
        cell.levels[1],
        lower=x,
        upper=x,
        lower_bound=root,
        upper_bound=root,
        delineability=None,
    )
    uncertified_cell = replace(cell, levels=(cell.levels[0], ylevel))
    uncertified_sol = replace(sol, cells=(uncertified_cell,))
    strat = stratify_intrinsic_solution(uncertified_sol, require_verified=False)
    assert len(strat.singular_strata) == 1
    assert not strat.regular_strata
    assert "regularity certificate" in strat.singular_strata[0].reasons[0]
