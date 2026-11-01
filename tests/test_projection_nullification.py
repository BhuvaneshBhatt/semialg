"""Order-invariance tests at nullified projection fibers."""

from unittest.mock import patch

import pytest
import sympy as sp

from semialg.algebraic.samples import RationalSample
from semialg.cache_control import clear_caches
from semialg.cad_algorithms.decomposition import (
    _build_stack,
    _stack_roots_over_point,
    decomp_collins_complete,
)
from semialg.cad_algorithms.lifting.stack import CADCell
from semialg.cad_algorithms.projection.collins import _squarefree_factors, build_collins_proj_set
from semialg.errors import ExactComputationUnavailable
from semialg.formula import parse_formula
from semialg.partial.qe import lazy_resolve_formula
from semialg.qe.complete import evaluate_formula_on_cell


def test_isolated_nullification_delineates_root_collision():
    x, y, z, w = sp.symbols("x y z w", real=True)
    variables = (x, y, z, w)
    polynomial = ((w - 1) ** 2 - (z - 1)) ** 2 + x * w + y
    tower = build_collins_proj_set((polynomial,), variables)
    zero = RationalSample(sp.S.Zero)
    parent = CADCell(level=2, index=(1, 1), sample=(zero, zero))
    roots = _stack_roots_over_point(
        tower.level(3).polynomials,
        variables,
        parent.sample,
        z,
        delineate_nullification=True,
    )
    assert tuple(root.as_expr() for root in roots) == (1, 2, 4)
    cells = _build_stack(parent, roots, 3, tower)
    collision = next(cell for cell in cells if cell.is_section and cell.sample[-1].as_expr() == 1)
    assert sp.expand(collision.section_polynomial).subs({x: 0, y: 0, z: 1}) == 0
    assert sp.expand(collision.section_polynomial).subs({x: 0, y: 0, z: 0}) != 0
    for cell in cells:
        value = cell.sample[-1].as_expr()
        fiber_roots = _stack_roots_over_point(tower.level(4).polynomials, variables, cell.sample, w)
        assert len(fiber_roots) == (0 if value < 1 else 1 if value == 1 else 2)


def test_positive_dimensional_nullification_requires_refinement():
    x, y, t, z = sp.symbols("x y t z", real=True)
    polynomial = sp.Poly(x**2 * (z - 1) ** 2 + y**3, x, y, t, z)
    with pytest.raises(ExactComputationUnavailable):
        _stack_roots_over_point(
            (polynomial,),
            (x, y, t, z),
            (RationalSample(sp.S.Zero),) * 3,
            z,
            delineate_nullification=True,
            base_is_point=False,
        )


def test_algebraic_isolated_nullification():
    from semialg.algebraic.roots import isolate_real_roots

    x, z = sp.symbols("x z", real=True)
    alpha = isolate_real_roots(sp.Poly(x**2 - 2, x))[-1]
    polynomial = sp.Poly((x**2 - 2) ** 2 * (z**2 - 3) + (x**2 - 2) ** 3, x, z)
    roots = _stack_roots_over_point(
        (polynomial,), (x, z), (alpha,), z, delineate_nullification=True
    )
    assert tuple(root.as_expr() for root in roots) == (-sp.sqrt(3), sp.sqrt(3))


def test_lazy_qe_finds_missing_collision_boundary():
    x, y, z, w = sp.symbols("x y z w", real=True)
    polynomial = ((w - 1) ** 2 - (z - 1)) ** 2 + x * w + y
    matrix = parse_formula(sp.And(sp.Eq(x, 0), sp.Eq(y, 0), sp.Eq(polynomial, 0)))
    result = lazy_resolve_formula(
        (x, y, z, w),
        (("exists", x), ("exists", y), ("forall", z), ("exists", w)),
        matrix,
    )
    # At x=y=0, real roots exist exactly for z>=1, so the universal
    # statement fails. Sampling z=1 for the old unsplit cell z<2 hid z<1.
    assert result.truth_value is False
    assert result.status == "complete"


def test_complete_cad_delineates_isolated_projection_zero():
    x, y, z, w = sp.symbols("x y z w", real=True)
    source = decomp_collins_complete((w**2 + y * z - x**2,), (x, y, z, w))
    cells = [cell for cell in source.cells_by_level[3] if cell.index[:2] == (1, 1)]
    assert [cell.sample_exprs[-1] for cell in cells] == [-1, 0, 1]
    assert [cell.kind for cell in cells] == ["sector", "section", "sector"]
    for value in (-2, 0, 2):
        section = source._point_locator.locate((0, 0, value, 0))
        assert section.is_section
        assert evaluate_formula_on_cell(
            parse_formula(sp.Eq(w**2 + y * z - x**2, 0)), section, (x, y, z, w)
        )


def test_factorization_failure_preserves_squarefreeness():
    x, y = sp.symbols("x y", real=True)
    polynomial = sp.Poly((y - x) ** 2 * (y + x - 1) ** 3, x, y)
    with patch("sympy.factor_list", side_effect=NotImplementedError):
        basis = _squarefree_factors(polynomial)
    expected = sp.Poly((y - x) * (y + x - 1), x, y)
    assert len(basis) == 1
    assert basis[0].monic() == expected.monic()


def test_squarefree_failure_is_explicit():
    x, y = sp.symbols("x y", real=True)
    polynomial = sp.Poly((y - x) ** 2, x, y)
    clear_caches(collect=False)
    with (
        patch.object(sp.Poly, "sqf_part", side_effect=NotImplementedError),
        patch("sympy.factor_list", side_effect=NotImplementedError),
        pytest.raises(ExactComputationUnavailable),
    ):
        _squarefree_factors(polynomial)
