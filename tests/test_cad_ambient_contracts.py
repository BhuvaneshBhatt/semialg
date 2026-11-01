"""CAD construction must retain and cover its declared real polynomial statement."""

from dataclasses import replace

import pytest
import sympy as sp

from semialg.cad_algorithms import decomp_collins_complete, decomp_lazard_complete, verify_cad
from semialg.cad_algorithms.projection.collins import build_collins_proj_set
from semialg.cad_algorithms.projection.lazard import build_lazard_proj_set
from semialg.cad_algorithms.projection.lazard_complete import build_lazard_complete_projection
from semialg.cad_algorithms.projection.mccallum import build_mccallum_proj_set


@pytest.mark.parametrize("backend", [decomp_collins_complete, decomp_lazard_complete])
@pytest.mark.parametrize("dimension", [0, 1])
def test_missing_coordinate_rejected(backend, dimension):
    x, y = sp.symbols("x y", real=True)
    with pytest.raises(ValueError, match="unassigned"):
        backend(iter((y,)), (x,) if dimension else ())


@pytest.mark.parametrize(
    "builder",
    [
        build_collins_proj_set,
        build_lazard_complete_projection,
        build_mccallum_proj_set,
        build_lazard_proj_set,
    ],
)
def test_projection_binds_coefficients(builder):
    x, y = sp.symbols("x y", real=True)
    with pytest.raises(ValueError, match="unassigned"):
        builder((y * x,), (x,))


@pytest.mark.parametrize("backend", [decomp_collins_complete, decomp_lazard_complete])
@pytest.mark.parametrize("dimension", [0, 1])
def test_nonreal_coefficients_rejected(backend, dimension):
    x = sp.Symbol("x", real=True)
    for coefficient in (sp.I, sp.oo):
        with pytest.raises(ValueError):
            backend((coefficient * x if dimension else coefficient,), (x,) if dimension else ())


@pytest.mark.parametrize("backend", [decomp_collins_complete, decomp_lazard_complete])
def test_constant_statement_retained(backend):
    source = backend(iter((sp.S(-2), sp.sqrt(2), sp.S.Zero)), ())
    assert tuple(source.tower.original_polynomials) == (-2, sp.sqrt(2), 0)
    assert source.cells_by_level[0][0].sample == ()
    assert source.complete and verify_cad(source)
    malformed = replace(source.tower, original_polynomials=(sp.I,))
    assert not verify_cad(replace(source, tower=malformed))


@pytest.mark.parametrize("backend", [decomp_collins_complete, decomp_lazard_complete])
def test_algebraic_coefficients_remain_supported(backend):
    x = sp.Symbol("x", real=True)
    source = backend((x * x - sp.sqrt(2),), (x,))
    roots = [cell.sample[-1].as_expr() for cell in source.cells if cell.is_section]
    assert len(roots) == 2
    assert all(sp.simplify(root**4 - 2) == 0 for root in roots)
    assert roots[0] < 0 and roots[1] > 0
    assert verify_cad(source)


def test_reduced_iterator_retains_reference_statement():
    from semialg.cad_algorithms.polynomial_utils import polynomial_family_key
    from semialg.cad_algorithms.projection.collins import squarefree_basis
    from semialg.cad_algorithms.reduced import decompose_reduced_safe

    x, y = sp.symbols("x y", real=True)
    family = (y - x, x * x - 2)
    result = decompose_reduced_safe(iter(family), (x, y), backend="mccallum")
    expected = squarefree_basis(sp.Poly(p, x, y) for p in family)
    assert polynomial_family_key(
        result.fallback_cad.tower.original_polynomials
    ) == polynomial_family_key(expected)
    assert result.complete and result.certificate.valid


def test_fallback_reports_actual_backend():
    from unittest.mock import patch

    from semialg.cad_algorithms.reduced import decompose_reduced_safe

    x = sp.Symbol("x", real=True)
    family = (x * x - 2,)
    original = decompose_reduced_safe(family, (x,), backend="mccallum")
    full = decomp_lazard_complete(family, (x,))
    bad_certificate = replace(original.certificate, valid=False)
    bad_side = replace(original.side_conditions, valid=False)
    with (
        patch("semialg.cad_algorithms.reduced.decomp_collins_complete", return_value=full),
        patch(
            "semialg.cad_algorithms.reduced._attempt_reduced_certifi",
            return_value=(original.cad, bad_side, bad_certificate, original.validity),
        ),
        patch("semialg.cad_algorithms.reduced.delin_repair_polys", return_value=()),
    ):
        result = decompose_reduced_safe(family, (x,), backend="mccallum")
    assert result.used_fallback
    assert result.effective_backend == result.cad.backend == "lazard-complete"
    assert verify_cad(result.cad, polynomials=family)


@pytest.mark.parametrize("powers", [(1, 2), (2, 1), (2, 3), (2, 2)])
@pytest.mark.parametrize("scale", [-1, 1])
def test_algebraic_factor_intervals_bind_full_family(powers, scale):
    from semialg.algebraic.roots import isolate_real_roots
    from semialg.algebraic.signs import sign_at_sample

    x = sp.Symbol("x", real=True)
    polynomial = sp.Poly(
        scale * (x + sp.sqrt(2)) ** powers[0] * (x - sp.sqrt(2)) ** powers[1], x, extension=True
    )
    roots = isolate_real_roots(polynomial)
    assert len(roots) == 2
    assert [root.root_index for root in roots] == [0, 1]
    assert [root.multiplicity for root in roots] == list(powers)
    for root in roots:
        assert polynomial.count_roots(root.interval.left, root.interval.right) == 1
    assert sign_at_sample(sp.Poly(x, x), (roots[0],)) == -1
    assert sign_at_sample(sp.Poly(x, x), (roots[1],)) == 1
    assert sign_at_sample(sp.Poly(x + sp.sqrt(2), x, extension=True), (roots[0],)) == 0
    assert sign_at_sample(sp.Poly(x - sp.sqrt(2), x, extension=True), (roots[1],)) == 0
