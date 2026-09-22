"""Production-path contracts for public types that users normally receive indirectly."""

from __future__ import annotations

import pytest
import sympy as sp

import semialg

pytestmark = pytest.mark.contract


def test_cad_region_is_obtained_through_as_cad_region() -> None:
    x = sp.Symbol("x", real=True)

    region = semialg.as_cad_region(sp.true, (x,))

    assert isinstance(region, semialg.CADRegion)
    assert region.variables == (x,)
    assert region.formula is sp.true
    signature = region.signature()
    assert signature.ambient_dimension == 1
    assert signature.region_dimension == 1
    assert signature.selected_cell_count == 1
    assert region.as_region().contains((sp.Integer(0),))


def test_affine_box_clip_is_obtained_through_clipping_operation() -> None:
    clip = semialg.clip_affine_subspace_to_box(
        (0, 0),
        ((1, 1),),
        ((0, 1), (0, 1)),
    )

    assert isinstance(clip, semialg.AffineBoxClip)
    assert clip.vertices == ((sp.Integer(0), sp.Integer(0)), (sp.Integer(1), sp.Integer(1)))
    assert clip.dimension == 1
    assert clip.parameter_dimension == 1


def test_geometry_and_standard_region_are_exercised_through_concrete_region() -> None:
    region = semialg.Interval(0, 1)

    assert isinstance(region, semialg.Geometry)
    assert isinstance(region, semialg.StandardRegion)
    assert region.dimension() == 1
    assert region.ambient_dimension() == 1
    assert region.contains((sp.Rational(1, 2),))
    assert not region.contains((sp.Integer(2),))


def test_resource_limit_raising_api_exercises_public_exception_bases() -> None:
    x = sp.Symbol("x", real=True)

    with pytest.raises(semialg.ResourceLimitError) as caught:
        semialg.cad(sp.Abs(x) <= 1, (x,), max_preprocess_aux_vars=0)

    error = caught.value
    assert isinstance(error, semialg.SemialgStrategyFailure)
    assert isinstance(error, semialg.SemialgError)
    assert "auxiliary" in str(error).lower() or "limit" in str(error).lower()


def test_unsupported_fragment_exception_keeps_explicit_no_raising_site_contract() -> None:
    # There is intentionally no package production raising site for this public
    # exception yet.  Keep construction/inheritance coverage explicit rather
    # than inventing an artificial production path.
    error = semialg.UnsupportedFragmentError("unsupported symbolic fragment")

    assert isinstance(error, semialg.SemialgStrategyFailure)
    assert isinstance(error, semialg.SemialgError)
    assert isinstance(error, ValueError)
    assert "unsupported symbolic fragment" in str(error)
