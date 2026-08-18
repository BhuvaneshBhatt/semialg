from __future__ import annotations

import pytest
import sympy as sp

from semialg.cad.reduced import ReducedCertificate, decomp_form_reduced_safe, decompose_reduced_safe
from semialg.formula import parse_formula
from semialg.tticad.safe import decompose_tticad_safe

pytestmark = pytest.mark.slow


def test_active_reduced_01():
    x, y = sp.symbols("x y")
    formula = parse_formula(sp.Eq(y**2 - x, 0))

    result = decomp_form_reduced_safe(formula, [x, y], backend="mccallum")

    assert result.complete is True
    assert result.used_fallback is False
    assert result.effective_backend == "mccallum-reduced-certified"
    assert isinstance(result.certificate, ReducedCertificate)
    assert result.certificate.valid is True
    assert result.validity.valid is True
    assert result.validity.complete_if_used is True
    assert result.fallback_cad is not None


def test_active_reduced_02():
    x, y = sp.symbols("x y")
    formula = parse_formula(sp.Eq(y**2 - x, 0))

    result = decomp_form_reduced_safe(formula, [x, y], backend="lazard")

    assert result.complete is True
    assert result.used_fallback is False
    assert result.effective_backend == "lazard-reduced-certified"
    assert result.certificate is not None
    assert result.certificate.valid is True
    assert result.certificate.invariant == "truth"
    assert result.fallback_cad is not None


def test_active_reduced_03():
    x = sp.symbols("x")

    result = decompose_reduced_safe([x**2 - 1], [x], backend="mccallum")

    assert result.complete is True
    assert result.used_fallback is False
    assert result.certificate is not None
    assert result.certificate.valid is True
    assert result.certificate.invariant == "sign"


def test_active_reduced_04():
    x, y = sp.symbols("x y")
    formula = parse_formula(sp.Eq(y**2 - x, 0))

    result = decompose_tticad_safe(formula, [x, y])

    assert result.complete is True
    assert result.used_fallback is False
    assert result.effective_backend == "tticad-reduced-certified"
    assert result.certificate is not None
    assert result.certificate.valid is True
    assert result.validity.valid is True
    assert result.projection_validity.valid is True
    assert result.fallback_cad is not None
