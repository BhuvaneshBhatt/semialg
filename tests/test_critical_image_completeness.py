"""An unavailable critical image is not an empty value set."""

from unittest.mock import patch

import pytest
import sympy as sp

from semialg import critical_value_image
from semialg.errors import ExactComputationUnavailable


def test_default_output_collision():
    t, y = sp.symbols("t y", real=True)
    result = critical_value_image(t, sp.Eq(y * y, 0), (t, y))
    assert result.value_symbol not in (t, y)
    assert result.formula is sp.true
    assert result.component_images == (sp.true,)


def test_explicit_output_collision():
    x, y = sp.symbols("x y", real=True)
    with pytest.raises(ValueError, match="distinct"):
        critical_value_image(x, sp.Eq(y * y, 0), (x, y), value_symbol=x)


def test_image_decline_is_not_empty():
    x, y = sp.symbols("x y", real=True)
    with patch(
        "semialg.optimization.function_range",
        side_effect=NotImplementedError("forced image decline"),
    ):
        with pytest.raises(ExactComputationUnavailable, match="component image"):
            critical_value_image(x, sp.Eq(y * y, 0), (x, y))


def test_locus_decline_is_not_empty():
    x, y = sp.symbols("x y", real=True)
    with patch("semialg._critical_loci.polynomial_locus_dimension", return_value=None):
        with pytest.raises(ExactComputationUnavailable, match="locus dimension"):
            critical_value_image(x, sp.Eq(y * y, 0), (x, y))
