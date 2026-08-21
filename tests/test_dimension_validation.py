from __future__ import annotations

import pytest
import sympy as sp

from semialg import (
    DimensionMismatchError,
    require_point_dimension,
    require_same_length,
    sign_at,
    zip_equal,
)


def test_require_same_length_reports_named_dimensions():
    with pytest.raises(
        DimensionMismatchError, match=r"coordinates dimension mismatch .*variables=2.*values=1"
    ):
        require_same_length(
            (sp.Symbol("x"), sp.Symbol("y")),
            (1,),
            context="coordinates",
            names=("variables", "values"),
        )


def test_require_point_dimension_is_public_and_contextual():
    with pytest.raises(DimensionMismatchError, match="optimizer point dimension mismatch"):
        require_point_dimension((1,), (sp.Symbol("x"), sp.Symbol("y")), context="optimizer point")


def test_zip_equal_translates_strict_zip_value_error():
    with pytest.raises(DimensionMismatchError, match="assignment dimension mismatch"):
        list(zip_equal((1, 2), (3,), context="assignment"))


def test_sampling_point_dimension_uses_semialg_error():
    x, y = sp.symbols("x y", real=True)
    with pytest.raises(DimensionMismatchError, match="sampling point dimension mismatch"):
        sign_at(x + y, (1,), variables=(x, y))
