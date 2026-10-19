"""Migration preserves certified identity versus pointwise decisions."""

from unittest.mock import patch

import sympy as sp
from exprtools import zerotest

import semialg._zero_testing as boundary


def test_boundary_uses_exprtools_and_retains_parameter_strata():
    assert boundary.zerotest is zerotest
    x = sp.Symbol("x", real=True)
    assert boundary.certified_zero(x) is False
    assert boundary.certified_pointwise_zero(x) is None
    assert boundary.certified_pointwise_zero(x, assumptions=sp.Q.positive(x)) is False


def test_unknown_is_preserved_and_certified_mode_is_forwarded():
    x = sp.Symbol("x")
    with patch.object(boundary, "zerotest", return_value=None) as oracle:
        assert boundary.certified_pointwise_zero(sp.sin(x)) is None
    oracle.assert_called_once_with(sp.sin(x), assumptions=True, confidence="certified")
