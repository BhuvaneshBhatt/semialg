from __future__ import annotations

import warnings

import sympy as sp
from sympy.utilities.exceptions import SymPyDeprecationWarning

from semialg import reduce_text


def test_reduce_text_forall_implies_no_boolean_mul_warning():
    b, x = sp.symbols("b x", real=True)
    with warnings.catch_warnings(record=True) as captured:
        warnings.simplefilter("always")
        result = reduce_text(
            "forall x. (-1 <= x and x <= 1) implies x^2 <= b",
            variable_order=[b, x],
        )
    assert not [w for w in captured if issubclass(w.category, SymPyDeprecationWarning)]
    assert result == (b >= 1)
