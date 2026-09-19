"""Semantic helpers for CAD tests.

These helpers deliberately avoid asserting that a source polynomial is stored
verbatim in a projection/sign table. CAD may replace a reducible polynomial by
an exact squarefree factor basis without changing its sign semantics.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence

import sympy as sp

from semialg.cad_algorithms.polynomial_utils import polynomial_key


def sign_from_factor_signs(
    polynomial: sp.Expr,
    variables: Sequence[sp.Symbol],
    signs: Mapping[str, int],
) -> int:
    """Recover a source polynomial sign from recorded squarefree-factor signs."""
    coeff, factors = sp.factor_list(sp.expand(polynomial), *tuple(variables))
    coeff_sign = int(sp.sign(coeff))
    if coeff_sign == 0:
        return 0
    result = coeff_sign
    for factor, multiplicity in factors:
        poly = sp.Poly(factor, *tuple(variables))
        key = polynomial_key(poly)
        factor_sign = signs[key]
        if factor_sign == 0:
            return 0
        if multiplicity % 2:
            result *= factor_sign
    return result
