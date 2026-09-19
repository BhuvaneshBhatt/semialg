from __future__ import annotations

from collections.abc import Mapping, Sequence

import sympy as sp


def _modular_reduce_value(value, modulus: int):
    value = sp.sympify(value)
    if bool(value.is_integer):
        return int(value) % modulus
    num, den = sp.fraction(sp.together(value))
    if not (num.is_integer and den.is_integer):
        raise ValueError("non-rational value in modular check")
    return (int(num) * pow(int(den), -1, modulus)) % modulus


def form_sat_by_assign(
    formula: sp.Expr,
    assignment: Mapping[sp.Symbol, object],
    *,
    domain=sp.Complexes,
    check_numeric_equalities: bool = False,
    modulus: int | None = None,
) -> bool:
    substituted = formula.subs(dict(assignment))
    if modulus is not None:
        substituted = substituted.replace(
            lambda e: isinstance(e, sp.Equality),
            lambda e: (
                sp.true
                if _modular_reduce_value(e.lhs, modulus) == _modular_reduce_value(e.rhs, modulus)
                else sp.false
            ),
        ).replace(
            lambda e: isinstance(e, sp.Unequality),
            lambda e: (
                sp.true
                if _modular_reduce_value(e.lhs, modulus) != _modular_reduce_value(e.rhs, modulus)
                else sp.false
            ),
        )
    if check_numeric_equalities:
        substituted = sp.simplify(substituted)
    return substituted is sp.true or substituted is True


def sample_assigns_sat_form(
    formula: sp.Expr,
    assignments: Sequence[Mapping[sp.Symbol, object]],
    *,
    domain=sp.Complexes,
    modulus: int | None = None,
    check_numeric_equalities: bool = True,
) -> list[bool]:
    return [
        form_sat_by_assign(
            formula,
            a,
            domain=domain,
            modulus=modulus,
            check_numeric_equalities=check_numeric_equalities,
        )
        for a in assignments
    ]


__all__ = ["form_sat_by_assign", "sample_assigns_sat_form"]
