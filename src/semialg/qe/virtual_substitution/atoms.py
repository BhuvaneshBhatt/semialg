from __future__ import annotations

from collections.abc import Iterable

import sympy as sp
from sympy.core.relational import Relational

from ...formulas.boolean import (
    canonical_relation as _shared_canonical_relation,
)
from ...formulas.boolean import (
    iter_relational_atoms as _shared_iter_relational_atoms,
)
from ...formulas.boolean import (
    make_and as _shared_and,
)
from ...formulas.boolean import (
    make_or as _shared_or,
)
from ...formulas.boolean import (
    negate_relation as _shared_negate_relation,
)
from ...formulas.boolean import (
    relation_from_residual as _shared_relation,
)
from ...formulas.boolean import (
    to_negation_normal_form as _shared_to_negation_normal_form,
)
from .types import VirtualSubstitutionError


def _and(*args: sp.Expr) -> sp.Expr:
    return _shared_and(*args)


def _or(*args: sp.Expr) -> sp.Expr:
    return _shared_or(*args)


def _relation(expr: sp.Expr, operator: str) -> sp.Expr:
    try:
        return _shared_relation(expr, operator)
    except ValueError as exc:
        raise VirtualSubstitutionError(str(exc)) from exc


def _canonical_atom(atom: Relational) -> tuple[sp.Expr, str]:
    try:
        return _shared_canonical_relation(atom)
    except ValueError as exc:
        raise VirtualSubstitutionError(str(exc)) from exc


def _negate_atom(atom: Relational) -> sp.Expr:
    try:
        return _shared_negate_relation(atom)
    except ValueError as exc:
        raise VirtualSubstitutionError(str(exc)) from exc


def _to_negation_normal_form(formula: sp.Expr, negate: bool = False) -> sp.Expr:
    try:
        return _shared_to_negation_normal_form(formula, negate=negate)
    except ValueError as exc:
        raise VirtualSubstitutionError(str(exc)) from exc


def _iter_atoms(formula: sp.Expr) -> Iterable[Relational]:
    try:
        yield from _shared_iter_relational_atoms(formula)
    except ValueError as exc:
        raise VirtualSubstitutionError(str(exc)) from exc


def _polynomial_degree(expr: sp.Expr, variable: sp.Symbol) -> int:
    try:
        poly = sp.Poly(sp.expand(expr), variable, domain="EX")
    except sp.PolynomialError as exc:
        raise VirtualSubstitutionError(f"non-polynomial atom in {variable}: {expr!r}") from exc
    return poly.degree()


def _coefficient_list(expr: sp.Expr, variable: sp.Symbol, length: int) -> list[sp.Expr]:
    poly = sp.Poly(sp.expand(expr), variable, domain="EX")
    coeffs = [sp.expand(poly.nth(i)) for i in range(length)]
    return coeffs
