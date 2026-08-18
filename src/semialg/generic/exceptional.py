from __future__ import annotations

from collections.abc import Iterable, Sequence
from dataclasses import dataclass

import sympy as sp

from ..cad.decomposition import CompleteCAD
from ..cad.projection.collins import ProjectionPolynomial
from ..formula import And, Atom, BoolConst, Formula, Not, Or
from ..simplify.result import simplify_qe_formula


@dataclass(frozen=True)
class ExceptionalCause:
    """A polynomial whose zero set can carry nongeneric CAD behavior.

    The source string is intentionally small and stable so callers and tests can
    distinguish input boundaries from projection-generated degeneracy loci.
    """

    polynomial: sp.Poly
    source: str
    variables: tuple[sp.Symbol, ...]
    parents: tuple[str, ...] = ()

    def equation(self) -> sp.Equality:
        return sp.Eq(sp.expand(self.polynomial.as_expr()), 0)


def _normalize_poly(expr: sp.Expr, variables: Sequence[sp.Symbol]) -> sp.Poly | None:
    if not variables:
        return None
    try:
        poly = sp.Poly(sp.expand(expr), *variables)
    except Exception:
        return None
    if poly.is_zero or poly.total_degree() == 0:
        return None
    primitive = poly.primitive()[1]
    if primitive.LC().could_extract_minus_sign():
        primitive = -primitive
    return primitive


def _atom_boundary_exprs(formula: Formula) -> tuple[sp.Expr, ...]:
    if isinstance(formula, Atom):
        if formula.op in {"=", "<", "<=", ">", ">=", "!="}:
            return (sp.expand(formula.expr),)
        return tuple()
    if isinstance(formula, BoolConst):
        return tuple()
    if isinstance(formula, (And, Or)):
        out: list[sp.Expr] = []
        for arg in formula.args:
            out.extend(_atom_boundary_exprs(arg))
        return tuple(out)
    if isinstance(formula, Not):
        return _atom_boundary_exprs(formula.arg)
    return tuple()


def input_boundary_causes(
    formula: Formula, variables: Sequence[sp.Symbol]
) -> tuple[ExceptionalCause, ...]:
    """Return atom-boundary polynomials from the user formula.

    These are the most important exceptional
    equations: inequalities contribute their boundary, equalities contribute the
    algebraic variety, and disequalities contribute their deleted locus.
    """

    seen: set[str] = set()
    causes: list[ExceptionalCause] = []
    vars_tuple = tuple(variables)
    for expr in _atom_boundary_exprs(formula):
        poly = _normalize_poly(expr, vars_tuple)
        if poly is None:
            continue
        key = sp.sstr(sp.expand(poly.as_expr()))
        if key in seen:
            continue
        seen.add(key)
        causes.append(ExceptionalCause(poly, "input_boundary", vars_tuple))
    return tuple(causes)


def projection_causes(
    cad: CompleteCAD, variables: Sequence[sp.Symbol]
) -> tuple[ExceptionalCause, ...]:
    """Return discriminant/resultant/nullification-style projection causes."""

    vars_tuple = tuple(variables)
    seen: set[tuple[str, str]] = set()
    causes: list[ExceptionalCause] = []
    for level in cad.tower.levels:
        for entry in level.entries:
            source = _projection_source(entry)
            if source is None:
                continue
            poly = _normalize_poly(entry.poly.as_expr(), vars_tuple[: level.level])
            if poly is None:
                continue
            key = (source, sp.sstr(sp.expand(poly.as_expr())))
            if key in seen:
                continue
            seen.add(key)
            causes.append(
                ExceptionalCause(
                    polynomial=poly,
                    source=source,
                    variables=tuple(vars_tuple[: level.level]),
                    parents=entry.parents,
                )
            )
    return tuple(causes)


def _projection_source(entry: ProjectionPolynomial) -> str | None:
    source = entry.source
    if "discriminant" in source:
        return "discriminant"
    if "resultant" in source:
        return "resultant"
    if "coefficient" in source or "content" in source:
        return "nullification"
    return None


def exceptional_formula_from_causes(causes: Iterable[ExceptionalCause]) -> sp.Expr:
    equations = [cause.equation() for cause in causes]
    if not equations:
        return sp.false
    return simplify_qe_formula(sp.Or(*equations), implication_minimize=False)


def relevant_causes_for_cells(
    causes: Iterable[ExceptionalCause],
    exceptional_formula: sp.Expr,
) -> tuple[ExceptionalCause, ...]:
    """Filter causes to those plausibly useful for the visible exception set.

    A full semantic minimal hitting set would require an additional QE pass. We
    keep the predicate intentionally conservative: when the exceptional cell
    formula is nonempty, retain the atom-boundary causes first, and let callers
    include projection causes in diagnostics for auditability.
    """

    if exceptional_formula is sp.false or exceptional_formula == sp.false:
        return tuple()
    kept = [cause for cause in causes if cause.source == "input_boundary"]
    return tuple(kept or tuple(causes))


__all__ = [
    "ExceptionalCause",
    "exceptional_formula_from_causes",
    "input_boundary_causes",
    "projection_causes",
    "relevant_causes_for_cells",
]
