from __future__ import annotations

import sympy as sp

from ..formula import (
    And,
    Atom,
    BoolConst,
    Formula,
    Not,
    Or,
    ParsedPrenexFormula,
    to_sympy,
)


def _sort_key(formula: Formula) -> tuple:
    if isinstance(formula, BoolConst):
        return (0, int(formula.value))
    if isinstance(formula, Atom):
        return (1, formula.op, sp.srepr(sp.expand(formula.expr)))
    if isinstance(formula, Not):
        return (2, _sort_key(formula.arg))
    if isinstance(formula, And):
        return (3, tuple(_sort_key(a) for a in formula.args))
    if isinstance(formula, Or):
        return (4, tuple(_sort_key(a) for a in formula.args))
    return (9, repr(formula))


def normalize_formula(formula: Formula) -> Formula:
    if isinstance(formula, (BoolConst, Atom)):
        return formula
    if isinstance(formula, Not):
        inner = normalize_formula(formula.arg)
        if isinstance(inner, BoolConst):
            return BoolConst(not inner.value)
        return Not(inner)
    if isinstance(formula, And):
        args = []
        for arg in formula.args:
            norm = normalize_formula(arg)
            if isinstance(norm, BoolConst):
                if not norm.value:
                    return BoolConst(False)
                continue
            if isinstance(norm, And):
                args.extend(norm.args)
            else:
                args.append(norm)
        if not args:
            return BoolConst(True)
        dedup = sorted({_sort_key(a): a for a in args}.values(), key=_sort_key)
        if len(dedup) == 1:
            return dedup[0]
        return And(tuple(dedup))
    if isinstance(formula, Or):
        args = []
        for arg in formula.args:
            norm = normalize_formula(arg)
            if isinstance(norm, BoolConst):
                if norm.value:
                    return BoolConst(True)
                continue
            if isinstance(norm, Or):
                args.extend(norm.args)
            else:
                args.append(norm)
        if not args:
            return BoolConst(False)
        dedup = sorted({_sort_key(a): a for a in args}.values(), key=_sort_key)
        if len(dedup) == 1:
            return dedup[0]
        return Or(tuple(dedup))
    raise TypeError(f"Unsupported formula node: {type(formula)!r}")


def normalize_parsed_formula(parsed: ParsedPrenexFormula) -> ParsedPrenexFormula:
    matrix = normalize_formula(parsed.matrix)
    return ParsedPrenexFormula(
        vars=parsed.vars,
        quantifiers=parsed.quantifiers,
        matrix=matrix,
        matrix_expr=to_sympy(matrix),
    )
