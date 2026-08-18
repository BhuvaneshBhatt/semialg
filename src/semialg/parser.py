from __future__ import annotations

import re

import sympy as sp
from sympy.parsing.sympy_parser import (
    convert_equals_signs,
    parse_expr,
    standard_transformations,
)
from sympy.parsing.sympy_parser import (
    implicit_multiplication_application as implicit_mul,
)

from .formula import And, Atom, BoolConst, Formula, Not, Or, ParsedPrenexFormula

Quantifier = tuple[str, sp.Symbol]
_TRANSFORMS = standard_transformations + (implicit_mul, convert_equals_signs)
_REL_RE = re.compile(r"(<=|>=|!=|=|<|>)")


def parse_quantified_formula(
    text: str, symbols: dict[str, sp.Symbol] | None = None, variable_order=None
):
    local_symbols = dict(symbols or {})
    quantifiers, matrix_text = _split_quantifier_prefix(text, local_symbols)
    matrix = _parse_form_order(matrix_text, local_symbols)
    matrix_expr = _formula_to_sympy(matrix)
    matrix_symbols = tuple(sorted(matrix_expr.free_symbols, key=lambda s: s.name))
    quantified_vars = tuple(var for _, var in quantifiers)
    free_vars = tuple(sym for sym in matrix_symbols if sym not in quantified_vars)
    if variable_order is not None:
        vars_ = tuple(variable_order)
    else:
        vars_ = tuple(list(free_vars) + [v for v in quantified_vars if v not in free_vars])
    return ParsedPrenexFormula(
        vars=vars_, quantifiers=tuple(quantifiers), matrix=matrix, matrix_expr=matrix_expr
    )


def _split_quantifier_prefix(
    text: str, symbols: dict[str, sp.Symbol]
) -> tuple[list[Quantifier], str]:
    s = text.strip()
    quantifiers: list[Quantifier] = []
    while True:
        m = re.match(r"^(exists|forall)\b", s, flags=re.IGNORECASE)
        if not m:
            break
        qname = m.group(1).lower()
        rest = s[m.end() :].lstrip()
        dot_idx = _find_top_level_char(rest, ".")
        if dot_idx < 0:
            raise ValueError("Expected '.' after quantified variable list")
        names_part = rest[:dot_idx].strip()
        s = rest[dot_idx + 1 :].strip()
        names = [piece.strip() for piece in names_part.split(",") if piece.strip()]
        for name in names:
            if name not in symbols:
                symbols[name] = sp.Symbol(name, real=True)
            quantifiers.append((qname, symbols[name]))
    return quantifiers, s


def _parse_form_order(text: str, symbols: dict[str, sp.Symbol]) -> Formula:
    s = text.strip()
    while s.startswith("(") and s.endswith(")") and _balanced_parens(s[1:-1]):
        s = s[1:-1].strip()
    # or
    parts = _split_top_level_regex(s, re.compile(r"\b(?:or)\b|\|", re.IGNORECASE))
    if len(parts) > 1:
        return Or(tuple(_parse_form_order(part, symbols) for part in parts))
    parts = _split_top_level_regex(s, re.compile(r"\b(?:and)\b|&", re.IGNORECASE))
    if len(parts) > 1:
        return And(tuple(_parse_form_order(part, symbols) for part in parts))
    if re.match(r"^(?:not\b|~)", s, flags=re.IGNORECASE):
        rest = re.sub(r"^(?:not\b|~)", "", s, count=1, flags=re.IGNORECASE).strip()
        return Not(_parse_form_order(rest, symbols))
    for op in ("<=", ">=", "!=", "=", "<", ">"):
        idx = _find_top_level_op(s, op)
        if idx >= 0:
            left = s[:idx].strip()
            right = s[idx + len(op) :].strip()
            lexpr = _parse_expr(left, symbols)
            rexpr = _parse_expr(right, symbols)
            return Atom(sp.expand(lexpr - rexpr), op)
    if s.lower() == "true":
        return BoolConst(True)
    if s.lower() == "false":
        return BoolConst(False)
    raise ValueError(f"Could not parse formula text: {text!r}")


def _parse_expr(text: str, symbols: dict[str, sp.Symbol]) -> sp.Expr:
    normalized = text.replace("^", "**")
    return parse_expr(normalized, local_dict=symbols, transformations=_TRANSFORMS, evaluate=False)


def _formula_to_sympy(formula: Formula) -> sp.Expr:
    if isinstance(formula, BoolConst):
        return sp.true if formula.value else sp.false
    if isinstance(formula, Atom):
        e = formula.expr
        if formula.op == "=":
            return sp.Eq(e, 0)
        if formula.op == "!=":
            return sp.Ne(e, 0)
        if formula.op == "<":
            return e < 0
        if formula.op == "<=":
            return e <= 0
        if formula.op == ">":
            return e > 0
        if formula.op == ">=":
            return e >= 0
    if isinstance(formula, And):
        return sp.And(*(_formula_to_sympy(arg) for arg in formula.args))
    if isinstance(formula, Or):
        return sp.Or(*(_formula_to_sympy(arg) for arg in formula.args))
    if isinstance(formula, Not):
        return sp.Not(_formula_to_sympy(formula.arg))
    raise TypeError(type(formula))


def _find_top_level_char(text: str, target: str) -> int:
    depth = 0
    for i, ch in enumerate(text):
        if ch == "(":
            depth += 1
        elif ch == ")":
            depth -= 1
        elif ch == target and depth == 0:
            return i
    return -1


def _split_top_level_regex(text: str, pattern: re.Pattern) -> list[str]:
    depth = 0
    start = 0
    parts: list[str] = []
    i = 0
    while i < len(text):
        ch = text[i]
        if ch == "(":
            depth += 1
            i += 1
            continue
        if ch == ")":
            depth -= 1
            i += 1
            continue
        if depth == 0:
            m = pattern.match(text, i)
            if m:
                parts.append(text[start:i].strip())
                i = m.end()
                start = i
                continue
        i += 1
    if start == 0:
        return [text]
    parts.append(text[start:].strip())
    return [p for p in parts if p]


def _find_top_level_op(text: str, op: str) -> int:
    depth = 0
    i = 0
    while i <= len(text) - len(op):
        ch = text[i]
        if ch == "(":
            depth += 1
        elif ch == ")":
            depth -= 1
        elif depth == 0 and text.startswith(op, i):
            return i
        i += 1
    return -1


def _balanced_parens(text: str) -> bool:
    depth = 0
    for ch in text:
        if ch == "(":
            depth += 1
        elif ch == ")":
            depth -= 1
            if depth < 0:
                return False
    return depth == 0
