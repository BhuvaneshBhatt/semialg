from __future__ import annotations

import re
from collections.abc import Sequence
from dataclasses import dataclass

import sympy as sp
from sympy.core.relational import (
    Equality,
    GreaterThan,
    LessThan,
    StrictGreaterThan,
    StrictLessThan,
    Unequality,
)
from sympy.logic.boolalg import (
    And as SymAnd,
)
from sympy.logic.boolalg import (
    BooleanFalse,
    BooleanTrue,
)
from sympy.logic.boolalg import (
    Implies as SymImplies,
)
from sympy.logic.boolalg import (
    Not as SymNot,
)
from sympy.logic.boolalg import (
    Or as SymOr,
)
from sympy.parsing.sympy_parser import (
    convert_equals_signs,
    parse_expr,
    standard_transformations,
)
from sympy.parsing.sympy_parser import (
    implicit_multiplication_application as implicit_mul,
)


class FormulaOps:
    def __and__(self, other):
        return And((self, other))

    def __or__(self, other):
        return Or((self, other))

    def __invert__(self):
        return Not(self)


@dataclass(frozen=True)
class Atom(FormulaOps):
    expr: sp.Expr
    op: str


@dataclass(frozen=True)
class BoolConst(FormulaOps):
    value: bool


@dataclass(frozen=True)
class And(FormulaOps):
    args: tuple[Formula, ...]


@dataclass(frozen=True)
class Or(FormulaOps):
    args: tuple[Formula, ...]


@dataclass(frozen=True)
class Not(FormulaOps):
    arg: Formula


Formula = Atom | BoolConst | And | Or | Not
Quantifier = tuple[str, sp.Symbol]


@dataclass(frozen=True)
class ParsedPrenexFormula:
    vars: tuple[sp.Symbol, ...]
    quantifiers: tuple[Quantifier, ...]
    matrix: Formula
    matrix_expr: sp.Expr

    @property
    def quantified_expr(self) -> sp.Expr:
        """Return this parsed prefix as semialg ``Exists``/``ForAll`` nodes."""

        from .quantifiers import apply_quantifiers

        return apply_quantifiers(self.matrix_expr, self.quantifiers)


_TRANSFORMS = standard_transformations + (implicit_mul, convert_equals_signs)


def parse_formula(expr: sp.Expr) -> Formula:
    if expr is True or isinstance(expr, BooleanTrue):
        return BoolConst(True)
    if expr is False or isinstance(expr, BooleanFalse):
        return BoolConst(False)
    if isinstance(expr, Equality):
        return Atom(sp.expand(expr.lhs - expr.rhs), "=")
    if isinstance(expr, Unequality):
        return Atom(sp.expand(expr.lhs - expr.rhs), "!=")
    if isinstance(expr, StrictLessThan):
        return Atom(sp.expand(expr.lhs - expr.rhs), "<")
    if isinstance(expr, LessThan):
        return Atom(sp.expand(expr.lhs - expr.rhs), "<=")
    if isinstance(expr, StrictGreaterThan):
        return Atom(sp.expand(expr.lhs - expr.rhs), ">")
    if isinstance(expr, GreaterThan):
        return Atom(sp.expand(expr.lhs - expr.rhs), ">=")
    if isinstance(expr, SymAnd):
        return And(tuple(parse_formula(arg) for arg in expr.args))
    if isinstance(expr, SymOr):
        return Or(tuple(parse_formula(arg) for arg in expr.args))
    if isinstance(expr, SymImplies):
        left, right = expr.args
        return Or((Not(parse_formula(left)), parse_formula(right)))
    if isinstance(expr, SymNot):
        return Not(parse_formula(expr.args[0]))
    raise TypeError(f"Unsupported formula expression: {expr!r}")


def to_sympy(formula: Formula) -> sp.Expr:
    if isinstance(formula, BoolConst):
        return sp.true if formula.value else sp.false
    if isinstance(formula, Atom):
        e = sp.expand(formula.expr)
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
        raise ValueError(f"Unknown operator {formula.op}")
    if isinstance(formula, And):
        return sp.And(*(to_sympy(arg) for arg in formula.args))
    if isinstance(formula, Or):
        return sp.Or(*(to_sympy(arg) for arg in formula.args))
    if isinstance(formula, Not):
        return sp.Not(to_sympy(formula.arg))
    raise TypeError(f"Unsupported formula node: {type(formula)}")


def formula_polynomials(formula: Formula) -> list[sp.Expr]:
    if isinstance(formula, Atom):
        return [sp.expand(formula.expr)]
    if isinstance(formula, BoolConst):
        return []
    if isinstance(formula, (And, Or)):
        out: list[sp.Expr] = []
        for arg in formula.args:
            out.extend(formula_polynomials(arg))
        return out
    if isinstance(formula, Not):
        return formula_polynomials(formula.arg)
    raise TypeError(f"Unsupported formula node: {type(formula)}")


def equational_constraints(formula: Formula) -> list[sp.Expr]:
    if isinstance(formula, Atom):
        return [sp.expand(formula.expr)] if formula.op == "=" else []
    if isinstance(formula, BoolConst):
        return []
    if isinstance(formula, And):
        out: list[sp.Expr] = []
        for arg in formula.args:
            out.extend(equational_constraints(arg))
        return list(dict.fromkeys(out))
    if isinstance(formula, Or):
        common: set[sp.Expr] | None = None
        for arg in formula.args:
            current = set(equational_constraints(arg))
            common = current if common is None else common & current
        return list(common or set())
    if isinstance(formula, Not):
        return []
    raise TypeError(f"Unsupported formula node: {type(formula)}")


def parse_formula_text(
    text: str, symbols: dict[str, sp.Symbol] | None = None
) -> tuple[sp.Expr, Formula]:
    expr = _parse_matrix_text(text, symbols=symbols)
    return expr, parse_formula(expr)


def parse_quantified_expr(
    expr: sp.Expr,
    variable_order: Sequence[sp.Symbol] | None = None,
) -> ParsedPrenexFormula:
    """Parse a semialg ``Exists``/``ForAll`` prenex expression.

    This is the expression-level counterpart of :func:`parse_quant_form_text`.
    The quantifier nodes are lowered to semialg's internal prefix representation
    while preserving ``quantified_expr`` for round-tripping.
    """

    from .quantifiers import split_quantifiers

    quantifiers, matrix_expr = split_quantifiers(expr)
    matrix = parse_formula(matrix_expr)
    matrix_symbols = tuple(sorted(matrix_expr.free_symbols, key=lambda s: s.name))
    quantified_vars = tuple(var for _, var in quantifiers)
    free_vars = tuple(sym for sym in matrix_symbols if sym not in quantified_vars)

    if variable_order is not None:
        vars_ = tuple(variable_order)
    else:
        vars_ = tuple(dict.fromkeys(free_vars + quantified_vars))

    return ParsedPrenexFormula(
        vars=vars_, quantifiers=quantifiers, matrix=matrix, matrix_expr=matrix_expr
    )


def parse_quant_form_text(
    text: str,
    symbols: dict[str, sp.Symbol] | None = None,
    variable_order: Sequence[sp.Symbol] | None = None,
) -> ParsedPrenexFormula:
    local_symbols = dict(symbols or {})
    quantifiers, matrix_text = _split_quantifier_prefix(text, local_symbols)
    matrix_expr = _parse_matrix_text(matrix_text, symbols=local_symbols)
    matrix = parse_formula(matrix_expr)
    matrix_symbols = tuple(sorted(matrix_expr.free_symbols, key=lambda s: s.name))
    quantified_vars = tuple(var for _, var in quantifiers)
    free_vars = tuple(sym for sym in matrix_symbols if sym not in quantified_vars)

    if variable_order is not None:
        vars_ = tuple(variable_order)
    else:
        ordered = []
        seen = set()
        for sym in free_vars + quantified_vars:
            if sym not in seen:
                ordered.append(sym)
                seen.add(sym)
        vars_ = tuple(ordered)

    return ParsedPrenexFormula(
        vars=vars_, quantifiers=tuple(quantifiers), matrix=matrix, matrix_expr=matrix_expr
    )


def _split_quantifier_prefix(
    text: str, symbols: dict[str, sp.Symbol] | None
) -> tuple[list[Quantifier], str]:
    s = text.strip()
    quantifiers: list[Quantifier] = []
    while True:
        match = re.match(r"^(exists|forall)\b", s, flags=re.IGNORECASE)
        if not match:
            break
        qname = match.group(1).lower()
        rest = s[match.end() :].lstrip()
        dot_idx = _find_top_level_dot(rest)
        if dot_idx < 0:
            raise ValueError("Expected '.' after quantified variable list")
        names_part = rest[:dot_idx].strip()
        s = rest[dot_idx + 1 :].strip()
        if not names_part:
            raise ValueError("Expected quantified variable name(s)")
        names = [piece.strip() for piece in names_part.split(",") if piece.strip()]
        for name in names:
            if not re.fullmatch(r"[A-Za-z_]\w*", name):
                raise ValueError(f"Invalid quantified variable name: {name!r}")
            quantifiers.append((qname, _get_or_create_symbol(name, symbols)))
    return quantifiers, s


def _find_top_level_dot(text: str) -> int:
    depth = 0
    for i, ch in enumerate(text):
        if ch == "(":
            depth += 1
        elif ch == ")":
            depth -= 1
        elif ch == "." and depth == 0:
            return i
    return -1


def _parse_matrix_text(text: str, symbols: dict[str, sp.Symbol] | None = None) -> sp.Expr:
    local_dict = dict(symbols or {})
    implication = _split_top_level_keyword(text, "implies")
    if implication is not None:
        left, right = implication
        return sp.Implies(
            _parse_matrix_text(left, symbols=local_dict),
            _parse_matrix_text(right, symbols=local_dict),
            evaluate=False,
        )
    normalized = text.replace("^", "**")
    normalized = re.sub(r"\band\b", " & ", normalized, flags=re.IGNORECASE)
    normalized = re.sub(r"\bor\b", " | ", normalized, flags=re.IGNORECASE)
    normalized = re.sub(r"\bnot\b", " ~ ", normalized, flags=re.IGNORECASE)
    expr = parse_expr(
        normalized, local_dict=local_dict, transformations=_TRANSFORMS, evaluate=False
    )
    return expr


def _split_top_level_keyword(text: str, keyword: str) -> tuple[str, str] | None:
    pattern = re.compile(rf"\b{re.escape(keyword)}\b", flags=re.IGNORECASE)
    for match in pattern.finditer(text):
        depth = 0
        for ch in text[: match.start()]:
            if ch == "(":
                depth += 1
            elif ch == ")":
                depth -= 1
        if depth == 0:
            left = text[: match.start()].strip()
            right = text[match.end() :].strip()
            if left and right:
                return left, right
    return None


def _get_or_create_symbol(name: str, symbols: dict[str, sp.Symbol] | None) -> sp.Symbol:
    if symbols is not None and name in symbols:
        return symbols[name]
    sym = sp.Symbol(name, real=True)
    if symbols is not None:
        symbols[name] = sym
    return sym
