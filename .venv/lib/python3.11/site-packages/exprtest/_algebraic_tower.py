"""Bounded triangular algebraic towers and sparse quotient reduction."""

from __future__ import annotations

from dataclasses import dataclass

import sympy as sp

from . import _config as cfg
from ._algebraic_model import _mod_poly, classify_algebraic_expression
from ._cost import stage_allowed
from ._errors import EXACT_METHOD_ERRORS
from ._fast import quick_reduce
from ._memo import minpoly_for
from ._result import Verdict, ZeroClassification


@dataclass(frozen=True)
class TowerStep:
    """One generator and defining relation in an algebraic tower."""

    value: sp.Expr
    symbol: sp.Symbol
    relation: sp.Poly


@dataclass(frozen=True)
class AlgebraicTower:
    """Triangular algebraic presentation of a closed exact expression."""

    expression: sp.Expr
    steps: tuple[TowerStep, ...]
    numerator: sp.Expr
    denominator: sp.Expr
    numerator_support: frozenset[sp.Symbol] = frozenset()
    denominator_support: frozenset[sp.Symbol] = frozenset()


def _tower_nodes(term: sp.Expr, out: list[sp.Expr]) -> None:
    """Collect algebraic generators after their dependencies."""
    if term.is_Rational:
        return
    for arg in term.args:
        _tower_nodes(arg, out)
    is_root = getattr(term, "is_CRootOf", False)
    is_radical = term.is_Pow and term.exp.is_Rational and not term.exp.is_Integer
    is_algnum = isinstance(term, sp.AlgebraicNumber)
    if (is_root or is_radical or is_algnum or term == sp.I) and term not in out:
        out.append(term)


def _tower_relation(term: sp.Expr, sym: sp.Symbol,
                    repl: dict[sp.Expr, sp.Symbol]) -> sp.Poly | None:
    """Build a defining relation whose coefficients use earlier symbols."""
    try:
        if term == sp.I:
            return sp.Poly(sym**2 + 1, sym, *repl.values(), domain=sp.QQ)
        if term.is_Pow and term.exp.is_Rational and not term.exp.is_Integer:
            power = sp.Rational(term.exp)
            num = int(power.p)
            den = int(power.q)
            base = term.base.xreplace(repl)
            rel = (sym**den - base**num if num >= 0
                   else base**(-num) * sym**den - 1)
            vars_ = (sym,) + tuple(repl.values())
            return sp.Poly(rel, *vars_, domain=sp.QQ)
        probe = sp.Dummy("z")
        rel = minpoly_for(term, probe).as_expr().xreplace({probe: sym})
        return sp.Poly(rel, sym, domain=sp.QQ)
    except EXACT_METHOD_ERRORS:
        return None


def _tower_parts(term: sp.Expr, syms: frozenset[sp.Symbol],
                 max_terms: int) -> tuple[sp.Expr, sp.Expr, frozenset[sp.Symbol]] | None:
    """Represent a rational DAG while tracking each subtree's tower support.

    Term-count checks use only variables that actually occur in the current
    subtree. This avoids repeatedly constructing multivariate ``Poly`` objects
    over unrelated tower levels.
    """
    if not syms or not term.has(*syms):
        return term, sp.Integer(1), frozenset()
    if term.is_Symbol and term in syms:
        return term, sp.Integer(1), frozenset((term,))

    def admissible(num: sp.Expr, den: sp.Expr, active: frozenset[sp.Symbol]) -> bool:
        if not active:
            return True
        try:
            vars_ = tuple(sorted(active, key=str))
            return (len(sp.Poly(num, *vars_, domain=sp.QQ).terms()) <= max_terms
                    and len(sp.Poly(den, *vars_, domain=sp.QQ).terms()) <= max_terms)
        except EXACT_METHOD_ERRORS:
            return False

    if term.is_Add:
        num = sp.Integer(0)
        den = sp.Integer(1)
        active: frozenset[sp.Symbol] = frozenset()
        for arg in term.args:
            part = _tower_parts(arg, syms, max_terms)
            if part is None:
                return None
            pn, pd, support = part
            num = num * pd + pn * den
            den *= pd
            active = active | support
            if not admissible(num, den, active):
                return None
        return num, den, active
    if term.is_Mul:
        num = sp.Integer(1)
        den = sp.Integer(1)
        active: frozenset[sp.Symbol] = frozenset()
        for arg in term.args:
            part = _tower_parts(arg, syms, max_terms)
            if part is None:
                return None
            pn, pd, support = part
            num *= pn
            den *= pd
            active = active | support
            if not admissible(num, den, active):
                return None
        return num, den, active
    if term.is_Pow and term.exp.is_Integer:
        part = _tower_parts(term.base, syms, max_terms)
        if part is None:
            return None
        num, den, active = part
        power = int(term.exp)
        if power < 0:
            num, den = den, num
            power = -power
        if power.bit_length() > cfg.EXACT_MAX_POW_BITS:
            return None
        num, den = num**power, den**power
        return (num, den, active) if admissible(num, den, active) else None
    return None


def _tower_support(active: frozenset[sp.Symbol],
                   steps: tuple[TowerStep, ...]) -> frozenset[sp.Symbol]:
    """Close direct symbol support over dependencies in defining relations."""
    if not active:
        return active
    known = set(active)
    by_symbol = {step.symbol: step for step in steps}
    changed = True
    while changed:
        changed = False
        for symbol in tuple(known):
            step = by_symbol.get(symbol)
            if step is None:
                continue
            deps = (step.relation.as_expr().free_symbols - {symbol}) & set(by_symbol)
            for dep in deps:
                if dep not in known:
                    known.add(dep)
                    changed = True
    return frozenset(known)


def build_algebraic_tower(term: sp.Expr) -> AlgebraicTower | None:
    """Build a dependency-ordered triangular model for nested algebraic values.

    The model is constructed structurally rather than through ``cancel`` or
    ``together``, keeping this representation suitable for the fast oracle.
    """
    term = sp.sympify(term)
    if not stage_allowed(term, "tower"):
        return None
    info = classify_algebraic_expression(term)
    if not info.is_algebraic:
        return None
    nodes: list[sp.Expr] = []
    _tower_nodes(term, nodes)
    if len(nodes) > min(cfg.ALG_TOWER_MAX_GENS,
                        cfg.TOWER_MAX_GENS):
        return None
    repl: dict[sp.Expr, sp.Symbol] = {}
    steps: list[TowerStep] = []
    for pos, value in enumerate(nodes):
        sym = sp.Dummy(f"t{pos}")
        relation = _tower_relation(value, sym, repl)
        if relation is None:
            return None
        repl[value] = sym
        steps.append(TowerStep(value, sym, relation))
    try:
        represented = term.xreplace(repl)
        syms = frozenset(step.symbol for step in steps)
        parts = _tower_parts(represented, syms, cfg.EXACT_MAX_POLY_TERMS)
        if parts is None:
            return None
        top, bottom, support = parts
        step_tuple = tuple(steps)
        top_direct = frozenset(top.free_symbols) & syms
        bottom_direct = frozenset(bottom.free_symbols) & syms
        top_support = _tower_support(top_direct or support, step_tuple)
        bottom_support = _tower_support(bottom_direct, step_tuple)
        if syms:
            if top_support:
                sp.Poly(top, *tuple(sorted(top_support, key=str)), domain=sp.QQ)
            if bottom_support:
                sp.Poly(bottom, *tuple(sorted(bottom_support, key=str)), domain=sp.QQ)
        return AlgebraicTower(term, step_tuple, top, bottom, top_support, bottom_support)
    except EXACT_METHOD_ERRORS:
        return None



def _sparse_poly(expr: sp.Expr, variables: tuple[sp.Symbol, ...]):
    """Convert a bounded rational polynomial to a sparse monomial mapping."""
    try:
        poly = sp.Poly(expr, *variables, domain=sp.QQ)
        terms = {tuple(mon): sp.Rational(coeff) for mon, coeff in poly.terms() if coeff != 0}
        if len(terms) > cfg.TOWER_SPARSE_MAX_TERMS:
            return None
        return terms
    except EXACT_METHOD_ERRORS:
        return None


def _sparse_expr(terms: dict[tuple[int, ...], sp.Rational],
                 variables: tuple[sp.Symbol, ...]) -> sp.Expr:
    pieces = []
    for powers, coeff in terms.items():
        item = sp.sympify(coeff)
        for variable, power in zip(variables, powers):
            if power:
                item *= variable ** power
        pieces.append(item)
    return sp.Add(*pieces) if pieces else sp.Integer(0)


def _sparse_reduce_var(terms: dict[tuple[int, ...], sp.Rational],
                       relation: dict[tuple[int, ...], sp.Rational],
                       index: int):
    """Reduce one monic triangular relation directly in sparse quotient form."""
    if not terms:
        return terms
    degree = max((powers[index] for powers in relation), default=0)
    lead_power = tuple(degree if pos == index else 0 for pos in range(len(next(iter(terms)))))
    if degree <= 0 or relation.get(lead_power) != 1:
        return None
    lower = [(powers, coeff) for powers, coeff in relation.items() if powers != lead_power]
    out = dict(terms)
    while True:
        target = next((powers for powers in out if powers[index] >= degree), None)
        if target is None:
            break
        coeff = out.pop(target)
        base = list(target)
        base[index] -= degree
        for rel_power, rel_coeff in lower:
            powers = tuple(a + b for a, b in zip(base, rel_power))
            value = out.get(powers, sp.Rational(0)) - coeff * rel_coeff
            if value:
                out[powers] = value
            else:
                out.pop(powers, None)
        if len(out) > cfg.TOWER_SPARSE_MAX_TERMS:
            return None
    return out


def _tower_sparse_rem(poly: sp.Expr, tower: AlgebraicTower,
                      support: frozenset[sp.Symbol]) -> sp.Expr | None:
    """Reduce a tower polynomial with sparse monomial arithmetic.

    Only monic triangular relations are admitted. Unsupported relations fall
    back to the existing bounded SymPy polynomial reducer.
    """
    active_steps = tuple(step for step in tower.steps if step.symbol in support)
    if not active_steps:
        return poly
    variables = tuple(step.symbol for step in active_steps)
    terms = _sparse_poly(poly, variables)
    if terms is None:
        return None
    for index in range(len(active_steps) - 1, -1, -1):
        relation = _sparse_poly(active_steps[index].relation.as_expr(), variables)
        if relation is None:
            return None
        terms = _sparse_reduce_var(terms, relation, index)
        if terms is None:
            return None
        if not terms:
            return sp.Integer(0)
    return quick_reduce(_sparse_expr(terms, variables))


def _tower_rem(poly: sp.Expr, tower: AlgebraicTower,
               support: frozenset[sp.Symbol] | None = None) -> sp.Expr:
    """Reduce only the dependency-closed tower levels used by ``poly``."""
    rem = poly
    active = support if support is not None else _tower_support(
        frozenset(rem.free_symbols), tower.steps)
    sparse = _tower_sparse_rem(rem, tower, active)
    if sparse is not None:
        return sparse
    for step in reversed(tower.steps):
        if step.symbol not in active or step.symbol not in rem.free_symbols:
            continue
        try:
            divisor = sp.Poly(step.relation.as_expr(), step.symbol, domain="EX")
            rem = _mod_poly(rem, step.symbol, divisor)
        except EXACT_METHOD_ERRORS:
            return quick_reduce(rem)
        if rem == 0:
            break
    return quick_reduce(rem)


def _restore_tower(term: sp.Expr, tower: AlgebraicTower) -> sp.Expr:
    repl = {step.symbol: step.value for step in tower.steps}
    return quick_reduce(term.xreplace(repl))


def _tower_nonzero(term: sp.Expr) -> bool:
    """Prove nonzero for a small closed algebraic value without heuristics."""
    term = quick_reduce(term)
    if term.is_zero is False:
        return True
    if term.is_Rational:
        return term != 0
    try:
        from ._exact_constants import algebraic_sign
        return algebraic_sign(term) in (-1, 1)
    except EXACT_METHOD_ERRORS:
        return False


def tower_algebraic_test(term: sp.Expr) -> ZeroClassification:
    """Use a bounded triangular tower for exact zero and nonzero proofs."""
    tower = build_algebraic_tower(term)
    if tower is None:
        return ZeroClassification(
            Verdict.UNKNOWN, "algebraic-tower",
            detail="could not build a bounded triangular algebraic model",
        )
    top = _tower_rem(tower.numerator, tower, tower.numerator_support)
    bottom = _tower_rem(tower.denominator, tower, tower.denominator_support)
    restored_den = _restore_tower(bottom, tower)
    if top == 0:
        if _tower_nonzero(restored_den):
            return ZeroClassification(
                Verdict.ZERO_PROVEN, "algebraic-tower",
                detail="triangular quotient reduction produced zero with a proven nonzero denominator",
                evidence="triangular-algebraic-reduction",
            )
        return ZeroClassification(
            Verdict.UNKNOWN, "algebraic-tower",
            detail="zero numerator found but denominator nonvanishing is unresolved",
        )

    # A nonzero quotient remainder alone is not a proof because the triangular
    # relations may describe extra branches. Restore the selected algebraic
    # values and require an independent exact nonzero/sign proof.
    if not _tower_nonzero(restored_den):
        return ZeroClassification(
            Verdict.UNKNOWN, "algebraic-tower",
            detail="tower denominator nonvanishing is unresolved",
        )
    restored_top = _restore_tower(top, tower)
    if _tower_nonzero(restored_top):
        return ZeroClassification(
            Verdict.NONZERO_PROVEN, "algebraic-tower",
            detail="reduced tower numerator has an independent exact nonzero proof",
            evidence="triangular-algebraic-nonzero",
        )
    return ZeroClassification(
        Verdict.UNKNOWN, "algebraic-tower",
        detail="tower remainder is nonzero but selected-branch nonvanishing is unresolved",
    )

