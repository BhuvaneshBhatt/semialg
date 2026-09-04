"""Exact algebraization helpers for function-domain/range computations."""

from __future__ import annotations

from dataclasses import dataclass
from math import lcm

import sympy as sp

from .function_graph import UnsupportedFunctionGraph, semialgebraic_function_graph
from .internal_symbols import fresh_real_dummy


@dataclass(frozen=True)
class BackSubstitutionResult:
    expression: sp.Expr
    constraints: sp.Expr
    variables: tuple[sp.Symbol, ...]
    substitutions: tuple[tuple[sp.Symbol, sp.Expr], ...] = ()


@dataclass(frozen=True)
class ExactAlgebraization:
    expression: sp.Expr
    constraints: sp.Expr
    variables: tuple[sp.Symbol, ...]
    original_variable: sp.Symbol
    kind: str
    substitutions: tuple[tuple[sp.Expr, sp.Expr], ...]
    side_conditions: tuple[sp.Expr, ...]
    diagnostics: tuple[str, ...] = ()


def _contains_transcendental(expr: sp.Expr) -> bool:
    heads = {
        sp.sin,
        sp.cos,
        sp.tan,
        sp.cot,
        sp.sec,
        sp.csc,
        sp.exp,
        sp.sinh,
        sp.cosh,
        sp.tanh,
        sp.coth,
        sp.sech,
        sp.csch,
        sp.log,
    }
    return any(item.func in heads for item in sp.preorder_traversal(sp.sympify(expr)))


def _safe_algebraic_rhs(expr: sp.Expr) -> bool:
    if _contains_transcendental(expr):
        return False
    probe = fresh_real_dummy("semialg_backsub_probe")
    try:
        semialgebraic_function_graph(expr, probe)
    except (UnsupportedFunctionGraph, TypeError, ValueError, NotImplementedError):
        return False
    return True


def algebraic_back_substitute(
    expression: sp.Expr,
    constraints: sp.Expr,
    variables: tuple[sp.Symbol, ...],
) -> BackSubstitutionResult:
    """Substitute uniquely solved algebraic equalities before harder analysis.

    Only conjunctions are rewritten.  An equality is used when it uniquely
    solves one listed variable as a supported algebraic/semialgebraic expression
    not containing that variable.  This preserves equivalence while reducing
    the number of variables seen by later transcendental algebraizers.
    """

    expr = sp.sympify(expression)
    condition = sp.sympify(constraints)
    vars_ = list(variables)
    substitutions: list[tuple[sp.Symbol, sp.Expr]] = []
    atoms = list(condition.args) if isinstance(condition, sp.And) else [condition]

    changed = True
    while changed:
        changed = False
        for atom in tuple(atoms):
            if not isinstance(atom, sp.Equality):
                continue
            current = sp.simplify(atom.subs(dict(substitutions)))
            candidates = [v for v in vars_ if current.has(v)]
            # Prefer eliminating a variable that occurs in the function being
            # analyzed.  This turns y=g(x), f(y) into f(g(x)) rather than
            # unnecessarily solving for the independent x variable.
            candidates.sort(key=lambda v: (not expr.has(v), vars_.index(v)))
            for variable in candidates:
                try:
                    solutions = sp.solve(current, variable, dict=False)
                except (NotImplementedError, TypeError, ValueError):
                    continue
                if len(solutions) != 1:
                    continue
                rhs = sp.simplify(solutions[0])
                if variable in rhs.free_symbols or not _safe_algebraic_rhs(rhs):
                    continue
                expr = sp.simplify(expr.subs(variable, rhs))
                atoms = [sp.simplify(a.subs(variable, rhs)) for a in atoms if a != atom]
                substitutions = [
                    (v, sp.simplify(value.subs(variable, rhs))) for v, value in substitutions
                ]
                substitutions.append((variable, rhs))
                vars_.remove(variable)
                changed = True
                break
            if changed:
                break

    cleaned = [a for a in atoms if a is not sp.true and a != sp.true]
    condition = sp.And(*cleaned) if cleaned else sp.true
    return BackSubstitutionResult(expr, condition, tuple(vars_), tuple(substitutions))


def _linear_rate(argument: sp.Expr, variable: sp.Symbol) -> sp.Expr | None:
    try:
        poly = sp.Poly(sp.expand(argument), variable)
    except (sp.PolynomialError, TypeError, ValueError):
        return None
    if poly.degree() != 1 or sp.simplify(poly.eval(0)) != 0:
        return None
    coeff = sp.simplify(poly.coeff_monomial(variable))
    if coeff == 0 or coeff.is_real is not True:
        return None
    return coeff


def _commensurate_base(rates: list[sp.Expr]) -> tuple[sp.Expr, tuple[int, ...]] | None:
    if not rates:
        return None
    first = rates[0]
    ratios: list[sp.Rational] = []
    for rate in rates:
        ratio = sp.simplify(rate / first)
        if not isinstance(ratio, sp.Rational):
            return None
        ratios.append(ratio)
    scale = lcm(*(int(r.q) for r in ratios))
    base = sp.simplify(first / scale)
    multiples = tuple(int(r * scale) for r in ratios)
    if base == 0 or base.is_real is not True:
        return None
    return base, multiples


def _trig_multiple(kind, n: int, s: sp.Symbol, c: sp.Symbol) -> sp.Expr:
    theta = sp.Dummy("theta", real=True)
    expanded = sp.expand_trig(kind(sp.Integer(n) * theta))
    return sp.expand(expanded.xreplace({sp.sin(theta): s, sp.cos(theta): c}))


def _algebraize_trigonometric(
    expression: sp.Expr, constraints: sp.Expr, variable: sp.Symbol
) -> ExactAlgebraization | None:
    combined = sp.Tuple(expression, constraints)
    atoms = sorted(combined.atoms(sp.sin, sp.cos), key=sp.default_sort_key)
    if not atoms:
        return None
    # Exact polynomial/rational trig layer intentionally excludes other
    # transcendental heads and tangent-family poles in this transformation.
    if combined.has(sp.exp, sp.sinh, sp.cosh, sp.tanh, sp.log, sp.tan, sp.cot, sp.sec, sp.csc):
        return None
    rates = [_linear_rate(atom.args[0], variable) for atom in atoms]
    if any(rate is None for rate in rates):
        return None
    commensurate = _commensurate_base([rate for rate in rates if rate is not None])
    if commensurate is None:
        return None
    base, multiples = commensurate
    s = fresh_real_dummy("semialg_sin")
    c = fresh_real_dummy("semialg_cos")
    replacements: list[tuple[sp.Expr, sp.Expr]] = []
    for atom, multiple in zip(atoms, multiples, strict=True):
        replacements.append((atom, _trig_multiple(atom.func, multiple, s, c)))
    mapping = dict(replacements)
    new_expr = sp.cancel(sp.expand(expression.xreplace(mapping)))
    new_condition = sp.simplify(constraints.xreplace(mapping))
    if new_expr.has(variable) or new_condition.has(variable):
        return None
    if any(
        item.func in {sp.sin, sp.cos}
        for item in sp.preorder_traversal(sp.Tuple(new_expr, new_condition))
    ):
        return None
    circle = sp.Eq(s**2 + c**2, 1)
    return ExactAlgebraization(
        new_expr,
        sp.And(new_condition, circle),
        (s, c),
        variable,
        "commensurate_trigonometric",
        tuple(replacements),
        (circle,),
        (f"fundamental_frequency={sp.sstr(base)}",),
    )


def _algebraize_exp_hyperbolic(
    expression: sp.Expr, constraints: sp.Expr, variable: sp.Symbol
) -> ExactAlgebraization | None:
    combined = sp.Tuple(expression, constraints)
    atoms = sorted(
        set(combined.atoms(sp.exp)) | set(combined.atoms(sp.sinh)) | set(combined.atoms(sp.cosh)),
        key=sp.default_sort_key,
    )
    if not atoms:
        return None
    if combined.has(
        sp.sin, sp.cos, sp.tan, sp.cot, sp.sec, sp.csc, sp.log, sp.tanh, sp.coth, sp.sech, sp.csch
    ):
        return None
    rates = [_linear_rate(atom.args[0], variable) for atom in atoms]
    if any(rate is None for rate in rates):
        return None
    commensurate = _commensurate_base([rate for rate in rates if rate is not None])
    if commensurate is None:
        return None
    base, multiples = commensurate
    t = fresh_real_dummy("semialg_exp")
    replacements: list[tuple[sp.Expr, sp.Expr]] = []
    for atom, n in zip(atoms, multiples, strict=True):
        power = t**n
        if atom.func is sp.exp:
            value = power
        elif atom.func is sp.sinh:
            value = sp.Rational(1, 2) * (power - t ** (-n))
        else:
            value = sp.Rational(1, 2) * (power + t ** (-n))
        replacements.append((atom, value))
    mapping = dict(replacements)
    new_expr = sp.cancel(expression.xreplace(mapping))
    new_condition = sp.simplify(constraints.xreplace(mapping))
    if new_expr.has(variable) or new_condition.has(variable):
        return None
    if _contains_transcendental(sp.Tuple(new_expr, new_condition)):
        return None
    positive = t > 0
    return ExactAlgebraization(
        new_expr,
        sp.And(new_condition, positive),
        (t,),
        variable,
        "commensurate_exponential_hyperbolic",
        tuple(replacements),
        (positive,),
        (f"fundamental_rate={sp.sstr(base)}",),
    )


def exact_algebraize_function_problem(
    expression: sp.Expr,
    constraints: sp.Expr,
    variables: tuple[sp.Symbol, ...],
) -> ExactAlgebraization | None:
    """Return an exact finite semialgebraic reduction when one is recognized."""

    if len(variables) != 1:
        return None
    variable = variables[0]
    for transform in (_algebraize_trigonometric, _algebraize_exp_hyperbolic):
        result = transform(expression, constraints, variable)
        if result is not None:
            return result
    return None


__all__ = [
    "BackSubstitutionResult",
    "ExactAlgebraization",
    "algebraic_back_substitute",
    "exact_algebraize_function_problem",
]
