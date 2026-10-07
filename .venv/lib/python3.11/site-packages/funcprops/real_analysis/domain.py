from __future__ import annotations

import sympy as sp
from sympy import S

from .. import _semialg
from ..engine import evaluate_property, simplify_condition
from ..exceptions import SYMBOLIC_QUERY_ERRORS
from ..results import PropertyResult
from ..stratification import (
    semantic_condition,
    stratify_boolean_by_parameters,
)
from .common import is_semialgebraic_expression


def _pow_real_domain(expr, var):
    base, exponent = expr.args
    inner = _function_domain_impl(base, var)
    if exponent.is_Integer:
        if exponent.is_negative:
            return sp.And(inner, sp.Ne(base, 0))
        return inner
    if exponent.is_Rational:
        p, q = int(exponent.p), int(exponent.q)
        if q % 2 == 0:
            # Even root: nonnegative base. A negative exponent also excludes 0.
            rel = base > 0 if p < 0 else base >= 0
            return sp.And(inner, rel)
        # Odd root is real for all real bases, but negative powers exclude zero.
        return sp.And(inner, sp.Ne(base, 0)) if p < 0 else inner
    # Irrational real powers use the principal real branch only for positive base.
    if exponent.is_real is True:
        return sp.And(inner, base > 0)
    return inner


def _function_domain_impl(expr, var, *, domain=S.Reals, stratified=False):
    """Largest supported real/complex domain condition for ``expr`` in ``var``.

    funcprops owns transcendental structural rules. semialg is invoked only for
    expressions that are structurally semialgebraic.
    """
    expr = sp.sympify(expr)
    var = sp.sympify(var)
    real_mode = domain in (S.Reals, "real", "reals", "R", "RR")
    complex_mode = domain in (S.Complexes, "complex", "complexes", "C", "CC")
    if not (real_mode or complex_mode):
        raise ValueError("domain must be Reals or Complexes")
    if expr == var or not expr.has(var):
        return S.true

    if real_mode and isinstance(expr, sp.Piecewise):
        clauses = []
        covered = S.false
        for value, cond in expr.args:
            active = (
                sp.And(cond, sp.Not(covered)) if cond is not S.true else sp.Not(covered)
            )
            clauses.append(
                sp.And(active, _function_domain_impl(value, var, domain=S.Reals))
            )
            covered = sp.Or(covered, cond)
        return sp.simplify(sp.Or(*clauses)) if clauses else S.false

    if real_mode and expr.func is sp.Pow:
        if is_semialgebraic_expression(expr, var):
            try:
                return sp.sympify(_semialg.domain(expr, var))
            except SYMBOLIC_QUERY_ERRORS:
                pass
        # Symbolic exponents have genuinely parameter-dependent real domains
        # (integer, negative-integer, and noninteger cases).  The catalogue is
        # substantially more informative than the conservative x>0 fallback.
        if expr.exp.free_symbols:
            reg = evaluate_property(expr, "RealDomain", var=var, mode="real")
            if reg.primary is not None:
                dom = semantic_condition(reg.primary)
                # SymPy defines x**0 as 1 even at x=0; the catalogue-derived
                # source treats 0^0 as indeterminate.  SymPy semantics take
                # precedence in funcprops.
                return sp.Or(sp.Eq(expr.exp, 0), dom)
        return simplify_condition(_pow_real_domain(expr, var))

    # Ordinary arithmetic is defined wherever its varying operands are.
    # This exact fallback keeps basic domain reasoning independent of semialg.
    if real_mode and expr.func in (sp.Add, sp.Mul):
        ds = []
        for a in expr.args:
            if a.has(var):
                d = _function_domain_impl(a, var, domain=S.Reals)
                if isinstance(d, PropertyResult):
                    return d
                ds.append(d)
        return simplify_condition(sp.And(*ds)) if ds else S.true

    # Trigonometric reciprocal/pole domains are clearer and more useful in
    # structural form than the equivalent integer-membership catalogue rule.
    if (
        real_mode
        and expr.func in (sp.tan, sp.sec, sp.cot, sp.csc)
        and len(expr.args) == 1
    ):
        a = expr.args[0]
        inner = (
            _function_domain_impl(a, var, domain=S.Reals)
            if a.has(var) and a != var
            else S.true
        )
        denom = sp.cos(a) if expr.func in (sp.tan, sp.sec) else sp.sin(a)
        return simplify_condition(sp.And(inner, sp.Ne(denom, 0)))

    constraints = []

    # semialg is an exact backend for semialgebraic expressions only.
    if real_mode and is_semialgebraic_expression(expr, var):
        try:
            constraints.append(sp.sympify(_semialg.domain(expr, var)))
        except _semialg.SemialgUnavailable:
            pass
        except SYMBOLIC_QUERY_ERRORS:
            pass

    mode = "real" if real_mode else "complex"
    prop = "RealDomain" if real_mode else "ComplexDomain"
    reg = evaluate_property(expr, prop, var=var, mode=mode)
    if reg.primary is not None:
        constraints.append(sp.sympify(reg.primary))
    elif reg.cases:
        # Preserve parameter-dependent domains as a Boolean piecewise equivalent:
        # OR(guard AND domain-value). Unknown default stays unasserted.
        terms = []
        for case in reg.cases:
            if case.value in (S.true, S.false) or isinstance(
                case.value, (sp.logic.boolalg.Boolean, sp.core.relational.Relational)
            ):
                terms.append(sp.And(case.condition, case.value))
        if terms:
            constraints.append(sp.Or(*terms))

    # Explicit elementary rules not represented as standalone SymPy heads in the registry.
    if real_mode:
        if expr.func is sp.exp:
            constraints.append(S.true)
        elif expr.func is sp.log and len(expr.args) == 1:
            constraints.append(expr.args[0] > 0)
        elif expr.func is sp.asin or expr.func is sp.acos:
            a = expr.args[0]
            constraints.append(sp.And(a >= -1, a <= 1))
        elif expr.func is sp.acosh:
            constraints.append(expr.args[0] >= 1)
        elif expr.func is sp.atanh:
            a = expr.args[0]
            constraints.append(sp.And(a > -1, a < 1))

    # Every varying argument must itself be defined.
    for arg in getattr(expr, "args", ()):
        if arg == var or not getattr(arg, "has", lambda *_: False)(var):
            continue
        d = _function_domain_impl(
            arg, var, domain=S.Reals if real_mode else S.Complexes
        )
        if not isinstance(d, PropertyResult):
            constraints.append(sp.sympify(d))

    if not constraints:
        return PropertyResult(primary="unknown_domain", source="funcprops")
    try:
        result = simplify_condition(sp.And(*constraints))
    except SYMBOLIC_QUERY_ERRORS:
        result = sp.And(*constraints)
    if stratified and real_mode:
        return stratify_boolean_by_parameters(
            result, var, source="domain-parameter-stratification"
        )
    return result


def domain(expr, var, *, domain=S.Reals, stratified=False):
    """Cached public wrapper around the lazy real/complex domain prover."""
    from ..proof_cache import property_proof_context

    complex_mode = domain in (S.Complexes, "complex", "complexes", "C", "CC")
    if not (complex_mode or domain in (S.Reals, "real", "reals", "R", "RR")):
        raise ValueError("domain must be Reals or Complexes")
    return property_proof_context(expr, var).natural_domain(
        complex_domain=complex_mode, stratified=stratified
    )
