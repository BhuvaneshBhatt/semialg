"""Shared parameter-stratification and three-valued proof helpers.

This module centralizes parameter-case normalization, semantic-predicate
lowering, and exact set comparisons so every property engine uses the same
conditional-result semantics.
"""

from __future__ import annotations

import contextlib
from collections.abc import Callable, Iterable
from typing import Any

import sympy as sp
from sympy import S

from .conditions import (
    AAlgebraicQ,
    AEqualQ,
    AEvenQ,
    AIntegerQ,
    ANegative,
    ANonIntegerQ,
    ANonNegative,
    ANonPositive,
    AOddQ,
    APositive,
    ARationalQ,
    ARealIfDefinedQ,
    ARealQ,
    AUnequalQ,
)
from .engine import simplify_condition
from .exceptions import SIMPLIFICATION_ERRORS, SYMBOLIC_QUERY_ERRORS
from .results import ConditionalExpression, PropertyResult, TruthValue


def semantic_condition(condition):
    """Lower catalogue semantic predicates to ordinary SymPy predicates.

    This is deliberately one-way.  Native inequalities/membership expressions
    are much easier for ``simplify_logic``, ``refine`` and set conversion to
    reason about, while the catalogue parser can continue to preserve source
    semantics internally.
    """
    c = sp.sympify(condition)

    def repl(node):
        if not isinstance(node, sp.Basic):
            return node
        f = node.func
        a = node.args
        if f is APositive and len(a) == 1:
            return a[0] > 0
        if f is ANegative and len(a) == 1:
            return a[0] < 0
        if f is ANonPositive and len(a) == 1:
            return a[0] <= 0
        if f is ANonNegative and len(a) == 1:
            return a[0] >= 0
        if f is AIntegerQ and len(a) == 1:
            return sp.Contains(a[0], S.Integers)
        if f is ANonIntegerQ and len(a) == 1:
            return sp.Not(sp.Contains(a[0], S.Integers))
        if f is ARationalQ and len(a) == 1:
            return sp.Contains(a[0], S.Rationals)
        if f in (ARealQ, ARealIfDefinedQ) and len(a) == 1:
            return sp.Contains(a[0], S.Reals)
        if f is AAlgebraicQ and len(a) == 1:
            # SymPy has no global Algebraics singleton that behaves as cleanly
            # as Reals/Integers in Boolean simplification.  Preserve this one.
            return node
        if f is AEvenQ and len(a) == 1:
            return sp.And(sp.Contains(a[0], S.Integers), sp.Eq(sp.Mod(a[0], 2), 0))
        if f is AOddQ and len(a) == 1:
            return sp.And(sp.Contains(a[0], S.Integers), sp.Eq(sp.Mod(a[0], 2), 1))
        if f is AUnequalQ and len(a) == 2:
            return sp.Ne(a[0], a[1])
        if f is AEqualQ and len(a) == 2:
            return sp.Eq(a[0], a[1])
        return node

    with contextlib.suppress(SYMBOLIC_QUERY_ERRORS):
        c = c.replace(
            lambda n: (
                isinstance(n, sp.Basic)
                and n.func
                in {
                    APositive,
                    ANegative,
                    ANonPositive,
                    ANonNegative,
                    AIntegerQ,
                    ANonIntegerQ,
                    ARationalQ,
                    ARealQ,
                    ARealIfDefinedQ,
                    AAlgebraicQ,
                    AEvenQ,
                    AOddQ,
                    AUnequalQ,
                    AEqualQ,
                }
            ),
            repl,
        )
    with contextlib.suppress(SYMBOLIC_QUERY_ERRORS):
        c = sp.refine(c)
    with contextlib.suppress(SYMBOLIC_QUERY_ERRORS):
        c = sp.simplify_logic(c)
    return c


def parameter_symbols(expr, var=None, *extra):
    syms = set(sp.sympify(expr).free_symbols)
    for obj in extra:
        with contextlib.suppress(SYMBOLIC_QUERY_ERRORS):
            syms.update(sp.sympify(obj).free_symbols)
    if var is not None:
        syms.discard(sp.sympify(var))
    return tuple(sorted(syms, key=sp.default_sort_key))


def normalize_cases(
    cases: Iterable[ConditionalExpression],
    *,
    source=None,
    primary=None,
    all_properties=(),
    evidence=(),
):
    """Canonicalize conditional cases and merge equal-valued strata."""
    merged: list[tuple[Any, sp.Basic]] = []
    for item in cases:
        cond = semantic_condition(simplify_condition(item.condition))
        if cond is S.false:
            continue
        found = False
        for i, (value, old) in enumerate(merged):
            if value == item.value:
                try:
                    new = sp.simplify(sp.simplify_logic(sp.Or(old, cond)))
                except SYMBOLIC_QUERY_ERRORS:
                    new = sp.Or(old, cond)
                merged[i] = (value, new)
                found = True
                break
        if not found:
            merged.append((item.value, cond))

    # Truth-valued strata are mutually exclusive claims.  If independently
    # derived TRUE and FALSE guards overlap, make that overlap explicitly
    # UNKNOWN.
    true_cond = S.false
    false_cond = S.false
    for value, cond in merged:
        if value is TruthValue.TRUE:
            true_cond = sp.Or(true_cond, cond)
        elif value is TruthValue.FALSE:
            false_cond = sp.Or(false_cond, cond)
    conflict = S.false
    if true_cond is not S.false and false_cond is not S.false:
        try:
            conflict = sp.simplify(sp.simplify_logic(sp.And(true_cond, false_cond)))
        except SYMBOLIC_QUERY_ERRORS:
            conflict = sp.And(true_cond, false_cond)
        if conflict is not S.false and not _guard_satisfiable(conflict):
            conflict = S.false
        if conflict is not S.false:
            adjusted = []
            for value, cond in merged:
                if value in (TruthValue.TRUE, TruthValue.FALSE):
                    try:
                        cond = sp.simplify(
                            sp.simplify_logic(sp.And(cond, sp.Not(conflict)))
                        )
                    except SYMBOLIC_QUERY_ERRORS:
                        cond = sp.And(cond, sp.Not(conflict))
                if cond is not S.false:
                    adjusted.append((value, cond))
            merged = [(TruthValue.UNKNOWN, conflict), *adjusted]

    out: list[ConditionalExpression] = []
    covered = S.false
    for value, cond in merged:
        try:
            effective = sp.simplify(sp.simplify_logic(sp.And(cond, sp.Not(covered))))
        except SYMBOLIC_QUERY_ERRORS:
            effective = sp.And(cond, sp.Not(covered))
        if effective is S.false:
            continue
        out.append(ConditionalExpression(value, effective))
        try:
            covered = sp.simplify(sp.simplify_logic(sp.Or(covered, cond)))
        except SYMBOLIC_QUERY_ERRORS:
            covered = sp.Or(covered, cond)

    if primary is None:
        for c in out:
            if c.condition is S.true:
                primary = c.value
                break
    return PropertyResult(
        primary=primary,
        cases=tuple(out),
        all_properties=tuple(all_properties),
        source=source,
        evidence=tuple(evidence),
    )


def resolve_result_assumptions(result, assumptions):
    """Filter a parameter-stratified result through known assumptions.

    A uniquely implied stratum is returned as its compact value.  Otherwise
    impossible strata are removed and the remaining conditional result keeps
    its original proof metadata.
    """
    if not isinstance(result, PropertyResult):
        return result
    assumptions = semantic_condition(assumptions)
    if assumptions is S.true or not result.cases:
        return result

    from .assumptions import entails

    remaining = []
    implied = []
    for case in result.cases:
        if entails(case.condition, assumptions) is True:
            implied.append(case.value)
            continue
        if entails(sp.Not(case.condition, evaluate=False), assumptions) is True:
            continue
        remaining.append(case)

    if implied:
        first = implied[0]
        if all(value == first for value in implied):
            return first
    if not remaining:
        return result.primary if result.primary is not None else result
    return normalize_cases(
        remaining,
        source=result.source,
        primary=result.primary
        if result.primary in {c.value for c in remaining}
        else None,
        all_properties=result.all_properties,
        evidence=result.evidence,
    )


def map_conditional(value, mapper: Callable[[Any], Any], *, source=None):
    if not isinstance(value, PropertyResult):
        return mapper(value)
    cases = [ConditionalExpression(mapper(c.value), c.condition) for c in value.cases]
    primary = mapper(value.primary) if value.primary is not None else None
    return normalize_cases(
        cases, primary=primary, source=source or value.source, evidence=value.evidence
    )


def implication_guard(domain, locus, var):
    """Return a parameter-only sufficient guard for ``domain => locus``.

    ``None`` means that the implication could not be reduced without the
    function variable.  The returned guard is exact for the structural forms
    handled here; we never existentially sample parameters.
    """
    var = sp.sympify(var)
    d = semantic_condition(domain)
    locus_cond = semantic_condition(locus)
    try:
        imp = sp.simplify_logic(sp.Implies(d, locus_cond))
    except SYMBOLIC_QUERY_ERRORS:
        imp = sp.Implies(d, locus_cond)
    with contextlib.suppress(SYMBOLIC_QUERY_ERRORS):
        imp = sp.refine(imp, d)
    with contextlib.suppress(SYMBOLIC_QUERY_ERRORS):
        imp = sp.simplify_logic(imp)
    if not getattr(imp, "has", lambda *_: False)(var):
        return imp

    # Common catalogue form: conjunction of domain facts and parameter facts.
    if isinstance(locus_cond, sp.And):
        residual = []
        for term in locus_cond.args:
            try:
                test = sp.simplify_logic(sp.Implies(d, term))
            except SYMBOLIC_QUERY_ERRORS:
                test = sp.Implies(d, term)
            if test is S.true:
                continue
            if not term.has(var):
                residual.append(term)
            else:
                return None
        return sp.And(*residual) if residual else S.true
    return None


def exact_set_equal(a, b):
    """Three-valued exact set equality test."""
    if a == b:
        return TruthValue.TRUE
    if not isinstance(a, sp.Set) or not isinstance(b, sp.Set):
        return TruthValue.UNKNOWN
    try:
        diff = a.symmetric_difference(b)
        if diff is S.EmptySet or diff == S.EmptySet:
            return TruthValue.TRUE
        # Nonempty finite/interval/union differences are exact disproofs.
        if isinstance(diff, (sp.FiniteSet, sp.Interval)):
            return TruthValue.FALSE
        emp = getattr(diff, "is_empty", None)
        if emp is True:
            return TruthValue.TRUE
        if emp is False:
            return TruthValue.FALSE
    except SYMBOLIC_QUERY_ERRORS:
        pass
    try:
        ab = a.is_subset(b)
        ba = b.is_subset(a)
        if ab is True and ba is True:
            return TruthValue.TRUE
        if ab is False or ba is False:
            return TruthValue.FALSE
    except SYMBOLIC_QUERY_ERRORS:
        pass
    return TruthValue.UNKNOWN


def combine_truth_results(a, b, *, op="and", source="combined"):
    """Combine scalar/conditional TruthValue results without losing strata."""

    def scalar(x, y):
        if op == "and":
            if x is TruthValue.FALSE or y is TruthValue.FALSE:
                return TruthValue.FALSE
            if x is TruthValue.TRUE and y is TruthValue.TRUE:
                return TruthValue.TRUE
            return TruthValue.UNKNOWN
        if x is TruthValue.TRUE or y is TruthValue.TRUE:
            return TruthValue.TRUE
        if x is TruthValue.FALSE and y is TruthValue.FALSE:
            return TruthValue.FALSE
        return TruthValue.UNKNOWN

    if not isinstance(a, PropertyResult) and not isinstance(b, PropertyResult):
        return scalar(a, b)
    acases = (
        a.cases
        if isinstance(a, PropertyResult)
        else (ConditionalExpression(a, S.true),)
    )
    bcases = (
        b.cases
        if isinstance(b, PropertyResult)
        else (ConditionalExpression(b, S.true),)
    )
    cases = []
    for ca in acases:
        for cb in bcases:
            try:
                cond = sp.simplify_logic(sp.And(ca.condition, cb.condition))
            except SYMBOLIC_QUERY_ERRORS:
                cond = sp.And(ca.condition, cb.condition)
            if cond is not S.false:
                cases.append(ConditionalExpression(scalar(ca.value, cb.value), cond))
    return normalize_cases(cases, source=source)


def _guard_satisfiable(guard):
    """Return False only when a parameter guard is certified impossible.

    SymPy's Boolean satisfier is propositional for inequalities, so polynomial
    real guards are sent to semialg first.  Unsupported guards fall back to
    conservative propositional reasoning.
    """
    guard = semantic_condition(guard)
    if guard is S.false:
        return False
    if guard is S.true:
        return True
    variables = tuple(sorted(guard.free_symbols, key=sp.default_sort_key))
    if variables:
        try:
            from . import _semialg

            return _semialg.is_satisfiable(guard, variables)
        except (
            _semialg.SemialgUnavailable,
            TypeError,
            ValueError,
            NotImplementedError,
            ArithmeticError,
        ):
            pass
    try:
        return sp.satisfiable(guard) is not False
    except SIMPLIFICATION_ERRORS:
        return True


def stratify_boolean_by_parameters(
    condition, var, *, max_atoms=8, source="boolean-parameter-stratification"
):
    """Split a Boolean condition into parameter guards and var-only values.

    This is exact for the finite Boolean atom partition it constructs.  It is
    especially useful for catalogue domains stored as a single mixed formula,
    e.g. ``IntegerQ(nu) | x>0 | (nu>=0 & x>=0)``.
    """
    var = sp.sympify(var)
    c = semantic_condition(condition)
    # Parameter-only Boolean atoms.  Do not split arbitrary function-valued
    # endpoint expressions such as gamma(a); only predicates/relations matter.
    try:
        atoms = set(c.atoms(sp.core.relational.Relational, sp.Contains))
    except SYMBOLIC_QUERY_ERRORS:
        atoms = set()
    atoms = {a for a in atoms if not a.has(var)}
    atoms = tuple(sorted(atoms, key=sp.default_sort_key))
    if not atoms or len(atoms) > max_atoms:
        return PropertyResult(primary=c, source=source)
    cases = []
    from itertools import product

    for bits in product((False, True), repeat=len(atoms)):
        guard_terms = []
        repl = {}
        for a, b in zip(atoms, bits, strict=False):
            repl[a] = S.true if b else S.false
            guard_terms.append(a if b else sp.Not(a))
        guard = semantic_condition(sp.And(*guard_terms))
        if not _guard_satisfiable(guard):
            continue
        try:
            value = c.xreplace(repl)
            value = semantic_condition(value)
            with contextlib.suppress(SYMBOLIC_QUERY_ERRORS):
                value = sp.simplify(value)
        except SYMBOLIC_QUERY_ERRORS:
            value = c
        # If parameters remain in predicates, this partition was insufficient.
        params = set(value.free_symbols)
        params.discard(var)
        if params:
            continue
        cases.append(ConditionalExpression(value, guard))
    if not cases:
        return PropertyResult(primary=c, source=source)
    return normalize_cases(cases, source=source)


def result_cases(result):
    """Return a uniform conditional-case view of a scalar or PropertyResult."""
    if isinstance(result, PropertyResult):
        if result.cases:
            return result.cases
        if result.primary is not None:
            return (ConditionalExpression(result.primary, S.true),)
        return ()
    return (ConditionalExpression(result, S.true),)


def guarded_result(result, guard, *, source=None):
    """Attach a parameter guard to every case in a property result."""
    guard = semantic_condition(guard)
    cases = []
    for case in result_cases(result):
        try:
            cond = sp.simplify_logic(sp.And(guard, case.condition))
        except SYMBOLIC_QUERY_ERRORS:
            cond = sp.And(guard, case.condition)
        if cond is not S.false:
            cases.append(ConditionalExpression(case.value, cond))
    if not cases:
        return PropertyResult(primary=None, source=source)
    evidence = result.evidence if isinstance(result, PropertyResult) else ()
    exceptions = result.exceptions if isinstance(result, PropertyResult) else ()
    all_props = result.all_properties if isinstance(result, PropertyResult) else ()
    out = normalize_cases(
        cases,
        source=source or getattr(result, "source", None),
        all_properties=all_props,
        evidence=evidence,
    )
    if exceptions:
        return PropertyResult(
            out.primary,
            cases=out.cases,
            all_properties=out.all_properties,
            source=out.source,
            exceptions=exceptions,
            evidence=out.evidence,
        )
    return out


def evaluate_parameter_strata(
    condition,
    var,
    evaluator: Callable[[sp.Basic], Any],
    *,
    source="parameter-stratification",
    max_atoms=8,
):
    """Evaluate one property separately on exact parameter-domain strata.

    ``condition`` may mix the function variable with parameter-only Boolean
    atoms.  When the finite parameter partition is exact, ``evaluator`` is
    called on each resulting variable-only domain and nested conditional
    results are flattened.  ``None`` means there was no useful split.
    """
    split = stratify_boolean_by_parameters(
        condition, var, max_atoms=max_atoms, source=source
    )
    if not split.cases or (
        len(split.cases) == 1 and split.cases[0].condition is S.true
    ):
        return None
    cases = []
    evidence = []
    for domain_case in split.cases:
        value = evaluator(domain_case.value)
        if isinstance(value, PropertyResult):
            evidence.extend(value.evidence)
        for item in result_cases(value):
            try:
                guard = sp.simplify_logic(sp.And(domain_case.condition, item.condition))
            except SYMBOLIC_QUERY_ERRORS:
                guard = sp.And(domain_case.condition, item.condition)
            if guard is not S.false:
                cases.append(ConditionalExpression(item.value, guard))
    if not cases:
        return None
    return normalize_cases(cases, source=source, evidence=tuple(evidence))
