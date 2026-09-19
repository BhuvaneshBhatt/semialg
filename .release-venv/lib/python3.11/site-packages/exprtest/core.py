"""Latency-sensitive public zero-testing oracle."""

from __future__ import annotations

import random

import sympy as sp

from . import _cache
from . import _config as cfg
from ._algebraic import exact_algebraic_number_test, tower_algebraic_test
from ._assumptions import (
    has_nontrivial_assumptions,
    normalize_assumptions,
    refine_with_assumptions,
)
from ._call_memo import OracleMemo
from ._cost import ExactBudget, within_budget
from ._cyclotomic import cyclotomic_zero_test
from ._defined import quick_defined
from ._domains import NumberKind, number_kind
from ._fast import literal_number_zero, quick_reduce
from ._features import ExprFeatures, expression_features
from ._finite_field import finite_field_identity_test
from ._flint_values import classify_flint_value
from ._identity import elementary_identity_normal_form
from ._negative_cache import (
    cyclotomic_inapplicable,
    log_inapplicable,
    pslq_inapplicable,
    special_inapplicable,
    sqrt_sum_inapplicable,
    tower_inapplicable,
)
from ._nonzero import quick_nonzero
from ._normal_forms import exp_log_normal_form, special_normal_form, trig_normal_form
from ._numeric import adaptive_precision_ball_test, random_witness_nonzero_check
from ._profile import StageProfiler, ZeroTestProfile
from ._pslq import pslq_zero_relation
from ._result import Verdict, ZeroClassification
from ._rewrite import exact_rewrite
from ._sqrt_sum import square_root_sum_test
from ._timing import run_with_time_budget
from ._transcendental import transcendental_zero_test


def _fast_unknown(detail: str) -> ZeroClassification:
    return ZeroClassification(Verdict.UNKNOWN, "fast-oracle", detail=detail)


def _verify_pslq(term: sp.Expr, assumptions) -> bool | None:
    """Verify a PSLQ candidate using only independent exact proof stages."""
    term = quick_reduce(sp.sympify(term))
    if term == 0 or term.is_zero is True:
        return True
    special = special_normal_form(term)
    if special == 0 or special.is_zero is True:
        return True
    normalized = exp_log_normal_form(special, assumptions)
    if normalized == 0 or normalized.is_zero is True:
        return True
    cyclo = cyclotomic_zero_test(normalized)
    if cyclo.verdict is Verdict.ZERO_PROVEN:
        return True
    exact_term = trig_normal_form(normalized)
    alg = run_with_time_budget(exact_algebraic_number_test, exact_term,
                               seconds=cfg.EXACT_STAGE_TIMEOUT, default=None)
    return True if alg is not None and alg.verdict is Verdict.ZERO_PROVEN else None


def _cheap_zero_property(term: sp.Expr, features: ExprFeatures | None = None) -> bool | None:
    """Use SymPy's zero property only on structurally tiny expressions."""
    features = features or expression_features(term)
    if features.has_float or features.has_nonfinite:
        return None
    if features.ops > 4 or features.nodes > 12:
        return None
    value = term.is_zero
    if value is None:
        return None
    if value is False and term.is_finite is not True:
        return None
    return bool(value)


def _stage_run(profiler: StageProfiler | None, stage: str, func, outcome=None):
    if profiler is None:
        return func()
    return profiler.run(stage, func, outcome)


def _verdict_name(result: ZeroClassification) -> str:
    return result.verdict.name.lower()


def _fast_exact(term: sp.Expr, assumptions,
                profiler: StageProfiler | None = None,
                memo: OracleMemo | None = None) -> ZeroClassification:
    """Classify a closed exact expression without generic simplification."""
    term = _stage_run(profiler, "quick-reduce", lambda: quick_reduce(term))
    literal_zero = literal_number_zero(term)
    if literal_zero is not None:
        return ZeroClassification(
            Verdict.ZERO_PROVEN if literal_zero else Verdict.NONZERO_PROVEN,
            "quick-arithmetic", detail="local exact arithmetic decided the expression",
            evidence="exact-arithmetic",
        )
    features = _stage_run(profiler, "features", lambda: expression_features(term))
    if features.has_domain_hazard and quick_defined(term, assumptions) is False:
        return ZeroClassification(
            Verdict.UNKNOWN, "definedness",
            detail="expression is undefined or has a pole under the supplied assumptions",
            evidence="undefined-expression",
        )
    identity = _stage_run(
        profiler, "elementary-identity", lambda: elementary_identity_normal_form(term, assumptions)
    )
    if identity != term:
        term = quick_reduce(identity)
        literal_zero = literal_number_zero(term)
        if literal_zero is not None:
            return ZeroClassification(
                Verdict.ZERO_PROVEN if literal_zero else Verdict.NONZERO_PROVEN,
                "elementary-identity",
                detail="bounded exact polynomial/rational or elementary identity normalization decided the expression",
                evidence="exact-identity",
            )
        features = _stage_run(profiler, "features-identity", lambda: expression_features(term))
    if memo is None and features.nodes >= cfg.CALL_MEMO_MIN_NODES:
        memo = OracleMemo()
        memo.feature_cache = {term: features}
    if features.has_rewrite:
        rewritten = _stage_run(profiler, "exact-rewrite", lambda: exact_rewrite(term))
        if rewritten != term:
            term = rewritten
            literal_zero = literal_number_zero(term)
            if literal_zero is not None:
                return ZeroClassification(
                    Verdict.ZERO_PROVEN if literal_zero else Verdict.NONZERO_PROVEN,
                    "exact-rewrite", detail="exact local rewrite produced a finite exact number",
                    evidence="exact-rewrite",
                )
            features = _stage_run(
                profiler, "features-rewrite",
                lambda: memo.features(term) if memo is not None else expression_features(term),
            )
    elif profiler is not None:
        profiler.note("exact-rewrite", "skipped", "structurally inapplicable")
    if not special_inapplicable(term):
        special = _stage_run(profiler, "special-values", lambda: special_normal_form(term))
        if special != term:
            term = special
            features = _stage_run(
                profiler, "features-special",
                lambda: memo.features(term) if memo is not None else expression_features(term),
            )
    elif profiler is not None:
        profiler.note("special-values", "skipped", "structurally inapplicable")
    literal_zero = literal_number_zero(term)
    if literal_zero is not None:
        return ZeroClassification(
            Verdict.ZERO_PROVEN if literal_zero else Verdict.NONZERO_PROVEN,
            "quick-arithmetic", detail="local exact arithmetic decided the expression",
            evidence="exact-arithmetic",
        )
    if not within_budget(term, features=features):
        return _fast_unknown("exact-expression budget exceeded")
    nonzero = _stage_run(profiler, "quick-nonzero", lambda: (memo.nonzero(term, assumptions) if memo is not None else quick_nonzero(term, assumptions)),
                         lambda value: "nonzero" if value is True else "zero" if value is False else "unknown")
    if nonzero is True:
        return ZeroClassification(Verdict.NONZERO_PROVEN, "quick-nonzero",
                                  detail="a cheap structural theorem proves the exact value is nonzero",
                                  evidence="structural-nonzero")
    if nonzero is False:
        return ZeroClassification(Verdict.ZERO_PROVEN, "quick-nonzero",
                                  detail="a cheap structural theorem proves the exact value is zero",
                                  evidence="structural-zero")
    zero_prop = _stage_run(profiler, "zero-property",
                           lambda: _cheap_zero_property(term, features))
    if zero_prop is True:
        return ZeroClassification(Verdict.ZERO_PROVEN, "zero-property",
                                  detail="cheap SymPy exact zero property is true", evidence="exact-property")
    if zero_prop is False:
        return ZeroClassification(Verdict.NONZERO_PROVEN, "zero-property",
                                  detail="cheap SymPy exact zero property is false", evidence="exact-property")
    kind = _stage_run(profiler, "number-kind", lambda: (memo.kind(term) if memo is not None else number_kind(term)),
                      lambda value: value.value)
    if kind is NumberKind.TRANSCENDENTAL:
        return ZeroClassification(Verdict.NONZERO_PROVEN, "number-kind",
                                  detail="a cheap transcendence theorem proves the exact constant is nonzero",
                                  evidence="proven-transcendental")
    if not features.has_log_exp or log_inapplicable(term):
        trans = _fast_unknown("no transcendental structural rule applies")
        if profiler is not None:
            profiler.note("transcendental", "skipped", "no log/exp structure")
    else:
        trans = _stage_run(profiler, "transcendental",
                           lambda: transcendental_zero_test(term, assumptions), _verdict_name)
    if trans.verdict is not Verdict.UNKNOWN:
        return trans
    normalized = term
    if features.has_log_exp and not log_inapplicable(term):
        normalized = _stage_run(profiler, "exp-log-normalize",
                                lambda: exp_log_normal_form(term, assumptions))
    elif profiler is not None:
        profiler.note("exp-log-normalize", "skipped", "structurally inapplicable")
    if normalized != term:
        term = quick_reduce(normalized)
        features = _stage_run(
            profiler, "features-normalized",
            lambda: memo.features(term) if memo is not None else expression_features(term),
        )
        literal_zero = literal_number_zero(term)
        if literal_zero is not None:
            return ZeroClassification(Verdict.ZERO_PROVEN if literal_zero else Verdict.NONZERO_PROVEN,
                                      "exp-log-normalize",
                                      detail="branch-safe exp/log normalization produced a finite exact number",
                                      evidence="exact-exp-log-identity")
        kind = _stage_run(profiler, "number-kind-normalized", lambda: (memo.kind(term) if memo is not None else number_kind(term)),
                          lambda value: value.value)
        if kind is NumberKind.TRANSCENDENTAL:
            return ZeroClassification(Verdict.NONZERO_PROVEN, "number-kind",
                                      detail="normalization exposed a provably transcendental exact constant",
                                      evidence="proven-transcendental")
    if not features.has_cyclotomic_shape or cyclotomic_inapplicable(term):
        cyclo = _fast_unknown("cyclotomic stage structurally inapplicable")
        if profiler is not None:
            profiler.note("cyclotomic", "skipped", "structurally inapplicable")
    else:
        cyclo = _stage_run(profiler, "cyclotomic", lambda: cyclotomic_zero_test(term), _verdict_name)
    if cyclo.verdict is not Verdict.UNKNOWN:
        return cyclo

    if sqrt_sum_inapplicable(term):
        sqrt_sum = _fast_unknown("square-root sum stage structurally inapplicable")
        if profiler is not None:
            profiler.note("sqrt-sum", "skipped", "structurally inapplicable")
    else:
        sqrt_sum = _stage_run(profiler, "sqrt-sum", lambda: square_root_sum_test(term), _verdict_name)
    if sqrt_sum.verdict is not Verdict.UNKNOWN:
        return sqrt_sum

    # The tower stage is cheaper than primitive-element construction for small
    # nested radicals. It is admitted only under a stricter fast-path budget.
    if (not features.has_tower_shape and not term.has(sp.I)) or tower_inapplicable(term):
        tower = _fast_unknown("tower stage structurally inapplicable")
        if profiler is not None:
            profiler.note("algebraic-tower", "skipped", "structurally inapplicable")
    else:
        tower = _stage_run(
            profiler, "algebraic-tower",
            lambda: run_with_time_budget(
                tower_algebraic_test, term,
                seconds=cfg.TOWER_SIGN_TIMEOUT,
                default=_fast_unknown("tower stage timed out"),
            ),
            _verdict_name,
        )
    if tower.verdict is not Verdict.UNKNOWN:
        return tower

    exact_term = _stage_run(profiler, "exact-trig", lambda: trig_normal_form(term))
    alg = _stage_run(
        profiler, "exact-algebraic",
        lambda: run_with_time_budget(exact_algebraic_number_test, exact_term,
                                     seconds=cfg.EXACT_STAGE_TIMEOUT, default=None),
        lambda value: "timeout" if value is None else _verdict_name(value),
    )
    if alg is not None and alg.verdict is not Verdict.UNKNOWN:
        return alg
    if pslq_inapplicable(term):
        relation = None
        if profiler is not None:
            profiler.note("pslq", "skipped", "structurally inapplicable")
    else:
        relation = _stage_run(
            profiler, "pslq",
            lambda: pslq_zero_relation(term, lambda candidate: _verify_pslq(candidate, assumptions), ExactBudget()),
            lambda value: "verified" if value is not None else "no-relation",
        )
    if relation is not None:
        return ZeroClassification(Verdict.ZERO_PROVEN, "pslq-verified",
                                  detail="PSLQ proposed the target relation and independent exact verification proved it",
                                  evidence="pslq-candidate-exact-verification")
    numerical = _stage_run(profiler, "arb", lambda: adaptive_precision_ball_test(term), _verdict_name)
    if numerical.verdict is not Verdict.UNKNOWN:
        return numerical
    return _fast_unknown("bounded exact and rigorous numeric stages were inconclusive")

def _fast_classify_uncached(term: sp.Expr, assumptions, rng: random.Random,
                            profiler: StageProfiler | None = None,
                            memo: OracleMemo | None = None) -> ZeroClassification:
    """Run the uncached bounded oracle after input normalization."""
    if assumptions is sp.false:
        return ZeroClassification(Verdict.UNKNOWN, "assumptions",
                                  detail="assumptions are inconsistent, so the requested domain is empty",
                                  evidence="empty-domain")
    refined = refine_with_assumptions(term, assumptions)
    if refined != term:
        term = refined
    literal_zero = literal_number_zero(term)
    if literal_zero is not None:
        return ZeroClassification(Verdict.ZERO_PROVEN if literal_zero else Verdict.NONZERO_PROVEN,
                                  "literal", detail="expression is a finite exact SymPy number",
                                  evidence="exact-arithmetic")
    has_free = bool(term.free_symbols)
    if has_free:
        free_features = expression_features(term)
        if free_features.has_domain_hazard and quick_defined(term, assumptions) is False:
            return ZeroClassification(
                Verdict.UNKNOWN, "definedness",
                detail="expression is undefined or has a pole under the supplied assumptions",
                evidence="undefined-expression",
            )
        if free_features.has_rewrite and free_features.ops <= cfg.EXACT_REWRITE_MAX_OPS:
            rewritten = exact_rewrite(term)
            if rewritten != term:
                term = rewritten
                literal_zero = literal_number_zero(term)
                if literal_zero is not None:
                    return ZeroClassification(Verdict.ZERO_PROVEN if literal_zero else Verdict.NONZERO_PROVEN,
                                              "exact-rewrite", detail="exact local rewrite produced a finite exact number",
                                              evidence="exact-rewrite")
                has_free = bool(term.free_symbols)
                free_features = expression_features(term) if has_free else None
        if has_free and memo is None and free_features is not None and free_features.nodes >= cfg.CALL_MEMO_MIN_NODES:
            memo = OracleMemo()
    if not has_free:
        return _fast_exact(term, assumptions, profiler, memo)
    identity = _stage_run(
        profiler, "elementary-identity", lambda: elementary_identity_normal_form(term, assumptions)
    )
    if identity != term:
        term = quick_reduce(identity)
        literal_zero = literal_number_zero(term)
        if literal_zero is not None:
            return ZeroClassification(
                Verdict.ZERO_PROVEN if literal_zero else Verdict.NONZERO_PROVEN,
                "elementary-identity",
                detail="bounded exact polynomial/rational or elementary identity normalization decided the expression",
                evidence="exact-identity",
            )
        has_free = bool(term.free_symbols)
        if not has_free:
            return _fast_exact(term, assumptions, profiler, memo)
    structural = (memo.nonzero(term, assumptions) if memo is not None
                  else quick_nonzero(term, assumptions))
    if structural is True:
        return ZeroClassification(Verdict.NONZERO_PROVEN, "quick-nonzero",
                                  detail="a cheap structural theorem proves nonzero",
                                  evidence="structural-nonzero")
    if structural is False:
        return ZeroClassification(Verdict.ZERO_PROVEN, "quick-nonzero",
                                  detail="a cheap structural theorem proves zero",
                                  evidence="structural-zero")
    if term.is_Mul:
        unresolved = []
        for factor in term.args:
            sub = _fast_classify_uncached(factor, assumptions, rng, profiler, memo)
            if sub.verdict is Verdict.ZERO_PROVEN:
                return ZeroClassification(Verdict.ZERO_PROVEN, "product-factor",
                                          detail="a product factor is exactly zero", evidence="factor-proof")
            if sub.verdict is not Verdict.NONZERO_PROVEN:
                unresolved.append(factor)
        if not unresolved:
            return ZeroClassification(Verdict.NONZERO_PROVEN, "product-factor",
                                      detail="every product factor is proven nonzero", evidence="factor-proofs")
    probe = random_witness_nonzero_check(term, term.free_symbols, assumptions=assumptions, rng=rng)
    if probe is not None and probe.verdict is Verdict.NONZERO_PROVEN:
        return probe
    if not has_nontrivial_assumptions(assumptions):
        sz = finite_field_identity_test(term, sorted(term.free_symbols, key=str), rng=rng)
        if sz.verdict in (Verdict.ZERO_PROVEN, Verdict.NONZERO_PROVEN):
            return sz
    if probe is not None and probe.verdict is Verdict.NONZERO_LIKELY:
        return probe
    return _fast_unknown("fast symbolic stages were inconclusive")


def zerotest(expr, assumptions=True, use_cache: bool = True, *,
             rng: random.Random | None = None, seed: int | None = None,
             confidence: str = "probable") -> bool | None:
    """Fast proof-oriented zero oracle returning ``True``, ``False``, or ``None``.

    The normal oracle path never invokes ``sympy.simplify`` or the general
    symbolic fallback. Potentially expensive exact methods are budgeted.
    """
    if rng is not None and seed is not None:
        raise ValueError("pass either rng or seed, not both")
    if confidence not in ("certified", "probable"):
        raise ValueError("confidence must be 'certified' or 'probable'")
    rng = rng or random.Random(seed)
    direct = classify_flint_value(expr)
    if direct is not None:
        result = direct
    else:
        term = sp.sympify(expr)
        assumptions = normalize_assumptions(assumptions)
        key = _cache.make_cache_key(term, assumptions, cfg.SZ_TARGET_FALSE_POSITIVE)
        result = _cache.cache_get(key) if use_cache else None
        if result is None:
            result = _fast_classify_uncached(term, assumptions, rng)
            if use_cache:
                _cache.cache_set(key, result)
    if result.verdict is Verdict.ZERO_PROVEN:
        return True
    if result.verdict is Verdict.NONZERO_PROVEN:
        return False
    if confidence == "probable" and result.verdict is Verdict.NONZERO_LIKELY:
        return False
    return None


def profile_zerotest(expr, assumptions=True, *,
                     rng: random.Random | None = None,
                     seed: int | None = None,
                     confidence: str = "probable") -> ZeroTestProfile:
    """Run the fast oracle with opt-in per-stage timings.

    This bypasses the final-result cache so the trace reflects actual stage
    work. The ordinary :func:`zerotest` API remains unchanged and does not
    allocate profiling records.
    """
    if rng is not None and seed is not None:
        raise ValueError("pass either rng or seed, not both")
    if confidence not in ("certified", "probable"):
        raise ValueError("confidence must be 'certified' or 'probable'")
    rng = rng or random.Random(seed)
    profiler = StageProfiler()
    direct = profiler.run("flint-direct", lambda: classify_flint_value(expr),
                          lambda value: "not-flint" if value is None else _verdict_name(value))
    if direct is not None:
        result = direct
    else:
        term = sp.sympify(expr)
        assumptions = normalize_assumptions(assumptions)
        result = _fast_classify_uncached(term, assumptions, rng, profiler)
    value = (True if result.verdict is Verdict.ZERO_PROVEN else
             False if result.verdict is Verdict.NONZERO_PROVEN else
             False if confidence == "probable" and result.verdict is Verdict.NONZERO_LIKELY else
             None)
    return profiler.finish(value, result)
