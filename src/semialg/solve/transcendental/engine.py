from __future__ import annotations

from dataclasses import dataclass, field

import sympy as sp
from sympy import S

from .cleanup import (
    finite_points_form,
    recon_solved_points,
    recon_univar_intv_form,
    remove_redundant_disjunc,
)
from .families import default_trans_handlers
from .periodic import (
    find_periodic_variables,
    recon_periodic_domain,
)
from .preprocess import prep_trans_problem
from .quantifier_elimination import eliminate_lead_block
from .roots import decomp_univar_inequality, isolate_univar_roots
from .state import TransProblemState, build_trans_state
from .system_roots import CompletenessCertificate, solve_bounded_trans_sys


@dataclass(frozen=True)
class TransReductionResult:
    state: TransProblemState
    formula: sp.Expr
    method: str
    complete: bool = False
    trace: tuple[str, ...] = ()
    completeness_certificate: CompletenessCertificate = CompletenessCertificate(
        False, "no_certificate", "none"
    )
    metadata: dict = field(default_factory=dict)


@dataclass(frozen=True)
class TranscendentalSolveTrace:
    steps: tuple[str, ...]
    metadata: dict = field(default_factory=dict)


def _domain_object(state: TransProblemState):
    return S.Reals if str(state.default_domain).lower().endswith("reals") else S.Complexes


def try_spec_fam_rewrites(state: TransProblemState) -> sp.Expr | None:
    variables = state.free_variables if state.free_variables else state.active_variable_order
    for handler in default_trans_handlers():
        if handler.rewrite_builder is None:
            continue
        try:
            rewritten = handler.rewrite_builder(state.formula, variables)
        except Exception:
            rewritten = None
        if rewritten is not None and rewritten != state.formula:
            return sp.simplify(rewritten)
    return None


def prior_univar_var(state: TransProblemState, prep) -> sp.Symbol | None:
    if prep.quantifier_plan and prep.quantifier_plan.prior_univar_vars:
        return prep.quantifier_plan.prior_univar_vars[0]
    return state.all_variables[0] if len(state.all_variables) == 1 else None


def reduce_trans_problem(state: TransProblemState) -> TransReductionResult:
    trace = []
    prep = prep_trans_problem(state, quantifier_aware=True)
    current = prep.state
    if prep.changed:
        trace.append("prep_trans_problem")
    if prep.quantifier_plan is not None:
        trace.append("quantifier_dispatch_plan")

    qe = eliminate_lead_block(current)
    if qe is not None:
        trace.append(qe.method)
        cert = CompletenessCertificate(
            qe.complete, "quantifier_elimination", qe.method, qe.metadata
        )
        return TransReductionResult(
            current,
            qe.resulting_formula,
            method=qe.method,
            complete=qe.complete,
            trace=tuple(trace),
            completeness_certificate=cert,
            metadata=qe.metadata,
        )

    rewritten = try_spec_fam_rewrites(current)
    if rewritten is not None:
        trace.append("special_family_rewrite")
        cleaned = remove_redundant_disjunc(rewritten)
        cert = CompletenessCertificate(True, "direct_family_rewrite", "special_family_rewrite")
        return TransReductionResult(
            current.with_formula(cleaned.cleaned),
            cleaned.cleaned,
            method="special_family_rewrite",
            complete=True,
            trace=tuple(trace),
            completeness_certificate=cert,
            metadata={"cleanup": cleaned.metadata},
        )

    uvar = prior_univar_var(current, prep)
    if uvar is not None:
        root_data = isolate_univar_roots(current.formula, uvar, domain=_domain_object(current))
        trace.append(root_data.method)
        if root_data.periodic_formula is not None:
            cleaned = remove_redundant_disjunc(root_data.periodic_formula)
            cert = CompletenessCertificate(
                root_data.complete,
                "periodic_root_reconstruction",
                root_data.method,
                root_data.metadata,
            )
            return TransReductionResult(
                current,
                cleaned.cleaned,
                method="univariate_periodic_roots",
                complete=root_data.complete,
                trace=tuple(trace),
                completeness_certificate=cert,
                metadata=root_data.metadata,
            )
        if root_data.roots:
            formula = finite_points_form((uvar,), [(r,) for r in root_data.roots])
            cleaned = remove_redundant_disjunc(formula)
            cert = CompletenessCertificate(
                root_data.complete,
                "root_isolation",
                root_data.method,
                {"certified_intervals": root_data.certified_intervals, **root_data.metadata},
            )
            return TransReductionResult(
                current,
                cleaned.cleaned,
                method="univariate_root_isolation",
                complete=root_data.complete,
                trace=tuple(trace),
                completeness_certificate=cert,
                metadata=root_data.metadata,
            )

        if _domain_object(current) == S.Reals:
            dec = decomp_univar_inequality(current.formula, uvar, domain=S.Reals)
            trace.append(dec.method)
            if dec.true_intervals:
                cleaned = recon_univar_intv_form(uvar, dec.true_intervals)
                cert = CompletenessCertificate(
                    False, "certified_interval_decomposition", dec.method, dec.metadata
                )
                return TransReductionResult(
                    current,
                    cleaned.cleaned,
                    method="univariate_inequality_decomposition",
                    complete=False,
                    trace=tuple(trace),
                    completeness_certificate=cert,
                    metadata={"intervals": dec.true_intervals, **dec.metadata},
                )
            periodic = find_periodic_variables(current)
            periodic_match = next(
                (p for p in periodic if p.variable == uvar and p.period is not None), None
            )
            if periodic_match is not None and dec.true_intervals:
                periodic_formula = recon_periodic_domain(
                    uvar, dec.true_intervals, periodic_match.period
                )
                cleaned = remove_redundant_disjunc(periodic_formula)
                cert = CompletenessCertificate(
                    False,
                    "periodic_interval_reconstruction",
                    dec.method,
                    {"period": periodic_match.period},
                )
                return TransReductionResult(
                    current,
                    cleaned.cleaned,
                    method="periodic_interval_reconstruction",
                    complete=False,
                    trace=tuple(trace),
                    completeness_certificate=cert,
                    metadata={"period": periodic_match.period},
                )

    periodic = find_periodic_variables(current)
    if periodic:
        trace.append("periodic_bounding_detected")

    fallback = solve_bounded_trans_sys(current.formula, current.all_variables)
    trace.append(fallback.method)
    if fallback.points:
        solved = recon_solved_points(current.all_variables, fallback.points)
        return TransReductionResult(
            current,
            solved.cleaned,
            method=fallback.method,
            complete=fallback.complete,
            trace=tuple(trace),
            completeness_certificate=fallback.completeness_certificate,
            metadata={"fallback": fallback.metadata, "certified_points": fallback.certified_points},
        )

    cleaned = remove_redundant_disjunc(current.formula)
    cert = CompletenessCertificate(
        False, "no_progress", "transcendental_no_progress", {"periodic_candidates": periodic}
    )
    return TransReductionResult(
        current,
        cleaned.cleaned,
        method="transcendental_no_progress",
        complete=False,
        trace=tuple(trace),
        completeness_certificate=cert,
        metadata={"periodic_candidates": periodic},
    )


__all__ = [
    "TransReductionResult",
    "TranscendentalSolveTrace",
    "reduce_trans_problem",
    "build_trans_state",
]
