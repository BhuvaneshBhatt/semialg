"""Unified replay and diagnostic views for exact semialg results."""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass, field

import sympy as sp


@dataclass(frozen=True)
class CertificateReplayResult:
    """Outcome of independently replaying an exact result certificate."""

    verified: bool | None
    kind: str
    method: str
    details: Mapping[str, object] = field(default_factory=dict)


@dataclass(frozen=True)
class ResultDiagnostics:
    """Uniform public diagnostic summary across solver result families."""

    kind: str
    method: str
    certified: bool | None
    counters: Mapping[str, int] = field(default_factory=dict)
    cache_stats: Mapping[str, object] = field(default_factory=dict)
    details: Mapping[str, object] = field(default_factory=dict)


def _bool_equal(left: object, right: object) -> bool:
    return left is right or left == right


def replay_certificate(result: object) -> CertificateReplayResult:
    """Replay a supported exact certificate using an independent public path.

    Replay deliberately prefers semantic rechecking over trusting stored flags.
    ``verified=None`` means that the result type does not carry enough
    public replay payload; it never means that a certificate failed.
    """

    from .model import CADResult, QEResult
    from .optimization_results import FunctionRangeResult, OptimizationResult

    if isinstance(result, CADResult):
        from .invariants import validate_cad_result

        issues = tuple(validate_cad_result(result))
        return CertificateReplayResult(
            not issues,
            "cad",
            result.projection_mode,
            {"invariant_issues": issues, "cell_levels": len(result.cells_by_level)},
        )

    if isinstance(result, QEResult):
        cad_replay = replay_certificate(result.cad)
        return CertificateReplayResult(
            cad_replay.verified,
            "quantifier_elimination",
            result.cad.projection_mode,
            {"cad": cad_replay, "is_sentence": result.is_sentence},
        )

    if isinstance(result, OptimizationResult):
        condition = result.diagnostics.get("constraint_formula")
        if condition is None:
            return CertificateReplayResult(
                None,
                "optimization",
                result.method,
                {"reason": "constraint_formula is unavailable for this result"},
            )
        from .optimization import semialgebraic_maximize, semialgebraic_minimize

        solver = semialgebraic_minimize if result.kind == "min" else semialgebraic_maximize
        replayed = solver(
            result.objective,
            condition,
            result.variables,
            return_result=True,
            certification="complete" if result.certified else "auto",
        )
        verified = sp.simplify(replayed.value - result.value) == 0
        if result.attained and result.points:
            verified = verified and all(
                sp.simplify(result.objective.subs(point) - result.value) == 0
                for point in result.points
            )
        return CertificateReplayResult(
            bool(verified),
            "optimization",
            result.method,
            {"replayed_method": replayed.method, "replayed_certified": replayed.certified},
        )

    if isinstance(result, FunctionRangeResult):
        # The range formula itself is a replayable certificate of the reported
        # endpoints.  Check endpoint membership/limit orientation where finite.
        value = result.value_symbol
        checks: list[bool] = []
        if result.infimum is not None and result.infimum not in (-sp.oo, sp.oo):
            at_lower = sp.simplify(result.formula.subs(value, result.infimum))
            if result.minimum_attained is True:
                checks.append(_bool_equal(at_lower, sp.true))
        if result.supremum is not None and result.supremum not in (-sp.oo, sp.oo):
            at_upper = sp.simplify(result.formula.subs(value, result.supremum))
            if result.maximum_attained is True:
                checks.append(_bool_equal(at_upper, sp.true))
        return CertificateReplayResult(
            all(checks) if checks else None,
            "function_range",
            result.method,
            {"endpoint_checks": tuple(checks)},
        )

    # Convexity certificate imports are intentionally lazy to avoid adding a
    # heavy CAD import to the common result path.
    try:
        from .convexity import (
            ConvexityCertificate,
            PolynomialConvexityCertificate,
            QuadraticConvexityCertificate,
            convexity_certificate,
            polynomial_convexity_certificate,
            quadratic_convexity_certificate,
        )

        if isinstance(result, ConvexityCertificate):
            fresh = convexity_certificate(result.formula, result.variables)
            return CertificateReplayResult(
                fresh.outcome == result.outcome,
                "convexity",
                result.method,
                {"replayed_method": fresh.method},
            )
        if isinstance(result, PolynomialConvexityCertificate):
            fresh = polynomial_convexity_certificate(
                result.expression,
                result.variables,
                domain=result.domain,
                sense=result.sense,
            )
            return CertificateReplayResult(
                fresh.outcome == result.outcome,
                "polynomial_convexity",
                result.method,
                {"replayed_method": fresh.method},
            )
        if isinstance(result, QuadraticConvexityCertificate):
            fresh = quadratic_convexity_certificate(result.formula, result.variables)
            return CertificateReplayResult(
                fresh.outcome == result.outcome,
                "quadratic_convexity",
                result.method,
                {"replayed_method": fresh.method},
            )
    except (ImportError, TypeError, ValueError, NotImplementedError, ArithmeticError):
        pass

    return CertificateReplayResult(None, type(result).__name__, "unknown", {})


def result_diagnostics(result: object) -> ResultDiagnostics:
    """Return one stable diagnostic schema for CAD/QE/optimization certificates."""

    from .model import CADResult, QEResult

    if isinstance(result, CADResult):
        return ResultDiagnostics(
            "cad",
            result.projection_mode,
            bool(result.lifting_validated and result.well_oriented),
            dict(result.diagnostics.counters),
            {
                name: {"hits": item.hits, "misses": item.misses}
                for name, item in result.diagnostics.cache_stats.items()
            },
            {
                "nullification_events": len(result.nullification_events),
                "fallback_projection": result.used_fallback_projection,
                "lifting_certificates": len(result.lifting_certificates),
            },
        )
    if isinstance(result, QEResult):
        cad = result_diagnostics(result.cad)
        return ResultDiagnostics(
            "quantifier_elimination",
            cad.method,
            cad.certified,
            cad.counters,
            cad.cache_stats,
            {
                **dict(cad.details),
                "cells_visited": result.cells_visited,
                "partial_evaluation": result.partial_evaluation,
                "guided_pruning": result.guided_pruning,
            },
        )
    method = str(getattr(result, "method", type(result).__name__))
    certified = getattr(result, "certified", None)
    details = dict(getattr(result, "diagnostics", {}) or {})
    return ResultDiagnostics(type(result).__name__, method, certified, {}, {}, details)


__all__ = [
    "CertificateReplayResult",
    "ResultDiagnostics",
    "replay_certificate",
    "result_diagnostics",
]
