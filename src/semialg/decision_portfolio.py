"""Certified portfolio dispatch for specialized polynomial decisions."""

from __future__ import annotations

from collections.abc import Callable, Sequence
from dataclasses import dataclass
from typing import Any

import sympy as sp

from .sos_certificates import search_sos_certificate


@dataclass(frozen=True)
class CertifiedDecisionAttempt:
    backend: str
    applicable: bool
    complete: bool
    decision: bool | None
    witness: Any = None
    certificate: Any = None
    notes: tuple[str, ...] = ()


@dataclass(frozen=True)
class CertifiedDecisionResult:
    decision: bool
    backend: str
    attempts: tuple[CertifiedDecisionAttempt, ...]
    witness: Any = None
    certificate: Any = None


def _polynomial_nonnegative_portfolio(
    polynomial: sp.Expr,
    variables: Sequence[sp.Symbol],
    *,
    strategy: str = "auto",
    sos_backend: str | Callable[[sp.Expr, list[sp.Symbol]], object] = "auto",
    return_result: bool = False,
    random_lines: int = 8,
    seed: int = 1234,
) -> bool | CertifiedDecisionResult:
    """Decide global polynomial nonnegativity using certified fallbacks.

    Return the mathematical Boolean by default. Set ``return_result=True`` to
    obtain backend attempts, witnesses, and certificates.

    ``auto`` tries exact-verified SOS search, then Zeng (including its ARS
    positive-dimensional critical-locus path), then complete semialgebraic CAD.
    ``random_lines`` and ``seed`` tune only the Zeng witness-search stage.

    This is the internal dispatcher behind :func:`polynomial_nonnegative`; it is
    intentionally not a second public spelling of the same mathematical query.
    """

    from .polynomial_positivity import zeng_negative_point

    expr = sp.expand(sp.sympify(polynomial))
    vars_ = tuple(variables)
    strategy_name = strategy.lower()
    allowed = {"auto", "sos", "zeng", "ars", "cad"}
    if strategy_name not in allowed:
        raise ValueError(f"unknown polynomial decision strategy: {strategy!r}")
    if random_lines < 0:
        raise ValueError("random_lines must be nonnegative")
    attempts: list[CertifiedDecisionAttempt] = []

    def finish(result: CertifiedDecisionResult) -> bool | CertifiedDecisionResult:
        return result if return_result else result.decision

    if strategy_name in {"auto", "sos"}:
        sos = search_sos_certificate(
            expr, vars_, backend=sos_backend, use_planner=(strategy_name == "auto")
        )
        attempts.append(
            CertifiedDecisionAttempt(
                "sos",
                sos.plan.applicable if sos.plan is not None else True,
                sos.complete,
                True if sos.certified else None,
                certificate=sos.certificate,
                notes=sos.notes,
            )
        )
        if sos.certified:
            return finish(
                CertifiedDecisionResult(True, "sos", tuple(attempts), certificate=sos.certificate)
            )
        if strategy_name == "sos":
            raise NotImplementedError("SOS search did not yield an exactly verified certificate")

    if strategy_name in {"auto", "zeng", "ars"}:
        try:
            zeng = zeng_negative_point(expr, vars_, random_lines=random_lines, seed=seed)
        except (sp.PolynomialError, ValueError, TypeError):
            zeng = None
        if zeng is not None:
            decision = None if not zeng.complete else not bool(zeng.has_negative_point)
            attempts.append(
                CertifiedDecisionAttempt(
                    "zeng",
                    True,
                    zeng.complete,
                    decision,
                    witness=zeng.assignment,
                    certificate=zeng,
                    notes=zeng.notes,
                )
            )
            if decision is not None:
                return finish(
                    CertifiedDecisionResult(
                        decision,
                        zeng.method,
                        tuple(attempts),
                        witness=zeng.assignment,
                        certificate=zeng,
                    )
                )
        else:
            attempts.append(CertifiedDecisionAttempt("zeng", False, False, None))
        if strategy_name in {"zeng", "ars"}:
            raise NotImplementedError("Zeng/ARS polynomial decision is incomplete")

    # Complete fallback: f is nonnegative iff there is no real point f < 0.
    from .decision import is_satisfiable

    sat_result = is_satisfiable(expr < 0, vars_, strategy="cad", return_result=True)
    decision = not bool(sat_result.satisfiable)
    attempts.append(
        CertifiedDecisionAttempt(
            "cad",
            True,
            True,
            decision,
            witness=sat_result.witness,
            certificate=sat_result,
        )
    )
    return finish(
        CertifiedDecisionResult(
            decision, "cad", tuple(attempts), witness=sat_result.witness, certificate=sat_result
        )
    )


__all__ = [
    "CertifiedDecisionAttempt",
    "CertifiedDecisionResult",
]
