from __future__ import annotations

from collections.abc import Callable, Iterable, Mapping
from dataclasses import dataclass

import sympy as sp

from .solution_checking import form_sat_by_assign


@dataclass(frozen=True)
class BackendCrosscheckResult:
    backend_name: str
    consistent: bool
    checked_assignments: tuple[dict[sp.Symbol, object], ...]
    failing_assignments: tuple[dict[sp.Symbol, object], ...]


def crosscheck_backend_pred(
    formula: sp.Expr,
    backend_name: str,
    predicate: Callable[[Mapping[sp.Symbol, object]], bool],
    assignments: Iterable[Mapping[sp.Symbol, object]],
    *,
    domain=sp.Complexes,
    modulus: int | None = None,
) -> BackendCrosscheckResult:
    checked = []
    failing = []
    for assignment in assignments:
        assignment = dict(assignment)
        checked.append(assignment)
        truth = form_sat_by_assign(
            formula, assignment, domain=domain, modulus=modulus, check_numeric_equalities=True
        )
        if bool(predicate(assignment)) != truth:
            failing.append(assignment)
    return BackendCrosscheckResult(backend_name, not failing, tuple(checked), tuple(failing))


__all__ = ["BackendCrosscheckResult", "crosscheck_backend_pred"]
