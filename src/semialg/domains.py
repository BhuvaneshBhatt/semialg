from __future__ import annotations

from collections.abc import Iterable, Mapping
from dataclasses import dataclass, field
from enum import Enum
from typing import Any

import sympy as sp


class SolveDomain(str, Enum):
    REALS = "reals"
    RATIONALS = "rationals"
    COMPLEXES = "complexes"
    INTEGERS = "integers"
    BOOLEANS = "booleans"


class UnsupportedDomainError(ValueError):
    """Raised when an API is asked to use a domain it cannot handle."""


@dataclass(frozen=True)
class SolveRequest:
    text: str
    variable_order: tuple[object, ...] = ()
    domain: SolveDomain = SolveDomain.REALS
    assumptions: tuple[sp.Expr, ...] = ()
    use_preprocess: bool = True


@dataclass(frozen=True)
class SemialgOptions:
    """Shared options accepted by high-level semialgebraic APIs."""

    domain: SolveDomain = SolveDomain.REALS
    assumptions: tuple[sp.Expr, ...] = ()
    strategy: str = "auto"
    output: str = "formula"
    count: int = 1
    max_cells: int | None = None
    timeout: float | None = None
    random_seed: int | None = None
    exact: bool = True
    diagnostics: bool = True
    strict: bool = False
    metadata: Mapping[str, Any] = field(default_factory=dict)

    @classmethod
    def from_values(
        cls,
        *,
        domain: str | SolveDomain | None = None,
        assumptions: Iterable[sp.Expr] | sp.Expr | None = None,
        strategy: str | None = None,
        output: str | None = None,
        count: int | None = None,
        max_cells: int | None = None,
        timeout: float | None = None,
        random_seed: int | None = None,
        exact: bool | None = None,
        diagnostics: bool | None = None,
        strict: bool = False,
        metadata: Mapping[str, Any] | None = None,
    ) -> SemialgOptions:
        return cls(
            domain=normalize_domain(domain),
            assumptions=normalize_assumptions(assumptions),
            strategy=strategy or "auto",
            output=output or "formula",
            count=1 if count is None else int(count),
            max_cells=max_cells,
            timeout=timeout,
            random_seed=random_seed,
            exact=True if exact is None else bool(exact),
            diagnostics=True if diagnostics is None else bool(diagnostics),
            strict=bool(strict),
            metadata=dict(metadata or {}),
        )


def normalize_domain(domain: str | SolveDomain | None) -> SolveDomain:
    if domain is None:
        return SolveDomain.REALS
    if isinstance(domain, SolveDomain):
        return domain
    key = str(domain).strip().lower()
    mapping = {
        "real": SolveDomain.REALS,
        "reals": SolveDomain.REALS,
        "r": SolveDomain.REALS,
        "rational": SolveDomain.RATIONALS,
        "rationals": SolveDomain.RATIONALS,
        "q": SolveDomain.RATIONALS,
        "complex": SolveDomain.COMPLEXES,
        "complexes": SolveDomain.COMPLEXES,
        "c": SolveDomain.COMPLEXES,
        "integer": SolveDomain.INTEGERS,
        "integers": SolveDomain.INTEGERS,
        "int": SolveDomain.INTEGERS,
        "z": SolveDomain.INTEGERS,
        "bool": SolveDomain.BOOLEANS,
        "bools": SolveDomain.BOOLEANS,
        "boolean": SolveDomain.BOOLEANS,
        "booleans": SolveDomain.BOOLEANS,
    }
    if key not in mapping:
        raise UnsupportedDomainError(f"unsupported solve domain: {domain}")
    return mapping[key]


def normalize_assumptions(assumptions: Iterable[sp.Expr] | sp.Expr | None) -> tuple[sp.Expr, ...]:
    if assumptions is None:
        return ()
    if isinstance(assumptions, (sp.Basic, sp.logic.boolalg.Boolean)):
        return (sp.sympify(assumptions),)
    return tuple(sp.sympify(item) for item in assumptions)


def apply_assumptions(
    expr: sp.Expr, assumptions: Iterable[sp.Expr] | sp.Expr | None = None
) -> sp.Expr:
    normalized = normalize_assumptions(assumptions)
    if not normalized:
        return expr
    return sp.And(*(normalized + (expr,)))


def domain_name(domain: str | SolveDomain | None) -> str:
    return normalize_domain(domain).value


__all__ = [
    "SemialgOptions",
    "SolveDomain",
    "SolveRequest",
    "UnsupportedDomainError",
    "apply_assumptions",
    "domain_name",
    "normalize_assumptions",
    "normalize_domain",
]
