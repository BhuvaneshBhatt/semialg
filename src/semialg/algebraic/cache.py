from __future__ import annotations

from dataclasses import dataclass, field

import sympy as sp

from ..cache_utils import BoundedLRU

CACHE_VERSION = 5


@dataclass
class RootIsolationStats:
    calls: int = 0
    cache_hits: int = 0
    cache_misses: int = 0
    refinements: int = 0
    comparison_refinements: int = 0
    sign_hits: int = 0
    sign_misses: int = 0
    comparison_hits: int = 0
    comparison_misses: int = 0
    specialization_hits: int = 0
    specialization_misses: int = 0
    rur_hits: int = 0
    rur_misses: int = 0


@dataclass
class RootIsolationCache:
    """Bounded process-local caches for exact algebraic computations."""

    roots: BoundedLRU[tuple[object, ...]] = field(
        default_factory=lambda: BoundedLRU(1024, "algebraic.roots")
    )
    signs: BoundedLRU[int] = field(default_factory=lambda: BoundedLRU(4096, "algebraic.signs"))
    comparisons: BoundedLRU[int] = field(
        default_factory=lambda: BoundedLRU(4096, "algebraic.comparisons")
    )
    specializations: BoundedLRU[object] = field(
        default_factory=lambda: BoundedLRU(2048, "algebraic.specializations")
    )
    rur: BoundedLRU[object] = field(default_factory=lambda: BoundedLRU(128, "algebraic.rur"))
    stats: RootIsolationStats = field(default_factory=RootIsolationStats)

    def clear(self) -> None:
        self.roots.clear()
        self.signs.clear()
        self.comparisons.clear()
        self.specializations.clear()
        self.rur.clear()
        self.stats = RootIsolationStats()


CACHE = RootIsolationCache()


def expr_key(expr: sp.Expr) -> str:
    return sp.srepr(sp.factor(sp.expand(expr)))


def poly_key(poly: sp.Poly) -> str:
    return f"{tuple(str(g) for g in poly.gens)}::{expr_key(poly.as_expr())}"


def sample_expr_key(value: object) -> str:
    from .samples import AlgebraicRoot, RationalSample, sample_to_expr

    if isinstance(value, (RationalSample, AlgebraicRoot)):
        value = sample_to_expr(value)
    return sp.srepr(sp.sympify(value))


def root_isolation_costs() -> RootIsolationStats:
    return RootIsolationStats(**CACHE.stats.__dict__)


def algebraic_cache_stats() -> RootIsolationStats:
    """Return a snapshot of root/sign/comparison/specialization/RUR cache counters."""
    return root_isolation_costs()


def clear_algebraic_caches() -> None:
    """Clear all process-local exact-algebraic caches and counters."""
    CACHE.clear()


__all__ = [
    "CACHE_VERSION",
    "CACHE",
    "BoundedLRU",
    "RootIsolationStats",
    "algebraic_cache_stats",
    "clear_algebraic_caches",
    "expr_key",
    "poly_key",
    "root_isolation_costs",
    "sample_expr_key",
]
