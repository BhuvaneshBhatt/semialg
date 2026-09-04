from __future__ import annotations

from dataclasses import dataclass, field

import sympy as sp

from ..cache_utils import BoundedLRU

CACHE_FORMAT = 6


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


@dataclass(frozen=True)
class AlgebraicCacheLimits:
    """Capacities for process-local exact-algebraic performance caches."""

    roots: int = 1024
    signs: int = 4096
    comparisons: int = 4096
    specializations: int = 2048
    rur: int = 128


def configure_algebraic_cache_limits(
    *,
    roots: int | None = None,
    signs: int | None = None,
    comparisons: int | None = None,
    specializations: int | None = None,
    rur: int | None = None,
) -> AlgebraicCacheLimits:
    """Resize process-local algebraic caches and return the resulting limits.

    Existing newest entries are retained up to each new capacity.  This is an
    expert performance-control API; it does not change mathematical results.
    """

    updates = {
        "roots": roots,
        "signs": signs,
        "comparisons": comparisons,
        "specializations": specializations,
        "rur": rur,
    }
    for name, value in updates.items():
        if value is not None:
            getattr(CACHE, name).resize(value)
    return AlgebraicCacheLimits(
        roots=CACHE.roots.maxsize,
        signs=CACHE.signs.maxsize,
        comparisons=CACHE.comparisons.maxsize,
        specializations=CACHE.specializations.maxsize,
        rur=CACHE.rur.maxsize,
    )


def expr_key(expr: sp.Expr) -> sp.Expr:
    """Return a cheap structural expression key without algebraic rewriting.

    Cache lookup must not trigger factorization, rational-function combination,
    or other symbolic normalization. Callers may therefore miss reuse when two
    algebraically equivalent expressions have different trees, but a cache miss
    is preferable to making key construction asymptotically more expensive than
    the exact operation being cached.
    """

    return sp.sympify(expr)


def poly_key(
    poly: sp.Poly,
) -> tuple[tuple[sp.Symbol, ...], object, tuple[tuple[tuple[int, ...], sp.Expr], ...]]:
    """Return an exact structural polynomial key preserving symbol identity."""

    return (tuple(poly.gens), poly.domain, tuple(poly.terms()))


def sample_expr_key(value: object) -> sp.Expr:
    from .samples import AlgebraicRoot, RationalSample, sample_to_expr

    if isinstance(value, (RationalSample, AlgebraicRoot)):
        value = sample_to_expr(value)
    return sp.sympify(value)


def root_isolation_costs() -> RootIsolationStats:
    return RootIsolationStats(**CACHE.stats.__dict__)


def algebraic_cache_stats() -> RootIsolationStats:
    """Return a snapshot of root/sign/comparison/specialization/RUR cache counters."""
    return root_isolation_costs()


def clear_algebraic_caches() -> None:
    """Clear all process-local exact-algebraic caches and counters."""
    CACHE.clear()


__all__ = [
    "CACHE_FORMAT",
    "CACHE",
    "AlgebraicCacheLimits",
    "BoundedLRU",
    "RootIsolationStats",
    "algebraic_cache_stats",
    "clear_algebraic_caches",
    "configure_algebraic_cache_limits",
    "expr_key",
    "poly_key",
    "root_isolation_costs",
    "sample_expr_key",
]
