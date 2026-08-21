from __future__ import annotations

from collections.abc import Hashable, Iterator
from contextlib import contextmanager
from contextvars import ContextVar
from dataclasses import dataclass, field
from functools import wraps
from typing import Any


@dataclass
class ExactComputationContext:
    """Per-operation cache shared by exact CAD/algebraic subroutines.

    Context entries live only for one top-level solve/optimization/range call.
    Process-local bounded LRUs remain a second-level cache, but repeated work
    inside a solve is resolved here first and disappears when the operation
    finishes.  This avoids global cache pollution while allowing nested solver
    calls to reuse projection, root, sign, comparison, specialization and RUR
    results.
    """

    caches: dict[str, dict[Hashable, Any]] = field(default_factory=dict)
    hits: dict[str, int] = field(default_factory=dict)
    misses: dict[str, int] = field(default_factory=dict)

    def get(self, namespace: str, key: Hashable) -> tuple[bool, Any]:
        bucket = self.caches.get(namespace)
        if bucket is not None and key in bucket:
            self.hits[namespace] = self.hits.get(namespace, 0) + 1
            return True, bucket[key]
        self.misses[namespace] = self.misses.get(namespace, 0) + 1
        return False, None

    def put(self, namespace: str, key: Hashable, value: Any) -> None:
        self.caches.setdefault(namespace, {})[key] = value

    def cache_size(self, namespace: str | None = None) -> int:
        if namespace is not None:
            return len(self.caches.get(namespace, {}))
        return sum(len(bucket) for bucket in self.caches.values())

    def stats(self) -> dict[str, dict[str, int]]:
        names = set(self.caches) | set(self.hits) | set(self.misses)
        return {
            name: {
                "hits": self.hits.get(name, 0),
                "misses": self.misses.get(name, 0),
                "size": len(self.caches.get(name, {})),
            }
            for name in sorted(names)
        }


_CURRENT_CONTEXT: ContextVar[ExactComputationContext | None] = ContextVar(
    "semialg_exact_computation_context", default=None
)


def current_computation_context() -> ExactComputationContext | None:
    """Return the active exact-computation context, if any."""

    return _CURRENT_CONTEXT.get()


@contextmanager
def computation_context(
    context: ExactComputationContext | None = None,
) -> Iterator[ExactComputationContext]:
    """Create/reuse a context for one complete exact solve operation.

    Nested calls automatically reuse the active context.  Passing an explicit
    context installs it only when no context is already active.
    """

    active = _CURRENT_CONTEXT.get()
    if active is not None:
        yield active
        return
    owned = context if context is not None else ExactComputationContext()
    token = _CURRENT_CONTEXT.set(owned)
    try:
        yield owned
    finally:
        _CURRENT_CONTEXT.reset(token)


def with_computation_context(function):
    """Decorator that gives a top-level solver call one reusable cache scope."""

    @wraps(function)
    def wrapped(*args, **kwargs):
        with computation_context():
            return function(*args, **kwargs)

    return wrapped


def context_cache_get(namespace: str, key: Hashable) -> tuple[bool, Any]:
    context = _CURRENT_CONTEXT.get()
    if context is None:
        return False, None
    return context.get(namespace, key)


def context_cache_put(namespace: str, key: Hashable, value: Any) -> None:
    context = _CURRENT_CONTEXT.get()
    if context is not None:
        context.put(namespace, key, value)


__all__ = [
    "ExactComputationContext",
    "computation_context",
    "current_computation_context",
    "with_computation_context",
]
