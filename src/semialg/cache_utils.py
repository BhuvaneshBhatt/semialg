from __future__ import annotations

from collections import OrderedDict
from collections.abc import Hashable
from threading import RLock
from typing import Generic, TypeVar

T = TypeVar("T")
_MISSING = object()


class BoundedLRU(Generic[T]):
    """Small thread-safe process-local LRU with computation-context mirroring."""

    def __init__(self, maxsize: int, namespace: str) -> None:
        self.maxsize = int(maxsize)
        if self.maxsize < 1:
            raise ValueError("maxsize must be a positive integer")
        self.namespace = namespace
        self._data: OrderedDict[Hashable, T] = OrderedDict()
        self._lock = RLock()
        self._generation = 0

    def _context_namespace(self) -> tuple[str, int]:
        return (self.namespace, self._generation)

    def get(self, key: Hashable) -> T | None:
        from .context import context_cache_get, context_cache_put

        with self._lock:
            namespace = self._context_namespace()
            found, value = context_cache_get(namespace, key)
            if found:
                return value  # type: ignore[return-value]
            value = self._data.get(key, _MISSING)
            if value is _MISSING:
                return None
            self._data.move_to_end(key)
            context_cache_put(namespace, key, value)
            return value  # type: ignore[return-value]

    def put(self, key: Hashable, value: T) -> None:
        from .context import context_cache_put

        with self._lock:
            namespace = self._context_namespace()
            self._data[key] = value
            self._data.move_to_end(key)
            while len(self._data) > self.maxsize:
                self._data.popitem(last=False)
            context_cache_put(namespace, key, value)

    def clear(self) -> None:
        from .context import current_computation_context

        with self._lock:
            self._data.clear()
            self._generation += 1
            generation = self._generation
        context = current_computation_context()
        if context is not None:
            context.prune_namespace(self.namespace, generation)

    def resize(self, maxsize: int) -> None:
        """Change the capacity while preserving the newest cached entries."""

        new_size = int(maxsize)
        if new_size < 1:
            raise ValueError("maxsize must be a positive integer")
        with self._lock:
            self.maxsize = new_size
            while len(self._data) > self.maxsize:
                self._data.popitem(last=False)

    def __len__(self) -> int:
        with self._lock:
            return len(self._data)


__all__ = ["BoundedLRU"]
