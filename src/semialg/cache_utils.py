from __future__ import annotations

from collections import OrderedDict
from collections.abc import Hashable
from threading import RLock
from typing import Generic, TypeVar

T = TypeVar("T")


class BoundedLRU(Generic[T]):
    """Small thread-safe process-local LRU with computation-context mirroring."""

    def __init__(self, maxsize: int, namespace: str) -> None:
        self.maxsize = int(maxsize)
        if self.maxsize < 1:
            raise ValueError("maxsize must be a positive integer")
        self.namespace = namespace
        self._data: OrderedDict[Hashable, T] = OrderedDict()
        self._lock = RLock()

    def get(self, key: Hashable) -> T | None:
        from .context import context_cache_get, context_cache_put

        found, value = context_cache_get(self.namespace, key)
        if found:
            return value
        with self._lock:
            value = self._data.get(key)
            if value is not None:
                self._data.move_to_end(key)
                context_cache_put(self.namespace, key, value)
            return value

    def put(self, key: Hashable, value: T) -> None:
        from .context import context_cache_put

        context_cache_put(self.namespace, key, value)
        with self._lock:
            self._data[key] = value
            self._data.move_to_end(key)
            while len(self._data) > self.maxsize:
                self._data.popitem(last=False)

    def clear(self) -> None:
        with self._lock:
            self._data.clear()

    def __len__(self) -> int:
        return len(self._data)


__all__ = ["BoundedLRU"]
