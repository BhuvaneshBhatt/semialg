"""Time-bounded orchestration for independent real-instance methods."""

from __future__ import annotations

import signal
from collections.abc import Callable, Sequence

from .real_utils import _RECOVERABLE_ERRORS, MethodSearchResult


class _Timeout(Exception):
    """Internal signal used to abandon one timed method attempt."""


def _timeout_handler(signum, frame):  # type: ignore[no-untyped-def]
    raise _Timeout()


def try_methods(
    methods: Sequence[Callable[..., object]],
    args: Sequence[object] = (),
    *,
    initial_seconds: int = 1,
    growth_factor: int = 10,
    failure_value: object = None,
) -> MethodSearchResult:
    """Run methods with increasing time budgets.

    A method result that does not contain ``failure_value`` is considered a full
    success. Partial results containing other data plus ``failure_value`` are
    recorded and returned if no method fully succeeds.
    """

    remaining = list(methods)
    partial: list[object] = []
    budget = initial_seconds
    while remaining:
        timed_out: list[Callable[..., object]] = []
        for method in remaining:
            old_handler = signal.signal(signal.SIGALRM, _timeout_handler)
            signal.alarm(max(1, int(budget)))
            try:
                result = method(*args)
            except _Timeout:
                timed_out.append(method)
                continue
            except _RECOVERABLE_ERRORS as exc:
                result = failure_value
                partial.append(exc)
            finally:
                signal.alarm(0)
                signal.signal(signal.SIGALRM, old_handler)
            if result is failure_value:
                continue
            if isinstance(result, (list, tuple, set)) and failure_value in result:
                partial.append(result)
                continue
            return MethodSearchResult(
                result, "satisfied", tuple(partial), getattr(method, "__name__", repr(method))
            )
        remaining = timed_out
        budget = budget * growth_factor if len(remaining) > 1 else 10**9
    if partial:
        return MethodSearchResult(tuple(partial), "partial", tuple(partial))
    return MethodSearchResult(failure_value, "unknown")


__all__ = ["try_methods"]
