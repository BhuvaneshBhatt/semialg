"""Cooperative wall-clock timeouts for bounding symbolic operations."""

from __future__ import annotations

import signal
import threading
import time
from contextlib import contextmanager


class TimeoutError_(BaseException):
    """Internal timeout signal that ordinary ``Exception`` handlers do not hide."""


def _alarm_available() -> bool:
    """Return whether SIGALRM timers are usable in the current thread."""
    return (
        hasattr(signal, "SIGALRM")
        and hasattr(signal, "ITIMER_REAL")
        and hasattr(signal, "setitimer")
        and threading.current_thread() is threading.main_thread()
    )


@contextmanager
def bounded_time_context(seconds: float | None):
    """Bound a block with SIGALRM when the platform/thread supports it.

    On platforms without ``SIGALRM`` (notably Windows), and in non-main
    threads where Python forbids installing signal handlers, this context is a
    no-op. Nested POSIX timers preserve the original outer deadline.
    """
    if seconds is None or seconds <= 0 or not _alarm_available():
        yield
        return

    def _handler(signum, frame):
        del signum, frame
        raise TimeoutError_()

    old_handler = signal.signal(signal.SIGALRM, _handler)
    start = time.monotonic()
    old_remaining, old_interval = signal.setitimer(signal.ITIMER_REAL, seconds)
    try:
        yield
    finally:
        elapsed = time.monotonic() - start
        signal.signal(signal.SIGALRM, old_handler)
        if old_remaining > 0:
            restored = old_remaining - elapsed
            signal.setitimer(signal.ITIMER_REAL, max(restored, 1e-6), old_interval)
        else:
            signal.setitimer(signal.ITIMER_REAL, 0)


def run_with_time_budget(func, *args, seconds: float = 1.0, default=None, **kwargs):
    """Run ``func(*args, **kwargs)`` within a wall-clock budget.

    Returns ``default`` only when the supported wall-clock timer expires.
    Exceptions raised by ``func`` propagate so implementation bugs are not
    silently converted into inconclusive oracle results.
    """
    try:
        with bounded_time_context(seconds):
            return func(*args, **kwargs)
    except TimeoutError_:
        return default
