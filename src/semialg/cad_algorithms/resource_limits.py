"""Cooperative resource limits for exact CAD projection and lifting."""

import math
from contextlib import contextmanager
from contextvars import ContextVar
from time import perf_counter

from ..errors import ResourceLimitError


class CADBudget:
    """Track total constructed cells and elapsed time between exact operations.

    One symbolic projection, root-isolation or stack operation may overrun a
    deadline. No partially constructed CAD is returned as a complete result.
    """

    def __init__(self, max_cells, timeout):
        if max_cells is not None and (
            isinstance(max_cells, bool) or not isinstance(max_cells, int) or max_cells < 1
        ):
            raise ValueError("max_cells must be a positive integer or None")
        if timeout is not None and (
            isinstance(timeout, bool) or not math.isfinite(timeout) or timeout <= 0
        ):
            raise ValueError("timeout must be a positive finite number or None")
        self.max_cells = max_cells
        self._previous_cells = 0
        self._attempt_peak = 0
        self.deadline = None if timeout is None else perf_counter() + timeout

    def check(self, count=0):
        """Raise before continuing after a count or cooperative deadline expires."""
        if self.deadline is not None and perf_counter() >= self.deadline:
            raise ResourceLimitError("CAD deadline exceeded")
        self._attempt_peak = max(self._attempt_peak, count)
        if self.max_cells is not None and self._previous_cells + count > self.max_cells:
            raise ResourceLimitError("CAD total cell limit exceeded")

    def complete(self, count):
        """Charge a successful construction or cache return before the next one."""
        self.check(count)
        self._previous_cells += self._attempt_peak
        self._attempt_peak = 0

    def restart(self):
        """Retain work charged to a failed lift before starting its replacement."""
        self._previous_cells += self._attempt_peak
        self._attempt_peak = 0
        self.check()


# A high-level operation can construct several decompositions (for example,
# algebraic graph elimination followed by the output CAD). Context-local
# ownership keeps nested solvers on one deadline without process-global state.
_ACTIVE_BUDGET = ContextVar("semialg_cad_budget", default=None)


def current_budget():
    """Return the budget owned by the current high-level operation, if any."""
    return _ACTIVE_BUDGET.get()


@contextmanager
def budget_scope(budget):
    """Share a cooperative budget across nested exact CAD computations."""
    token = _ACTIVE_BUDGET.set(budget)
    try:
        yield budget
    finally:
        _ACTIVE_BUDGET.reset(token)
