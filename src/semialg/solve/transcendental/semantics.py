from __future__ import annotations

from enum import Enum


class ResultSemantics(str, Enum):
    """Coverage relationship between a returned formula and the true solution set.

    ``EXACT`` means equivalence on the stated domain. ``SUBSET`` and ``SUPERSET``
    are global one-sided guarantees. ``WINDOW_*`` variants are guarantees only
    inside ``validity_window``. ``WINDOW_APPROXIMATION`` and
    ``PERIODIC_WINDOW_APPROX`` are heuristic/two-sided approximations and
    must not be treated as one-sided certificates. ``WINDOW_NO_WITNESS`` means
    only that the numerical search found no witness in the stated window.
    """

    UNKNOWN = "unknown"
    EXACT = "exact"
    SUBSET = "subset"
    SUPERSET = "superset"
    WITNESS_SUBSET = "witness_subset"
    WINDOW_SUBSET = "window_subset"
    WINDOW_SUPERSET = "window_superset"
    WINDOW_APPROXIMATION = "window_approximation"
    PERIODIC_WINDOW_APPROX = "periodic_window_approximation"
    WINDOW_NO_WITNESS = "window_no_witness"

    @property
    def is_exact(self) -> bool:
        return self is ResultSemantics.EXACT

    @property
    def is_subset(self) -> bool:
        return self in {
            ResultSemantics.SUBSET,
            ResultSemantics.WITNESS_SUBSET,
            ResultSemantics.WINDOW_SUBSET,
        }

    @property
    def is_superset(self) -> bool:
        return self in {ResultSemantics.SUPERSET, ResultSemantics.WINDOW_SUPERSET}

    @property
    def is_window_scoped(self) -> bool:
        return self in {
            ResultSemantics.WINDOW_SUBSET,
            ResultSemantics.WINDOW_SUPERSET,
            ResultSemantics.WINDOW_APPROXIMATION,
            ResultSemantics.PERIODIC_WINDOW_APPROX,
            ResultSemantics.WINDOW_NO_WITNESS,
        }


__all__ = ["ResultSemantics"]
