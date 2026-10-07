"""Sound declarative composition rules for already-certified properties.

The module deliberately separates value algebra from interface validation.  Callers
that know domains/images should pass them; a rule is not returned when containment
cannot be proved.
"""

from __future__ import annotations

from dataclasses import dataclass

import sympy as sp

from .results import Convexity, Monotonicity, TruthValue


@dataclass(frozen=True)
class CompositionFact:
    property: str
    value: object
    theorem_id: str
    complete: bool = True
    bound: object = None


def _subset(a, b):
    if a is None or b is None:
        return None
    try:
        q = sp.sympify(a).is_subset(sp.sympify(b))
    except (TypeError, ValueError, NotImplementedError, AttributeError):
        return None
    return True if q is True else False if q is False else None


def _interface_ok(inner_image, outer_domain):
    return (
        True
        if inner_image is None and outer_domain is None
        else _subset(inner_image, outer_domain)
    )


def compose_property(
    prop,
    outer_value,
    inner_value,
    *,
    outer_increasing=True,
    outer_domain=None,
    inner_image=None,
    outer_complete=True,
    inner_complete=True,
    outer_bound=None,
    inner_bound=None,
):
    """Compose certified facts, refusing unresolved map-interface obligations."""
    if _interface_ok(inner_image, outer_domain) is not True and (
        inner_image is not None or outer_domain is not None
    ):
        return None
    complete = bool(outer_complete and inner_complete)
    if prop == "monotonicity":
        inc, dec = Monotonicity.STRICTLY_INCREASING, Monotonicity.STRICTLY_DECREASING
        table = {(inc, inc): inc, (inc, dec): dec, (dec, inc): dec, (dec, dec): inc}
        v = table.get((outer_value, inner_value))
        return (
            None
            if v is None
            else CompositionFact(prop, v, "compose.monotonicity", complete)
        )
    if prop in {"injective", "continuous", "differentiable", "analytic", "bijective"}:
        if outer_value is TruthValue.TRUE and inner_value is TruthValue.TRUE:
            return CompositionFact(prop, TruthValue.TRUE, f"compose.{prop}", complete)
        return None
    if prop == "convexity":
        if (
            outer_value is Convexity.CONVEX
            and inner_value is Convexity.CONVEX
            and outer_increasing
        ):
            return CompositionFact(
                prop, Convexity.CONVEX, "compose.convex-increasing-convex", complete
            )
        if (
            outer_value is Convexity.CONVEX
            and inner_value is Convexity.CONCAVE
            and not outer_increasing
        ):
            return CompositionFact(
                prop, Convexity.CONVEX, "compose.convex-decreasing-concave", complete
            )
        if (
            outer_value is Convexity.CONCAVE
            and inner_value is Convexity.CONCAVE
            and outer_increasing
        ):
            return CompositionFact(
                prop, Convexity.CONCAVE, "compose.concave-increasing-concave", complete
            )
        if (
            outer_value is Convexity.CONCAVE
            and inner_value is Convexity.CONVEX
            and not outer_increasing
        ):
            return CompositionFact(
                prop, Convexity.CONCAVE, "compose.concave-decreasing-convex", complete
            )
        return None
    if prop == "lipschitz":
        if (
            outer_value is TruthValue.TRUE
            and inner_value is TruthValue.TRUE
            and outer_bound is not None
            and inner_bound is not None
        ):
            b = sp.simplify(sp.sympify(outer_bound) * sp.sympify(inner_bound))
            if b.is_nonnegative is True:
                return CompositionFact(
                    prop, TruthValue.TRUE, "compose.lipschitz", complete, b
                )
        return None
    if prop == "bounded":
        if outer_value is TruthValue.TRUE and inner_value is TruthValue.TRUE:
            return CompositionFact(
                prop, TruthValue.TRUE, "compose.bounded-image", complete
            )
    return None
