from __future__ import annotations

from dataclasses import dataclass

from ..derived_geometry import is_compact
from ..geometry_queries import euler_characteristic
from ..normalization import normalize_formula, normalize_problem_variables
from ..regions.operations import region_dimension
from ..roadmaps import connected_component_count


@dataclass(frozen=True)
class TopologySummary:
    """Certified low-degree topological invariants of a semialgebraic set."""

    dimension: int
    connected_components: int
    euler_characteristic: int
    compact_support: bool
    betti_numbers: tuple[int | None, ...]

    @property
    def betti_0(self) -> int:
        return self.connected_components


def topology_summary(region, variables=None, *, compact_support: bool = True) -> TopologySummary:
    """Return exact component, Euler, and supported Betti-number information.

    ``b_0`` is available in every dimension. For compact sets of dimension at
    most one, ``b_1`` follows exactly from ``chi = b_0 - b_1``. Higher Betti
    numbers are left unspecified until an oriented cellular homology backend is
    available.
    """

    formula = normalize_formula(region)
    vars_ = normalize_problem_variables(variables, formula)
    dimension = int(region_dimension(formula, vars_))
    components = connected_component_count(formula, vars_)
    chi = euler_characteristic(formula, vars_, compact_support=compact_support)
    compact = is_compact(formula, vars_)
    betti: list[int | None] = [components]
    convex = False
    if compact and components == 1:
        try:
            from ..convexity import is_convex

            convex = is_convex(formula, vars_)
        except (ArithmeticError, TypeError, ValueError, NotImplementedError):
            convex = False
    if convex:
        betti.extend([0] * max(0, dimension))
    elif dimension >= 1:
        if compact and dimension == 1:
            betti.append(components - chi)
        else:
            betti.append(None)
        betti.extend([None] * max(0, dimension + 1 - len(betti)))
    return TopologySummary(dimension, components, chi, compact_support, tuple(betti))


def betti_number(region, degree: int, variables=None) -> int:
    """Return a certified Betti number when the current exact backend supports it."""

    if degree < 0:
        raise ValueError("Betti degree must be nonnegative")
    if degree == 0:
        return connected_component_count(region, variables)
    from ..standard_regions import StandardRegion

    if isinstance(region, StandardRegion):
        from .semialgebraic import triangulation_betti_numbers

        betti = triangulation_betti_numbers(region, variables)
        return 0 if degree >= len(betti) else int(betti[degree])
    if not is_compact(region, variables):
        raise NotImplementedError("positive-degree Betti numbers currently require a compact set")
    summary = topology_summary(region, variables, compact_support=True)
    if degree < len(summary.betti_numbers) and summary.betti_numbers[degree] is not None:
        return int(summary.betti_numbers[degree])
    # Supported certified triangulations give a complete finite chain complex.
    # Use exact rational boundary-matrix ranks before declaring the invariant
    # unsupported; this covers compact polyhedral regions in every dimension.
    from .semialgebraic import triangulation_betti_numbers

    try:
        betti = triangulation_betti_numbers(region, variables)
    except (ArithmeticError, TypeError, ValueError, NotImplementedError) as exc:
        raise NotImplementedError(
            "this Betti number requires a certified triangulation or a broader cellular-homology backend"
        ) from exc
    return 0 if degree >= len(betti) else int(betti[degree])


__all__ = ["TopologySummary", "topology_summary", "betti_number"]
