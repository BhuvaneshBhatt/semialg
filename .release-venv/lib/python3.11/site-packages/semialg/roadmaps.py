from __future__ import annotations

from dataclasses import dataclass

import sympy as sp

from .connectivity import CADConnectivityGraph, extract_cad_connectivity
from .normalization import normalize_formula, normalize_problem_variables
from .regions.operations import region_dimension


@dataclass(frozen=True)
class SemialgebraicRoadmap:
    """Certified one-dimensional roadmap of a semialgebraic set.

    Sets of dimension at most one use the set itself. Compact convex sets use
    a projection-spanning segment whose inclusion and fiber-intersection
    properties follow from convexity. General higher-dimensional roadmap
    construction remains a separate critical-point algorithm.
    """

    variables: tuple[sp.Symbol, ...]
    region_formula: sp.Expr
    roadmap_formula: sp.Expr
    dimension: int
    connectivity: CADConnectivityGraph
    construction: str = "identity-low-dimensional"
    rm1_certified: bool = True
    rm2_certified: bool = True

    @property
    def component_count(self) -> int:
        return self.connectivity.component_count

    @property
    def sample_points(self):
        return self.connectivity.sample_points

    @property
    def components(self):
        return self.connectivity.components


def roadmap(region, variables=None) -> SemialgebraicRoadmap:
    """Construct an exact Basu--Pollack--Roy roadmap when currently supported.

    Zero- and one-dimensional semialgebraic sets are returned as their own
    roadmap. Compact convex sets use an exact projection-spanning segment.
    Other higher-dimensional sets require the pseudo-critical-value recursion
    and are not approximated by a CAD-cell adjacency graph.
    """

    formula = normalize_formula(region)
    vars_ = normalize_problem_variables(variables, formula)
    dimension = int(region_dimension(formula, vars_))
    graph = extract_cad_connectivity(formula, vars_)
    if dimension <= 1:
        return SemialgebraicRoadmap(
            variables=vars_,
            region_formula=formula,
            roadmap_formula=formula,
            dimension=dimension,
            connectivity=graph,
        )

    # A compact convex set admits a particularly cheap exact roadmap. Every
    # fiber of the first-coordinate projection is convex and hence connected.
    # A segment joining points attaining the minimum and maximum first
    # coordinate stays inside the set and meets every nonempty fiber.
    from .convexity import is_convex
    from .derived_geometry import is_compact
    from .optimization import semialgebraic_maximize, semialgebraic_minimize
    from .standard_regions import Simplex

    if is_compact(formula, vars_) and is_convex(formula, vars_):
        coordinate = vars_[0]
        minimum = semialgebraic_minimize(coordinate, formula, vars_, return_result=True)
        maximum = semialgebraic_maximize(coordinate, formula, vars_, return_result=True)
        if minimum.attained and maximum.attained and minimum.point and maximum.point:
            left = tuple(minimum.point[var] for var in vars_)
            right = tuple(maximum.point[var] for var in vars_)
            vertices = (left,) if minimum.value == maximum.value else (left, right)
            section = Simplex(vertices)
            return SemialgebraicRoadmap(
                variables=vars_,
                region_formula=formula,
                roadmap_formula=section.as_formula(vars_, eliminate=True),
                dimension=section.dimension(),
                connectivity=graph,
                construction="convex-projection-section",
            )

    raise NotImplementedError(
        "general higher-dimensional roadmaps require the recursive pseudo-critical-value construction"
    )


def connected_component_count(region, variables=None) -> int:
    """Return the exact number of semialgebraically connected components."""

    formula = normalize_formula(region)
    vars_ = normalize_problem_variables(variables, formula)
    return extract_cad_connectivity(formula, vars_).component_count


def connected_component_samples(region, variables=None):
    """Return one exact CAD sample point from every connected component."""

    formula = normalize_formula(region)
    vars_ = normalize_problem_variables(variables, formula)
    return extract_cad_connectivity(formula, vars_).sample_points


__all__ = [
    "SemialgebraicRoadmap",
    "roadmap",
    "connected_component_count",
    "connected_component_samples",
]
