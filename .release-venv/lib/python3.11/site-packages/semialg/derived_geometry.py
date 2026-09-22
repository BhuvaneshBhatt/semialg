"""Derived exact semialgebraic set and geometry operations.

These functions compose the package's existing certified CAD/QE,
optimization, topology, and integration primitives instead of introducing
parallel algorithms.
"""

from __future__ import annotations

import sympy as sp

from ._zero_testing import certified_zero
from .connectivity import extract_cad_connectivity, factorized_equality_components
from .decision import is_satisfiable
from .geometry_queries import (
    distance_between_regions,
    distance_to_region,
    is_path_connected,
)
from .internal_symbols import fresh_real_dummy
from .moments import region_moment
from .normalization import normalize_formula, normalize_point, normalize_problem_variables
from .optimization import function_range, semialgebraic_maximize, semialgebraic_minimize
from .optimization_results import FunctionRangeResult, OptimizationResult
from .reasoning import (
    region_bounded,
    region_closed,
    region_compact,
    region_disjoint,
    region_equal,
    region_subset,
)
from .regions.operations import (
    explicit_region_components,
    region_closure,
    region_dimension,
    region_interior,
)


def _optimization_result(value: object, context: str) -> OptimizationResult:
    if not isinstance(value, OptimizationResult):
        raise TypeError(f"{context} returned a non-OptimizationResult")
    return value


def argmin_set(expression, region=sp.true, variables=None) -> sp.Expr:
    """Return the exact global minimizer set as a semialgebraic formula."""
    expr = sp.sympify(expression)
    formula = normalize_formula(region)
    vars_ = normalize_problem_variables(variables, sp.Tuple(expr, formula))
    result = _optimization_result(
        semialgebraic_minimize(expr, formula, vars_, return_result=True), "minimization"
    )
    if not result.attained:
        return sp.false
    return sp.And(formula, sp.Eq(expr, result.value))


def argmax_set(expression, region=sp.true, variables=None) -> sp.Expr:
    """Return the exact global maximizer set as a semialgebraic formula."""
    expr = sp.sympify(expression)
    formula = normalize_formula(region)
    vars_ = normalize_problem_variables(variables, sp.Tuple(expr, formula))
    result = _optimization_result(
        semialgebraic_maximize(expr, formula, vars_, return_result=True), "maximization"
    )
    if not result.attained:
        return sp.false
    return sp.And(formula, sp.Eq(expr, result.value))


def extrema_set(expression, region=sp.true, variables=None) -> sp.Expr:
    """Return the union of the exact global minimum and maximum sets."""
    return sp.Or(
        argmin_set(expression, region, variables), argmax_set(expression, region, variables)
    )


def level_set(expression, value=0, region=sp.true) -> sp.Expr:
    """Return ``region ∩ {expression = value}``."""
    return sp.And(normalize_formula(region), sp.Eq(sp.sympify(expression), sp.sympify(value)))


def sublevel_set(expression, value=0, region=sp.true, *, strict: bool = False) -> sp.Expr:
    """Return ``region ∩ {expression <= value}`` (or strict variant)."""
    expr, val = sp.sympify(expression), sp.sympify(value)
    rel = expr < val if strict else expr <= val
    return sp.And(normalize_formula(region), rel)


def superlevel_set(expression, value=0, region=sp.true, *, strict: bool = False) -> sp.Expr:
    """Return ``region ∩ {expression >= value}`` (or strict variant)."""
    expr, val = sp.sympify(expression), sp.sympify(value)
    rel = expr > val if strict else expr >= val
    return sp.And(normalize_formula(region), rel)


def is_empty(region, variables=None, *, strategy: str | None = None) -> bool:
    """Return whether the semialgebraic region is empty."""
    formula = normalize_formula(region)
    vars_ = normalize_problem_variables(variables, formula)
    return not is_satisfiable(formula, vars_, strategy=strategy)


def is_bounded(region, variables=None, *, strategy: str | None = None) -> bool:
    """Return whether the semialgebraic region is bounded."""
    return region_bounded(region, variables, strategy=strategy)


def is_compact(region, variables=None, *, strategy: str | None = None) -> bool:
    """Return whether the semialgebraic region is compact."""
    return region_compact(region, variables, strategy=strategy)


def is_closed(region, variables=None, *, strategy: str | None = None) -> bool:
    """Return whether the semialgebraic region is closed."""
    return region_closed(region, variables, strategy=strategy)


def is_open(region, variables=None, *, strategy: str | None = None) -> bool:
    """Return whether the semialgebraic region is open."""
    formula = normalize_formula(region)
    vars_ = normalize_problem_variables(variables, formula)
    return region_equal(formula, region_interior(formula, vars_), vars_, strategy=strategy)


def is_subset(left, right, variables=None, *, strategy: str | None = None) -> bool:
    """Return whether one semialgebraic region is contained in another."""
    return region_subset(left, right, variables, strategy=strategy)


def is_equal(left, right, variables=None, *, strategy: str | None = None) -> bool:
    """Return whether two semialgebraic regions define the same set."""
    return region_equal(left, right, variables, strategy=strategy)


def is_disjoint(left, right, variables=None, *, strategy: str | None = None) -> bool:
    """Return whether two semialgebraic regions are disjoint."""
    return region_disjoint(left, right, variables, strategy=strategy)


def intersects(left, right, variables=None, *, strategy: str | None = None) -> bool:
    """Return whether two semialgebraic regions have nonempty intersection."""
    return not is_disjoint(left, right, variables, strategy=strategy)


def is_interior_disjoint(left, right, variables=None, *, strategy: str | None = None) -> bool:
    """Return whether two regions have disjoint ambient interiors.

    Boundary contact is permitted: two closed regions that meet only along
    their boundaries are interior-disjoint.
    """
    left_interior = region_interior(left, variables, strategy=strategy)
    right_interior = region_interior(right, variables, strategy=strategy)
    return region_disjoint(left_interior, right_interior, variables, strategy=strategy)


def is_dense_in(subset, ambient, variables=None, *, strategy: str | None = None) -> bool:
    """Return whether ``subset`` is dense in ``ambient`` in the ambient Euclidean topology."""
    left = normalize_formula(subset)
    right = normalize_formula(ambient)
    vars_ = normalize_problem_variables(variables, sp.And(left, right))
    return region_subset(right, region_closure(left, vars_), vars_, strategy=strategy)


def contains_point(region, point, variables=None) -> bool:
    """Return whether an exact point belongs to the semialgebraic region."""
    formula = normalize_formula(region)
    vars_ = normalize_problem_variables(variables, formula)
    point_map = normalize_point(point, vars_, context=(formula,))
    value = sp.simplify(formula.subs(point_map))
    if value in (sp.true, True):
        return True
    if value in (sp.false, False):
        return False
    raise ValueError("point membership did not reduce to an exact Boolean value")


def connected_components(region, variables=None) -> tuple[sp.Expr, ...]:
    """Return exact CAD-connected-component formulas."""
    formula = normalize_formula(region)
    vars_ = normalize_problem_variables(variables, formula)
    factorized = factorized_equality_components(formula, vars_)
    if factorized is not None:
        return factorized
    if len(vars_) == 1:
        # The one-dimensional reducer returns exact maximal interval/point
        # components, so use it even when the set is connected. Rebuilding a
        # single component from CAD cells needlessly fragments the formula.
        return explicit_region_components(formula, vars_)
    graph = extract_cad_connectivity(formula, vars_)
    return tuple(component.as_formula(closed=False) for component in graph.components)


def is_connected(region, variables=None) -> bool:
    """Decide connectedness; for semialgebraic sets this equals path connectedness."""
    return is_path_connected(region, variables)


def is_full_dimensional(region, variables=None) -> bool:
    """Return whether the region has full dimension in its ambient variables."""
    formula = normalize_formula(region)
    vars_ = normalize_problem_variables(variables, formula)
    return region_dimension(formula, vars_) == len(vars_)


def has_empty_interior(region, variables=None) -> bool:
    """Return whether the region has empty ambient interior."""
    formula = normalize_formula(region)
    vars_ = normalize_problem_variables(variables, formula)
    return region_dimension(formula, vars_) < len(vars_)


def coordinate_range(region, variable, variables=None, *, return_result: bool = False):
    """Return the exact range of one coordinate over a region."""
    formula = normalize_formula(region)
    vars_ = normalize_problem_variables(variables, formula)
    # Resolve through the normalized variables without weakening assumptions.
    matches = [
        v for v in vars_ if variable == v or (isinstance(variable, str) and v.name == variable)
    ]
    if len(matches) != 1:
        raise ValueError("coordinate must resolve uniquely to a problem variable")
    return function_range(matches[0], formula, vars_, return_result=return_result)


def nearest_point(point, region, variables=None):
    """Return all exact nearest points when the distance is attained."""
    result = distance_to_region(point, region, variables, return_result=True)
    return result.nearest_points


def closest_points(left, right, variables=None):
    """Return all exact closest point pairs when the distance is attained."""
    result = distance_between_regions(left, right, variables, return_result=True)
    return result.closest_pairs


def diameter(region, variables=None, *, return_result: bool = False):
    """Return the exact Euclidean diameter (supremal pairwise distance)."""
    formula = normalize_formula(region)
    vars_ = normalize_problem_variables(variables, formula)
    other = tuple(fresh_real_dummy(f"{v.name}_diam") for v in vars_)
    copied = formula.xreplace(dict(zip(vars_, other, strict=True)))
    squared = sp.Add(*((x - y) ** 2 for x, y in zip(vars_, other, strict=True)))
    opt = _optimization_result(
        semialgebraic_maximize(
            squared, sp.And(formula, copied), (*vars_, *other), return_result=True
        ),
        "diameter optimization",
    )
    if return_result:
        return opt
    return sp.sqrt(opt.value)


def support_function(region, direction, variables=None, *, return_result: bool = False):
    """Return ``sup(x·direction)`` over the region."""
    formula = normalize_formula(region)
    vars_ = normalize_problem_variables(variables, formula)
    direction = tuple(map(sp.sympify, direction))
    if len(direction) != len(vars_):
        raise ValueError("direction must have the same dimension as variables")
    objective = sp.Add(*(u * x for u, x in zip(direction, vars_, strict=True)))
    result = _optimization_result(
        semialgebraic_maximize(objective, formula, vars_, return_result=True),
        "support optimization",
    )
    return result if return_result else result.value


def width(region, direction, variables=None) -> sp.Expr:
    """Return exact directional width ``max u·x - min u·x``."""
    formula = normalize_formula(region)
    vars_ = normalize_problem_variables(variables, formula)
    direction = tuple(map(sp.sympify, direction))
    if len(direction) != len(vars_):
        raise ValueError("direction must have the same dimension as variables")
    objective = sp.Add(*(u * x for u, x in zip(direction, vars_, strict=True)))
    hi = _optimization_result(
        semialgebraic_maximize(objective, formula, vars_, return_result=True), "width maximum"
    )
    lo = _optimization_result(
        semialgebraic_minimize(objective, formula, vars_, return_result=True), "width minimum"
    )
    return sp.simplify(hi.value - lo.value)


def moment_matrix(region, variables, **kwargs) -> sp.Matrix:
    """Return the normalized raw second-moment matrix ``E[x x.T]``."""
    formula = normalize_formula(region)
    vars_ = normalize_problem_variables(variables, formula)
    measure = region_moment(formula, vars_, tuple(0 for _ in vars_), **kwargs)
    if certified_zero(measure) is True:
        raise ValueError("region has zero ambient measure")
    entries = []
    for x in vars_:
        row = []
        for y in vars_:
            raw = region_moment(formula, vars_, integrand=x * y, **kwargs)
            row.append(sp.simplify(raw / measure))
        entries.append(row)
    return sp.Matrix(entries)


def inertia_tensor(region, variables, **kwargs) -> sp.Matrix:
    """Return the unit-density second moment-of-inertia tensor about the origin."""
    formula = normalize_formula(region)
    vars_ = normalize_problem_variables(variables, formula)
    norm_sq = sp.Add(*(v**2 for v in vars_))
    entries = []
    for i, x in enumerate(vars_):
        row = []
        for j, y in enumerate(vars_):
            integrand = (norm_sq if i == j else 0) - x * y
            row.append(region_moment(formula, vars_, integrand=integrand, **kwargs))
        entries.append(row)
    return sp.Matrix(entries)


def translate(region, vector, variables=None) -> sp.Expr:
    """Translate a region by ``vector`` while preserving coordinate symbols."""
    formula = normalize_formula(region)
    vars_ = normalize_problem_variables(variables, formula)
    vector = tuple(map(sp.sympify, vector))
    if len(vector) != len(vars_):
        raise ValueError("translation vector must have the same dimension as variables")
    return sp.simplify(formula.xreplace({v: v - a for v, a in zip(vars_, vector, strict=True)}))


def scale(region, factor, variables=None) -> sp.Expr:
    """Scale a region about the origin by a scalar factor."""
    formula = normalize_formula(region)
    vars_ = normalize_problem_variables(variables, formula)
    factor = sp.sympify(factor)
    if factor == 0:
        if is_empty(formula, vars_):
            return sp.false
        return sp.And(*(sp.Eq(v, 0) for v in vars_))
    return sp.simplify(formula.xreplace({v: v / factor for v in vars_}))


def linear_image(region, matrix, variables=None):
    """Return the exact linear image ``A*x``.

    This is the zero-offset specialization of :func:`semialg.affine_image` and
    uses the same representation-dispatching contract.
    """
    from .region_transformations import affine_image

    return affine_image(region, matrix, None, variables)


def squared_distance_range(left, right=None, variables=None, *, return_result: bool = False):
    """Return the exact range of squared pairwise distances.

    When ``right`` is omitted, both points range over ``left``.
    """
    left_formula = normalize_formula(left)
    right_formula = left_formula if right is None else normalize_formula(right)
    vars_ = normalize_problem_variables(variables, sp.And(left_formula, right_formula))
    rv = tuple(fresh_real_dummy(f"{v.name}_dist_r") for v in vars_)
    right_copy = right_formula.xreplace(dict(zip(vars_, rv, strict=True)))
    squared = sp.Add(*((x - y) ** 2 for x, y in zip(vars_, rv, strict=True)))
    return function_range(
        squared, sp.And(left_formula, right_copy), (*vars_, *rv), return_result=return_result
    )


def distance_set(left, right=None, variables=None, *, distance_symbol=None) -> sp.Expr:
    """Return the exact set of Euclidean pairwise distances as a formula."""
    d = sp.Symbol("d", real=True) if distance_symbol is None else sp.sympify(distance_symbol)
    if not isinstance(d, sp.Symbol):
        raise TypeError("distance_symbol must be a Symbol")
    result = squared_distance_range(left, right, variables, return_result=True)
    if not isinstance(result, FunctionRangeResult):
        raise TypeError("squared-distance range returned an unexpected result type")
    return sp.And(d >= 0, result.formula.xreplace({result.value_symbol: d**2}))


def minkowski_sum(left, right, variables=None):
    """Return the exact Minkowski sum of two semialgebraic regions.

    Canonical points, intervals, boxes, and filled balls retain structural
    representations.  General cases use exact semialgebraic image elimination.
    """
    from .geometry_queries import semialgebraic_image
    from .region_structural import UNRESOLVED, structural_minkowski_sum
    from .standard_regions import StandardRegion
    from .symbolic_regions import SemialgebraicRegion

    structural = structural_minkowski_sum(left, right)
    if structural is not UNRESOLVED:
        return structural

    object_input = isinstance(left, (StandardRegion, SemialgebraicRegion)) or isinstance(
        right, (StandardRegion, SemialgebraicRegion)
    )
    left_formula = normalize_formula(left)
    right_formula = normalize_formula(right)
    vars_ = normalize_problem_variables(variables, sp.And(left_formula, right_formula))
    rv = tuple(fresh_real_dummy(f"{v.name}_mink_r") for v in vars_)
    right_copy = right_formula.xreplace(dict(zip(vars_, rv, strict=True)))
    mapping = tuple(x + y for x, y in zip(vars_, rv, strict=True))
    targets = tuple(fresh_real_dummy(f"{v.name}_mink_t") for v in vars_)
    image = semialgebraic_image(
        mapping,
        sp.And(left_formula, right_copy),
        (*vars_, *rv),
        image_variables=targets,
    )
    formula = sp.simplify(image.xreplace(dict(zip(targets, vars_, strict=True))))
    return SemialgebraicRegion(formula, vars_) if object_input else formula


__all__ = [
    "argmin_set",
    "argmax_set",
    "extrema_set",
    "level_set",
    "sublevel_set",
    "superlevel_set",
    "is_empty",
    "is_bounded",
    "is_compact",
    "is_open",
    "is_closed",
    "is_subset",
    "is_equal",
    "is_disjoint",
    "intersects",
    "is_interior_disjoint",
    "is_dense_in",
    "contains_point",
    "connected_components",
    "is_connected",
    "is_full_dimensional",
    "has_empty_interior",
    "coordinate_range",
    "nearest_point",
    "closest_points",
    "diameter",
    "support_function",
    "width",
    "moment_matrix",
    "inertia_tensor",
    "translate",
    "scale",
    "linear_image",
    "minkowski_sum",
    "squared_distance_range",
    "distance_set",
]
