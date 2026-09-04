"""High-level exact semialgebraic geometry queries built from CAD/QE and optimization."""

from __future__ import annotations

from collections import deque
from collections.abc import Mapping, Sequence
from dataclasses import dataclass

import sympy as sp

from ._critical_loci import iter_positive_dimensional_critical_loci
from .cad_algorithms.cells import extract_cylindrical_solution
from .connectivity import (
    CADConnectivityGraph,
    extract_cad_connectivity,
    factorized_equality_components,
)
from .convexity import (
    _is_affine_polyhedron,  # noqa: F401 - shared implementation used by focused tests
    _is_convex_by_definition,  # noqa: F401 - shared implementation used by focused tests
    _is_convex_polynomial_intersection,  # noqa: F401 - shared implementation used by focused tests
    is_convex,
)
from .errors import ResourceLimitError
from .formula import parse_formula
from .formulas.boolean import is_true_expr
from .internal_symbols import fresh_real_dummy
from .normalization import (
    normalize_formula,
    normalize_point,
    normalize_problem_variables,
    normalize_variables,
)
from .optimization import (
    _constraint_data,
    _kkt_candidate_points,
    semialgebraic_maximize,
    semialgebraic_minimize,
)
from .optimization_results import OptimizationResult
from .qe import qe_by_complete_cad
from .reasoning import region_compact
from .solve.find_instance import find_instance
from .structural_keys import symbol_identity_key
from .symbol_resolution import resolve_symbol


def _mapping_tuple(mapping) -> tuple[sp.Expr, ...]:
    if (
        isinstance(mapping, sp.Basic)
        or not isinstance(mapping, Sequence)
        or isinstance(mapping, (str, bytes))
    ):
        return (sp.sympify(mapping),)
    return tuple(map(sp.sympify, mapping))


def semialgebraic_projection(
    region,
    eliminate: Sequence[sp.Symbol | str],
    variables: Sequence[sp.Symbol | str] | None = None,
) -> sp.Expr:
    """Project ``region`` by existentially eliminating the requested variables."""

    formula = normalize_formula(region)
    all_vars = normalize_problem_variables(variables, formula)
    elim = tuple(resolve_symbol(v, context=(formula,), known_symbols=all_vars) for v in eliminate)
    missing = tuple(v for v in elim if v not in all_vars)
    if missing:
        raise ValueError(f"projection variables are not problem variables: {missing!r}")
    kept = tuple(v for v in all_vars if v not in set(elim))
    if not elim:
        return formula
    result = qe_by_complete_cad(
        (*kept, *elim),
        tuple(("exists", v) for v in elim),
        parse_formula(formula),
        return_result=True,
    )
    return sp.simplify(result.formula)


def semialgebraic_image(
    mapping: Sequence[sp.Expr] | sp.Expr,
    domain,
    variables: Sequence[sp.Symbol | str] | None = None,
    *,
    image_variables: Sequence[sp.Symbol | str] | None = None,
    parameters: Sequence[sp.Symbol | str] | None = None,
) -> sp.Expr:
    """Return the exact semialgebraic image of a polynomial/rational map."""

    maps = _mapping_tuple(mapping)
    condition = normalize_formula(domain)
    context = (condition, *maps)
    params = normalize_variables(parameters, *context, append_context_symbols=False)
    source = normalize_variables(
        variables,
        *context,
        append_context_symbols=variables is None,
        exclude=params,
    )
    if image_variables is None:
        targets = tuple(sp.Symbol(f"y{i}", real=True) for i in range(len(maps)))
    else:
        targets = tuple(resolve_symbol(v) for v in image_variables)
    if len(targets) != len(maps):
        raise ValueError("image_variables must have the same length as mapping")
    if set(targets) & set(source):
        raise ValueError("image variables must be distinct from source variables")
    graph = sp.And(condition, *(sp.Eq(y, f) for y, f in zip(targets, maps, strict=True)))
    result = qe_by_complete_cad(
        (*targets, *source),
        tuple(("exists", v) for v in source),
        parse_formula(graph),
        return_result=True,
    )
    return sp.simplify(result.formula)


def semialgebraic_preimage(
    mapping: Sequence[sp.Expr] | sp.Expr,
    target,
    variables: Sequence[sp.Symbol | str] | None = None,
    *,
    target_variables: Sequence[sp.Symbol | str] | None = None,
) -> sp.Expr:
    """Return the preimage of a semialgebraic target under a symbolic map."""

    maps = _mapping_tuple(mapping)
    target_formula = normalize_formula(target)
    if target_variables is None:
        target_vars = tuple(sorted(target_formula.free_symbols, key=symbol_identity_key))
    else:
        target_vars = tuple(
            resolve_symbol(v, context=(target_formula,), known_symbols=target_formula.free_symbols)
            for v in target_variables
        )
    if len(target_vars) != len(maps):
        raise ValueError("target_variables must have the same length as mapping")
    source = normalize_problem_variables(variables, sp.Tuple(*maps))
    extra = target_formula.free_symbols - set(target_vars)
    if extra & set(source):
        raise ValueError(
            "target formula contains ambiguous source symbols; pass target_variables explicitly"
        )
    return sp.simplify(target_formula.xreplace(dict(zip(target_vars, maps, strict=True))))


def fiber(region, substitutions: Mapping[sp.Symbol | str, sp.Expr]) -> sp.Expr:
    """Specialize a semialgebraic family at fixed parameter/coordinate values.

    String keys are resolved against the actual symbols in ``region`` and an
    ambiguous same-name symbol is rejected rather than guessed.
    """

    formula = normalize_formula(region)
    known = tuple(formula.free_symbols)
    resolved: dict[sp.Symbol, sp.Expr] = {}
    for key, value in substitutions.items():
        symbol = resolve_symbol(key, context=(formula,), known_symbols=known)
        if symbol not in formula.free_symbols:
            raise ValueError(f"fiber substitution symbol {symbol!r} is not present in the region")
        resolved[symbol] = sp.sympify(value)
    return sp.simplify(formula.subs(resolved))


@dataclass(frozen=True)
class BoundingBoxResult:
    variables: tuple[sp.Symbol, ...]
    bounds: Mapping[sp.Symbol, tuple[sp.Expr, sp.Expr]]
    attained: Mapping[sp.Symbol, tuple[bool, bool]]
    certified: bool = True


def bounding_box(region, variables=None, *, return_result: bool = False):
    """Compute the exact axis-aligned bounding box by coordinate optimization."""

    formula = normalize_formula(region)
    vars_ = normalize_problem_variables(variables, formula)
    bounds: dict[sp.Symbol, tuple[sp.Expr, sp.Expr]] = {}
    attained: dict[sp.Symbol, tuple[bool, bool]] = {}
    for var in vars_:
        lo = semialgebraic_minimize(var, formula, vars_, return_result=True)
        hi = semialgebraic_maximize(var, formula, vars_, return_result=True)
        if not isinstance(lo, OptimizationResult) or not isinstance(hi, OptimizationResult):
            raise TypeError("coordinate optimization returned a non-OptimizationResult")
        bounds[var] = (lo.value, hi.value)
        attained[var] = (lo.attained, hi.attained)
    result = BoundingBoxResult(vars_, bounds, attained)
    return result if return_result else bounds


@dataclass(frozen=True)
class DistanceToRegionResult:
    distance: sp.Expr
    squared_distance: sp.Expr
    point: Mapping[sp.Symbol, sp.Expr] | None
    nearest_points: tuple[Mapping[sp.Symbol, sp.Expr], ...]
    attained: bool
    optimization: OptimizationResult


def _distance_via_smooth_algebraic_variety(
    point_map: Mapping[sp.Symbol, sp.Expr],
    formula: sp.Expr,
    variables: tuple[sp.Symbol, ...],
) -> OptimizationResult | None:
    """Solve a smooth equality-variety distance problem through its critical system."""

    from .decision import is_satisfiable
    from .exact_arithmetic import compare_exact_reals, exact_truth
    from .region_analysis import _jacobian_rank_singular_formula
    from .solve.zero_dimensional import solve_zero_dimensional_system

    if formula.free_symbols - set(variables):
        return None
    if any(value.free_symbols for value in point_map.values()):
        return None
    try:
        equalities, inequalities, disequalities = _constraint_data(formula, variables)
    except (TypeError, ValueError, NotImplementedError):
        return None
    if not equalities or inequalities or disequalities:
        return None
    # This specialization is deliberately limited to pure algebraic varieties.
    if any(
        not isinstance(atom, sp.Equality)
        for atom in (formula.args if isinstance(formula, sp.And) else (formula,))
    ):
        return None
    singular = _jacobian_rank_singular_formula(equalities, variables)
    if singular not in (False, sp.false):
        try:
            if is_satisfiable(singular, variables, strategy="auto"):
                return None
        except (TypeError, ValueError, ArithmeticError, NotImplementedError):
            return None

    if exact_truth(formula.subs(point_map)):
        return OptimizationResult(
            sp.Integer(0),
            variables,
            sp.Integer(0),
            (dict(point_map),),
            True,
            "min",
            method="distance_critical_points",
            certified=True,
        )

    objective = sp.Add(*((var - point_map[var]) ** 2 for var in variables))
    multipliers = tuple(sp.Dummy(f"lambda{i + 1}", real=True) for i in range(len(equalities)))
    gradient_equations = []
    for var in variables:
        gradient_equations.append(
            sp.expand(
                sp.diff(objective, var)
                - sum(
                    mult * sp.diff(eq, var)
                    for mult, eq in zip(multipliers, equalities, strict=True)
                )
            )
        )
    solve_vars = (*variables, *multipliers)
    try:
        solutions = solve_zero_dimensional_system(
            (*equalities, *gradient_equations), vars=solve_vars, real=True
        )
    except (TypeError, ValueError, ArithmeticError, NotImplementedError):
        return None
    if not solutions:
        return None
    candidates: list[tuple[sp.Expr, dict[sp.Symbol, sp.Expr]]] = []
    for solution in solutions:
        assignment = dict(zip(solve_vars, solution, strict=True))
        point_assignment = {var: assignment[var] for var in variables}
        value = sp.simplify(objective.subs(point_assignment))
        candidates.append((value, point_assignment))
    best_value = candidates[0][0]
    for value, _assignment in candidates[1:]:
        if compare_exact_reals(value, best_value) < 0:
            best_value = value
    best_points = tuple(
        assignment
        for value, assignment in candidates
        if compare_exact_reals(value, best_value) == 0
    )
    return OptimizationResult(
        objective,
        variables,
        sp.simplify(best_value),
        best_points,
        True,
        "min",
        method="distance_critical_points",
        certified=True,
    )


def distance_to_region(point, region, variables=None, *, return_result: bool = False):
    """Compute exact Euclidean distance from a point to a semialgebraic region."""

    formula = normalize_formula(region)
    vars_ = normalize_problem_variables(variables, formula)
    point_map = normalize_point(point, vars_, context=(formula,))
    p = tuple(point_map[v] for v in vars_)
    objective = sp.Add(*((v - c) ** 2 for v, c in zip(vars_, p, strict=True)))
    opt = _distance_via_smooth_algebraic_variety(point_map, formula, vars_)
    if opt is None:
        opt = semialgebraic_minimize(objective, formula, vars_, return_result=True)
    if not isinstance(opt, OptimizationResult):
        raise TypeError("distance optimization returned a non-OptimizationResult")
    distance = sp.sqrt(opt.value)
    result = DistanceToRegionResult(distance, opt.value, point_map, opt.points, opt.attained, opt)
    return result if return_result else distance


@dataclass(frozen=True)
class DistanceBetweenRegionsResult:
    distance: sp.Expr
    squared_distance: sp.Expr
    closest_pairs: tuple[tuple[Mapping[sp.Symbol, sp.Expr], Mapping[sp.Symbol, sp.Expr]], ...]
    attained: bool
    optimization: OptimizationResult


def distance_between_regions(left, right, variables=None, *, return_result: bool = False):
    """Compute exact Euclidean distance between two semialgebraic regions."""

    left_formula = normalize_formula(left)
    right_formula = normalize_formula(right)
    vars_ = normalize_problem_variables(variables, sp.And(left_formula, right_formula))
    right_vars = tuple(fresh_real_dummy(f"{v.name}_right") for v in vars_)
    right_sub = dict(zip(vars_, right_vars, strict=True))
    right_copy = right_formula.xreplace(right_sub)
    objective = sp.Add(*((x - y) ** 2 for x, y in zip(vars_, right_vars, strict=True)))
    opt = semialgebraic_minimize(
        objective, sp.And(left_formula, right_copy), (*vars_, *right_vars), return_result=True
    )
    if not isinstance(opt, OptimizationResult):
        raise TypeError("distance optimization returned a non-OptimizationResult")
    pairs = []
    for assignment in opt.points:
        a = {v: assignment[v] for v in vars_ if v in assignment}
        b = {v: assignment[y] for v, y in zip(vars_, right_vars, strict=True) if y in assignment}
        pairs.append((a, b))
    result = DistanceBetweenRegionsResult(
        sp.sqrt(opt.value), opt.value, tuple(pairs), opt.attained, opt
    )
    return result if return_result else result.distance


def _constant_values_on_positive_dimensional_critical_components(
    objective: sp.Expr,
    condition: sp.Expr,
    variables: tuple[sp.Symbol, ...],
    equalities: tuple[sp.Expr, ...],
    inequalities: tuple[sp.Expr, ...],
) -> tuple[sp.Expr, ...]:
    """Certify objective values on positive-dimensional KKT components.

    Zero-dimensional KKT points are handled by ``_kkt_candidate_points``.
    Here we project each positive-dimensional KKT active-set locus back to the
    original variables, split its feasible part into exact CAD-connected
    components, sample each component exactly, and prove that no point on that
    component has a different objective value.
    """

    from .decision import is_satisfiable
    from .derived_geometry import connected_components

    certified_values: list[sp.Expr] = []
    for locus in iter_positive_dimensional_critical_loci(
        objective, condition, variables, equalities, inequalities, include_singular=False
    ):
        feasible_locus = sp.And(
            condition, *(sp.Eq(equation, 0) for equation in locus.equations), evaluate=False
        )
        try:
            components = connected_components(feasible_locus, variables)
        except (ArithmeticError, TypeError, ValueError, NotImplementedError, sp.PolynomialError):
            continue
        for candidate_component in components:
            try:
                witness_result = find_instance(
                    candidate_component, variables, count=1, strategy="cad", return_result=True
                )
                witness = witness_result.first()
                if witness is None:
                    continue
                candidate_value = sp.simplify(objective.subs(witness))
                different_value_exists = is_satisfiable(
                    sp.And(candidate_component, sp.Ne(objective, candidate_value), evaluate=False),
                    variables,
                    strategy="cad",
                )
            except (
                ArithmeticError,
                TypeError,
                ValueError,
                NotImplementedError,
                sp.PolynomialError,
            ):
                continue
            if different_value_exists is False and candidate_value not in certified_values:
                certified_values.append(candidate_value)
    return tuple(certified_values)


def critical_values(expression, region=sp.true, variables=None) -> tuple[sp.Expr, ...]:
    """Return exact objective values from isolated and constant KKT components."""

    expr = sp.sympify(expression)
    condition = normalize_formula(region)
    vars_ = normalize_problem_variables(variables, sp.Tuple(expr, condition))
    equalities, inequalities, _ = _constraint_data(condition, vars_)
    points = _kkt_candidate_points(expr, condition, vars_, equalities, inequalities)
    values: list[sp.Expr] = []
    for point in points:
        value = sp.simplify(expr.subs(point))
        if value not in values:
            values.append(value)
    for value in _constant_values_on_positive_dimensional_critical_components(
        expr, condition, vars_, equalities, inequalities
    ):
        if value not in values:
            values.append(value)
    try:
        lo = semialgebraic_minimize(expr, condition, vars_, return_result=True)
        hi = semialgebraic_maximize(expr, condition, vars_, return_result=True)
        for result in (lo, hi):
            if (
                isinstance(result, OptimizationResult)
                and result.attained
                and result.value not in values
            ):
                values.append(result.value)
    except (ArithmeticError, TypeError, ValueError, NotImplementedError, sp.PolynomialError):
        pass
    return tuple(sorted(values, key=sp.default_sort_key))


def is_path_connected(region, variables=None, *, max_pair_checks: int | None = None) -> bool:
    """Decide path connectedness via exact CAD connectivity.

    Semialgebraic connected sets are semialgebraically path connected, so the
    existing certified CAD connected-component graph suffices for the decision.
    """

    formula = normalize_formula(region)
    vars_ = normalize_problem_variables(variables, formula)
    if max_pair_checks is None:
        factorized = factorized_equality_components(formula, vars_)
        if factorized is not None:
            return len(factorized) <= 1
    graph = extract_cad_connectivity(formula, vars_, max_pair_checks=max_pair_checks)
    if graph.component_count > 1 and max_pair_checks is not None:
        raise ResourceLimitError(
            "limited CAD adjacency checks cannot certify disconnection; rerun without max_pair_checks"
        )
    return graph.component_count <= 1


@dataclass(frozen=True)
class CADPathResult:
    connected: bool
    explicit: bool
    waypoints: tuple[Mapping[sp.Symbol, sp.Expr], ...]
    cell_indices: tuple[tuple[int, ...], ...]
    graph: CADConnectivityGraph
    note: str = ""


def _point_cell_index(
    graph: CADConnectivityGraph, point: Mapping[sp.Symbol, sp.Expr]
) -> int | None:
    for i, cell in enumerate(graph.cells):
        formula = cell.as_formula(closed=False)
        val = sp.simplify(formula.subs(point))
        if is_true_expr(val):
            return i
    return None


def path_between(
    region, start, end, variables=None, *, max_pair_checks: int | None = None
) -> CADPathResult:
    """Return a certified CAD cell-chain connecting two points in a region.

    For a one-dimensional connected component the waypoint chain is an explicit
    piecewise-linear path. In higher dimensions the current implementation
    returns the certified cell/connector chain supplied by CAD connectivity;
    this is intentionally not mislabeled as a full Canny-style roadmap path.
    """

    formula = normalize_formula(region)
    vars_ = normalize_problem_variables(variables, formula)
    start_map = normalize_point(start, vars_, context=(formula,))
    end_map = normalize_point(end, vars_, context=(formula,))
    graph = extract_cad_connectivity(formula, vars_, max_pair_checks=max_pair_checks)
    s = _point_cell_index(graph, start_map)
    e = _point_cell_index(graph, end_map)
    if s is None or e is None:
        raise ValueError("both endpoints must lie in the region")
    if s == e:
        return CADPathResult(
            True,
            len(vars_) == 1 or is_convex(formula, vars_),
            (start_map, end_map),
            (graph.cells[s].index,),
            graph,
        )
    adjacency: dict[int, list[tuple[int, object]]] = {i: [] for i in range(len(graph.cells))}
    for edge in graph.edges:
        adjacency[edge.left].append((edge.right, edge))
        adjacency[edge.right].append((edge.left, edge))
    queue = deque([s])
    previous: dict[int, tuple[int, object] | None] = {s: None}
    while queue and e not in previous:
        current = queue.popleft()
        for nxt, edge in adjacency[current]:
            if nxt not in previous:
                previous[nxt] = (current, edge)
                queue.append(nxt)
    if e not in previous:
        if max_pair_checks is not None:
            raise ResourceLimitError(
                "limited CAD adjacency checks cannot certify that the endpoints are disconnected; rerun without max_pair_checks"
            )
        return CADPathResult(
            False, False, (), (), graph, "endpoints lie in different CAD connected components"
        )
    chain = [e]
    edges = []
    while chain[-1] != s:
        prev, edge = previous[chain[-1]]  # type: ignore[misc]
        edges.append(edge)
        chain.append(prev)
    chain.reverse()
    edges.reverse()
    waypoints: list[Mapping[sp.Symbol, sp.Expr]] = [start_map]
    for edge in edges:
        witness = find_instance(
            edge.connector_formula, vars_, count=1, strategy="cad", return_result=True
        )
        if witness.found and witness.first() is not None:
            waypoints.append(witness.first())
    waypoints.append(end_map)
    explicit = len(vars_) == 1
    note = (
        "explicit piecewise-linear path"
        if explicit
        else "certified CAD cell/connector chain; full explicit algebraic roadmap parameterization is outside this API"
    )
    return CADPathResult(
        True, explicit, tuple(waypoints), tuple(graph.cells[i].index for i in chain), graph, note
    )


def euler_characteristic(region, variables=None, *, compact_support: bool = True) -> int:
    """Compute the semialgebraic Euler characteristic from selected CAD cells.

    The additive cell formula ``sum (-1)^dim`` is the compactly-supported Euler
    characteristic. For compact regions it equals the ordinary Euler
    characteristic. ``compact_support=False`` is therefore accepted only when
    the region is certified bounded (closedness is left to the caller's formula
    semantics; compact semialgebraic inputs are the intended use).
    """

    formula = normalize_formula(region)
    vars_ = normalize_problem_variables(variables, formula)
    if not compact_support and not region_compact(formula, vars_):
        raise ValueError(
            "ordinary Euler characteristic is only exposed here for compact regions; use compact_support=True"
        )
    solution = extract_cylindrical_solution(formula, vars_, selected_only=True)
    return int(sum((-1) ** cell.dimension for cell in solution.cells if cell.selected))


__all__ = [
    "BoundingBoxResult",
    "DistanceToRegionResult",
    "DistanceBetweenRegionsResult",
    "CADPathResult",
    "semialgebraic_projection",
    "semialgebraic_image",
    "semialgebraic_preimage",
    "fiber",
    "critical_values",
    "distance_to_region",
    "distance_between_regions",
    "bounding_box",
    "is_convex",
    "is_path_connected",
    "path_between",
    "euler_characteristic",
]
