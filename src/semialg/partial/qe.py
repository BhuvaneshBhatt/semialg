from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field
from time import perf_counter

import sympy as sp

from .._immutable import FrozenMappingPickle, freeze_mapping
from .._zero_testing import certified_pointwise_zero
from ..algebraic.comparison import sort_samples
from ..algebraic.roots import isolate_real_roots
from ..algebraic.samples import AlgebraicRoot, FiberRootContext, sample_to_expr
from ..cad_algorithms.decomposition import _build_stack, _stack_roots_over_point
from ..cad_algorithms.lifting.stack import CADCell
from ..cad_algorithms.polynomial_utils import exact_univariate_poly
from ..cad_algorithms.projection.collins import ProjectionTower, build_collins_proj_set
from ..cad_algorithms.resource_limits import CADBudget, budget_scope, current_budget
from ..errors import InternalInvariantError, ProjectionNullification
from ..exact_arithmetic import exact_truth
from ..formula import Formula, equational_constraints, formula_polynomials, parse_formula, to_sympy
from ..qe.complete import cells_to_formula, evaluate_formula_on_cell, norm_internal_order
from .pruning import evaluate_pruning_status


@dataclass(frozen=True)
class LazyCADStats:
    """Small execution record for lazy CAD/QE traversals."""

    variables: tuple[sp.Symbol, ...]
    quantifiers: tuple[tuple[str, sp.Symbol], ...]
    visited_cells_by_level: Mapping[int, int]
    lifted_stacks: int
    evaluated_leaf_cells: int
    stopped_early: bool
    elapsed_seconds: float
    pruned_prefix_cells: int = 0
    ec_pruned_cells: int = 0
    ec_section_lifts: int = 0
    ec_skipped_sector_cells: int = 0
    derived_ec_count: int = 0
    notes: tuple[str, ...] = ()


@dataclass(frozen=True)
class LazyResolveResult:
    """Truth result produced by lazy quantifier-guided CAD lifting."""

    truth_value: bool
    status: str
    backend: str
    stats: LazyCADStats
    witness: Mapping[sp.Symbol, sp.Expr] | None = None
    counterexample: Mapping[sp.Symbol, sp.Expr] | None = None


@dataclass(frozen=True)
class LazyInstanceResult:
    """Satisfying sample returned by lazy existential CAD search."""

    instance: Mapping[sp.Symbol, sp.Expr] | None
    status: str
    backend: str
    stats: LazyCADStats
    found: bool = False


@dataclass
class _TraversalState:
    tower: ProjectionTower
    variables: tuple[sp.Symbol, ...]
    matrix: Formula
    cells_by_level: dict[int, list[CADCell]] = field(default_factory=dict)
    lifted_stacks: int = 0
    evaluated_leaf_cells: int = 0
    stopped_early: bool = False
    pruned_prefix_cells: int = 0
    ec_pruned_cells: int = 0
    ec_section_lifts: int = 0
    ec_skipped_sector_cells: int = 0
    equational_constraints: tuple[sp.Expr, ...] = ()
    ec_by_level: Mapping[int, tuple[sp.Expr, ...]] = field(default_factory=dict)
    derived_ec_count: int = 0

    def add_cells(self, level: int, cells: Sequence[CADCell]) -> None:
        self.cells_by_level.setdefault(level, []).extend(cells)
        budget = current_budget()
        if budget is not None:
            budget.check(sum(map(len, self.cells_by_level.values())) + self.ec_skipped_sector_cells)

    def stats(
        self,
        quantifiers: Sequence[tuple[str, sp.Symbol]],
        start: float,
        *,
        notes: Sequence[str] = (),
    ) -> LazyCADStats:
        return LazyCADStats(
            variables=self.variables,
            quantifiers=tuple(quantifiers),
            visited_cells_by_level={
                level: len(cells) for level, cells in sorted(self.cells_by_level.items())
            },
            lifted_stacks=self.lifted_stacks,
            evaluated_leaf_cells=self.evaluated_leaf_cells,
            stopped_early=self.stopped_early,
            elapsed_seconds=perf_counter() - start,
            pruned_prefix_cells=self.pruned_prefix_cells,
            ec_pruned_cells=self.ec_pruned_cells,
            ec_section_lifts=self.ec_section_lifts,
            ec_skipped_sector_cells=self.ec_skipped_sector_cells,
            derived_ec_count=self.derived_ec_count,
            notes=tuple(notes),
        )


def _main_level(expr: sp.Expr, variables: Sequence[sp.Symbol]) -> int:
    positions = [i + 1 for i, var in enumerate(variables) if var in sp.sympify(expr).free_symbols]
    return max(positions, default=0)


def _normalize_ec(expr: sp.Expr) -> sp.Expr | None:
    expr = sp.expand(sp.sympify(expr))
    if expr == 0:
        return None
    try:
        _, primitive = sp.primitive(expr)
        expr = sp.expand(primitive)
    except (sp.PolynomialError, ValueError, TypeError):
        pass
    if expr.could_extract_minus_sign():
        expr = -expr
    return expr


def _propagate_equational_constraints(
    constraints: Sequence[sp.Expr], variables: Sequence[sp.Symbol]
) -> tuple[dict[int, tuple[sp.Expr, ...]], int]:
    """Derive lower-level ECs that are logically necessary for all solutions.

    Pairwise resultants of conjunctively required equations are necessary for a
    common root in the eliminated variable.  We recursively propagate only such
    necessary equalities; discriminants/coefficients are not used
    as ECs because their vanishing is not necessary for a single equation.
    """
    vars_tuple = tuple(variables)
    known: dict[sp.Expr, sp.Expr] = {}
    original_keys: set[sp.Expr] = set()
    for expr in constraints:
        norm = _normalize_ec(expr)
        if norm is not None:
            key = norm
            known[key] = norm
            original_keys.add(key)
    changed = True
    while changed:
        changed = False
        current = tuple(known.values())
        for level in range(len(vars_tuple), 1, -1):
            var = vars_tuple[level - 1]
            active = [
                expr
                for expr in current
                if _main_level(expr, vars_tuple) == level and var in expr.free_symbols
            ]
            if len(active) < 2:
                continue
            for i, left in enumerate(active):
                for right in active[i + 1 :]:
                    try:
                        resultant = sp.resultant(left, right, var)
                    except (sp.PolynomialError, ValueError, TypeError, NotImplementedError):
                        continue
                    norm = _normalize_ec(resultant)
                    if norm is None or norm.is_number:
                        continue
                    key = norm
                    if key not in known and _main_level(norm, vars_tuple) < level:
                        known[key] = norm
                        changed = True
        # New resultants can themselves combine at lower levels.
    by_level: dict[int, list[sp.Expr]] = {}
    for expr in known.values():
        level = _main_level(expr, vars_tuple)
        if level > 0:
            by_level.setdefault(level, []).append(expr)
    ordered = {
        level: tuple(sorted(items, key=lambda e: (sp.count_ops(e), sp.default_sort_key(e))))
        for level, items in by_level.items()
    }
    derived = sum(1 for key in known if key not in original_keys)
    return ordered, derived


def _quantifier_map(quantifiers: Sequence[tuple[str, sp.Symbol]]) -> dict[sp.Symbol, str]:
    out: dict[sp.Symbol, str] = {}
    for q, sym in quantifiers:
        q = q.lower()
        if q not in {"exists", "forall"}:
            raise ValueError(f"unsupported quantifier {q!r}")
        out[sym] = q
    return out


def _specialized_ec_roots(
    state: _TraversalState, parent: CADCell | None, level: int
) -> tuple[object, ...] | None:
    ecs = state.ec_by_level.get(level, ())
    if not ecs:
        return None
    variables = state.variables
    var = variables[level - 1]
    substitutions = {} if parent is None else _prefix_mapping(variables, parent)
    # One necessary EC is sufficient to restrict an existential lift; choose
    # the structurally cheapest.  Other ECs are checked by prefix pruning.
    for ec in ecs:
        specialized = sp.expand(ec.subs(substitutions))
        if specialized == 0:
            continue
        if specialized.free_symbols - {var}:
            continue
        try:
            poly = exact_univariate_poly(specialized, var, algebraic_extension=False)
        except (sp.PolynomialError, ValueError, TypeError):
            continue
        if poly.degree() <= 0:
            return tuple()
        roots = tuple(isolate_real_roots(poly))
        if parent is None:
            return roots
        context = FiberRootContext(
            polynomial=sp.Poly(ec, *variables[:level]),
            variables=tuple(variables[:level]),
            parent_samples=parent.sample,
            fiber_variable=var,
        )
        return tuple(
            AlgebraicRoot(
                root.polynomial,
                root.interval,
                root.root_index,
                root.multiplicity_hint,
                root.root_expr,
                context,
            )
            for root in roots
        )
    return None


def _stack_for_parent(
    state: _TraversalState, parent: CADCell | None, level: int, *, ec_sections_only: bool = False
) -> tuple[CADCell, ...]:
    tower = state.tower
    variables = state.variables
    ec_roots = _specialized_ec_roots(state, parent, level) if ec_sections_only else None
    if ec_roots is not None:
        if (
            parent is not None
            and level < len(variables)
            and any(position % 2 == 0 for position in parent.index)
        ):
            # Restricting an existential lift to EC sections does not discharge
            # order-invariance requirements for other projection polynomials.
            for poly in tower.level(level).polynomials:
                if poly.degree(variables[level - 1]) <= 0:
                    continue
                context = FiberRootContext(
                    polynomial=poly,
                    variables=variables[:level],
                    parent_samples=parent.sample,
                    fiber_variable=variables[level - 1],
                )
                if context.is_nullified_at_parent():
                    raise ProjectionNullification(
                        "EC-restricted CAD needs refinement of a positive-dimensional "
                        "projection-nullification cell"
                    )
        roots = sort_samples(tuple(ec_roots))
        full = _build_stack(parent, roots, level, tower)
        stack = tuple(cell for cell in full if cell.kind == "section")
        state.ec_section_lifts += 1
        skipped = len(full) - len(stack)
        state.ec_skipped_sector_cells += skipped
        state.ec_pruned_cells += skipped
    elif level == 1:
        roots = []
        for poly in tower.level(1).polynomials:
            roots.extend(isolate_real_roots(poly))
        stack = _build_stack(None, sort_samples(tuple(roots)), 1, tower)
    else:
        if parent is None:
            raise InternalInvariantError("higher-level CAD lifting requires a parent cell")
        roots = _stack_roots_over_point(
            tower.level(level).polynomials,
            variables,
            parent.sample,
            variables[level - 1],
            use_lazard=tower.metadata.get("projection") in {"lazard", "lazard-complete"},
            delineate_nullification=level < len(variables),
            base_is_point=all(position % 2 == 1 for position in parent.index),
        )
        stack = _build_stack(parent, roots, level, tower)
    state.lifted_stacks += 1
    state.add_cells(level, stack)
    return stack


def _sample_mapping(variables: Sequence[sp.Symbol], cell: CADCell) -> dict[sp.Symbol, sp.Expr]:
    return {
        sym: sample_to_expr(sample) for sym, sample in zip(variables, cell.sample, strict=False)
    }


def _prefix_mapping(variables: Sequence[sp.Symbol], cell: CADCell) -> dict[sp.Symbol, sp.Expr]:
    return {
        sym: sample_to_expr(sample)
        for sym, sample in zip(variables[: cell.level], cell.sample, strict=True)
    }


def _ec_false_at_prefix(state: _TraversalState, cell: CADCell) -> bool:
    if not state.ec_by_level:
        return False
    subs = _prefix_mapping(state.variables, cell)
    for level, ecs in state.ec_by_level.items():
        if level > cell.level:
            continue
        for ec in ecs:
            value = sp.simplify(sp.expand(ec).subs(subs))
            if certified_pointwise_zero(value) is False:
                return True
    return False


def _prefix_truth(state: _TraversalState, cell: CADCell) -> bool | None:
    """Truth determined on a prefix CAD cell, if Boolean short-circuiting allows it.

    Projection sign-invariance makes sample truth valid for atoms whose variables
    are all present in the prefix. Atoms involving later variables remain unknown.
    """
    if _ec_false_at_prefix(state, cell):
        state.pruned_prefix_cells += 1
        state.ec_pruned_cells += 1
        return False
    decision = evaluate_pruning_status(state.matrix, _prefix_mapping(state.variables, cell))
    if decision.should_prune:
        state.pruned_prefix_cells += 1
        if decision.status is False and _ec_false_at_prefix(state, cell):
            state.ec_pruned_cells += 1
        return decision.status
    return None


def lazy_resolve_formula(
    vars_: Sequence[sp.Symbol],
    quantifiers: Sequence[tuple[str, sp.Symbol]],
    matrix: Formula,
) -> LazyResolveResult:
    """Resolve a closed real-polynomial formula by quantifier-guided lazy CAD.

    Projection is still built up front, but stacks are lifted only along paths
    needed by the quantifier search. Existential branches stop on the first
    true child; universal branches stop on the first false child.
    """

    start = perf_counter()
    quantifiers = tuple(quantifiers)
    variables, free, _quantified, notes = norm_internal_order(vars_, quantifiers, matrix, None)
    if free:
        raise ValueError(
            "lazy_resolve_formula requires a sentence; use reduce for formulas with free variables"
        )
    if not variables:
        truth = exact_truth(to_sympy(matrix))
        state = _TraversalState(ProjectionTower((), (), ()), (), matrix)
        return LazyResolveResult(
            truth,
            "complete",
            "partial-collins-lazy",
            state.stats(quantifiers, start, notes=notes),
            witness={} if truth else None,
            counterexample=None if truth else {},
        )
    quantifiers = tuple((q.lower(), sym) for q, sym in quantifiers)
    qmap = _quantifier_map(quantifiers)
    if variables and len(set(qmap.values())) == 1:
        # Reordering inside a single quantifier block preserves semantics and can
        # dramatically reduce Collins projection size. The planner measures exact
        # projection proxies and is deterministic.
        try:
            from ..planner.features import extract_problem_features
            from ..planner.heuristics import choose_formula_variable_order

            features = extract_problem_features(
                matrix, variables=variables, quantifiers=quantifiers
            )
            variables = tuple(choose_formula_variable_order(matrix, features))
            quantifiers = tuple((qmap[var], var) for var in variables)
            qmap = _quantifier_map(quantifiers)
            notes = (*notes, "single-block CAD variable order optimized")
        except (sp.PolynomialError, ValueError, TypeError):
            pass
    polys = formula_polynomials(matrix) or [sp.Integer(1)]
    budget = current_budget()
    if budget is not None:
        budget.check()
    tower = build_collins_proj_set(polys, variables)
    if budget is not None:
        budget.check()
    ecs = tuple(equational_constraints(matrix))
    ecs_by_level, derived_ec_count = _propagate_equational_constraints(ecs, variables)
    state = _TraversalState(
        tower=tower,
        variables=variables,
        matrix=matrix,
        equational_constraints=ecs,
        ec_by_level=ecs_by_level,
        derived_ec_count=derived_ec_count,
    )
    witness: dict[sp.Symbol, sp.Expr] | None = None
    counterexample: dict[sp.Symbol, sp.Expr] | None = None

    def visit_truth_path(level: int, parent: CADCell | None) -> tuple[bool, CADCell | None]:
        """Evaluate quantified truth while lifting only logically needed cells."""

        nonlocal witness, counterexample
        if parent is not None:
            prefix_truth = _prefix_truth(state, parent)
            if prefix_truth is not None:
                return prefix_truth, parent
        if level > len(variables):
            if parent is None:
                raise InternalInvariantError(
                    "quantified truth evaluation reached a leaf without a cell"
                )
            state.evaluated_leaf_cells += 1
            truth = evaluate_formula_on_cell(matrix, parent, variables)
            if truth and witness is None:
                witness = {k: v for k, v in _sample_mapping(variables, parent).items()}
            if not truth and counterexample is None:
                counterexample = {k: v for k, v in _sample_mapping(variables, parent).items()}
            return truth, parent

        q = qmap[variables[level - 1]]
        stack = _stack_for_parent(state, parent, level, ec_sections_only=(q == "exists"))
        if q == "exists":
            saw_any = False
            last_cell = None
            for child in stack:
                value, leaf = visit_truth_path(level + 1, child)
                saw_any = saw_any or value
                last_cell = leaf or child
                if value:
                    state.stopped_early = True
                    return True, leaf or child
            return saw_any, last_cell
        if q == "forall":
            last_cell = None
            for child in stack:
                value, leaf = visit_truth_path(level + 1, child)
                last_cell = leaf or child
                if not value:
                    state.stopped_early = True
                    return False, leaf or child
            return True, last_cell
        raise ValueError(f"unsupported quantifier: {q!r}")

    try:
        truth, leaf = visit_truth_path(1, None)
    except ProjectionNullification:
        if budget is not None:
            budget.restart()
        from ..cad_algorithms.projection.lazard_complete import build_lazard_complete_projection

        state = _TraversalState(
            build_lazard_complete_projection(polys, variables), variables, matrix
        )
        witness = counterexample = None
        notes = (*notes, "restarted with full Lazard projection and native valuation lifting")
        truth, leaf = visit_truth_path(1, None)
    if leaf is not None:
        mapping = {k: v for k, v in _sample_mapping(variables, leaf).items()}
        if truth:
            witness = witness or mapping
        else:
            counterexample = counterexample or mapping
    if budget is not None:
        budget.complete(
            sum(map(len, state.cells_by_level.values())) + state.ec_skipped_sector_cells
        )
    return LazyResolveResult(
        truth_value=truth,
        status="complete",
        backend="partial-lazard-lazy"
        if state.tower.metadata.get("projection") == "lazard-complete"
        else "partial-collins-lazy",
        stats=state.stats(quantifiers, start, notes=notes),
        witness=witness,
        counterexample=counterexample,
    )


def lazy_find_inst_form(
    vars_: Sequence[sp.Symbol],
    matrix: Formula,
    *,
    quantifiers: Sequence[tuple[str, sp.Symbol]] | None = None,
) -> LazyInstanceResult:
    """Find one satisfying real sample by lazy existential CAD search."""

    start = perf_counter()
    quantifiers = tuple(quantifiers) if quantifiers is not None else ()
    if not quantifiers:
        variables = tuple(vars_)
        quantifiers = tuple(("exists", sym) for sym in variables)
    else:
        variables, free, _quantified, _notes = norm_internal_order(vars_, quantifiers, matrix, None)
        if free:
            # For instance finding, treat free variables as existential search
            # variables ahead of explicitly quantified ones.
            quantifiers = tuple(("exists", sym) for sym in free) + tuple(
                (q.lower(), sym) for q, sym in quantifiers
            )
        else:
            quantifiers = tuple((q.lower(), sym) for q, sym in quantifiers)
    qmap = _quantifier_map(quantifiers)
    if not variables:
        truth = exact_truth(to_sympy(matrix))
        state = _TraversalState(ProjectionTower((), (), ()), (), matrix)
        return LazyInstanceResult(
            {} if truth else None,
            "complete",
            "partial-collins-lazy",
            state.stats(quantifiers, start),
            found=truth,
        )
    if variables and len(set(qmap.values())) == 1:
        try:
            from ..planner.features import extract_problem_features
            from ..planner.heuristics import choose_formula_variable_order

            features = extract_problem_features(
                matrix, variables=variables, quantifiers=quantifiers
            )
            variables = tuple(choose_formula_variable_order(matrix, features))
            quantifiers = tuple((qmap[var], var) for var in variables)
            qmap = _quantifier_map(quantifiers)
        except (sp.PolynomialError, ValueError, TypeError):
            pass
    polys = formula_polynomials(matrix) or [sp.Integer(1)]
    budget = current_budget()
    if budget is not None:
        budget.check()
    tower = build_collins_proj_set(polys, variables)
    if budget is not None:
        budget.check()
    ecs = tuple(equational_constraints(matrix))
    ecs_by_level, derived_ec_count = _propagate_equational_constraints(ecs, variables)
    state = _TraversalState(
        tower=tower,
        variables=variables,
        matrix=matrix,
        equational_constraints=ecs,
        ec_by_level=ecs_by_level,
        derived_ec_count=derived_ec_count,
    )
    found_cell: CADCell | None = None

    def visit_instance_path(level: int, parent: CADCell | None) -> bool:
        nonlocal found_cell
        if parent is not None:
            prefix_truth = _prefix_truth(state, parent)
            if prefix_truth is not None:
                if prefix_truth:
                    found_cell = parent
                return prefix_truth
        if level > len(variables):
            if parent is None:
                raise InternalInvariantError("instance search reached a leaf without a cell")
            state.evaluated_leaf_cells += 1
            truth = evaluate_formula_on_cell(matrix, parent, variables)
            if truth:
                found_cell = parent
            return truth
        q = qmap.get(variables[level - 1], "exists")
        stack = _stack_for_parent(state, parent, level, ec_sections_only=(q == "exists"))
        if q == "exists":
            for child in stack:
                if visit_instance_path(level + 1, child):
                    state.stopped_early = True
                    return True
            return False
        if q == "forall":
            for child in stack:
                if not visit_instance_path(level + 1, child):
                    return False
            return True
        raise ValueError(f"unsupported quantifier: {q!r}")

    try:
        found = visit_instance_path(1, None) if variables else False
    except ProjectionNullification:
        if budget is not None:
            budget.restart()
        from ..cad_algorithms.projection.lazard_complete import build_lazard_complete_projection

        state = _TraversalState(
            build_lazard_complete_projection(polys, variables), variables, matrix
        )
        found_cell = None
        found = visit_instance_path(1, None)
    instance = _sample_mapping(variables, found_cell) if found_cell is not None and found else None
    if instance is not None:
        for variable in variables:
            instance.setdefault(variable, sp.S.Zero)
    if budget is not None:
        budget.complete(
            sum(map(len, state.cells_by_level.values())) + state.ec_skipped_sector_cells
        )
    return LazyInstanceResult(
        instance=instance,
        status="complete",
        backend="partial-lazard-lazy"
        if state.tower.metadata.get("projection") == "lazard-complete"
        else "partial-collins-lazy",
        stats=state.stats(quantifiers, start, notes=("lazy existential instance search",)),
        found=bool(found and instance is not None),
    )


__all__ = [
    "LazyCADStats",
    "LazyInstanceResult",
    "LazyResolveResult",
    "lazy_find_inst_form",
    "lazy_resolve_formula",
]


@dataclass(frozen=True)
class LazyQEResult(FrozenMappingPickle):
    """Free-cell truth table with a partially lifted quantified suffix.

    Every free-coordinate cell is present. Quantified stacks may be absent
    after a decisive branch, so ``cells_by_level`` is not a complete ambient
    CAD. The projection tower and original matrix retain the exact statement.
    """

    formula: sp.Expr
    variables: tuple[sp.Symbol, ...]
    free_variables: tuple[sp.Symbol, ...]
    quantifiers: tuple[tuple[str, sp.Symbol], ...]
    matrix: sp.Expr
    tower: ProjectionTower
    cells_by_level: Mapping[int, tuple[CADCell, ...]]
    truth_by_index: Mapping[tuple[int, ...], bool]
    stats: LazyCADStats
    status: str = "complete"
    backend: str = "partial-collins-free"

    def __post_init__(self):
        object.__setattr__(self, "cells_by_level", freeze_mapping(self.cells_by_level))
        object.__setattr__(self, "truth_by_index", freeze_mapping(self.truth_by_index))


def _quantified_truth(state, level, parent, quantifiers):
    """Resolve one fiber without dropping any cells of its free-coordinate base."""
    if parent is not None:
        truth = _prefix_truth(state, parent)
        if truth is not None:
            return truth
    if level > len(state.variables):
        state.evaluated_leaf_cells += 1
        return evaluate_formula_on_cell(state.matrix, parent, state.variables)
    existential = quantifiers[state.variables[level - 1]] == "exists"
    stack = _stack_for_parent(state, parent, level, ec_sections_only=existential)
    for child in stack:
        truth = _quantified_truth(state, level + 1, child, quantifiers)
        if truth == existential:
            state.stopped_early = True
            return truth
    return not existential


def lazy_quantifier_eliminate(
    variables, quantifiers, matrix, *, max_cells=None, timeout=None, return_result=False
):
    """Eliminate an alternating real prefix with complete free-cell truth coverage.

    Projection is complete; only the quantified suffix is lifted partially.
    Ordinary projection nullification restarts the entire computation with
    full Lazard projection and valuation lifting. Budgets include both attempts.
    Set ``return_result=True`` to retain the projection and free-cell truth table.
    """
    quantifiers = tuple((kind.lower(), variable) for kind, variable in quantifiers)
    matrix = matrix if isinstance(matrix, Formula) else parse_formula(matrix)
    variables, free, _quantified, notes = norm_internal_order(variables, quantifiers, matrix)
    if not free:
        raise ValueError("use lazy_resolve_formula for a closed sentence")
    budget = (
        CADBudget(max_cells, timeout)
        if max_cells is not None or timeout is not None
        else current_budget()
    )
    start = perf_counter()
    polynomials = formula_polynomials(matrix) or [sp.Integer(1)]
    qmap = _quantifier_map(quantifiers)
    ecs = tuple(equational_constraints(matrix))
    ecs_by_level, derived_count = _propagate_equational_constraints(ecs, variables)

    def construct(builder):
        if budget is not None:
            budget.check()
        tower = builder(polynomials, variables)
        if budget is not None:
            budget.check()
        ordinary = tower.metadata.get("projection") != "lazard-complete"
        state = _TraversalState(
            tower,
            variables,
            matrix,
            equational_constraints=ecs if ordinary else (),
            ec_by_level=ecs_by_level if ordinary else {},
            derived_ec_count=derived_count if ordinary else 0,
        )
        parents = (None,)
        for level in range(1, len(free) + 1):
            parents = tuple(
                child for parent in parents for child in _stack_for_parent(state, parent, level)
            )
        truth = {
            parent.index: _quantified_truth(state, len(free) + 1, parent, qmap)
            for parent in parents
        }
        selected = tuple(parent for parent in parents if truth[parent.index])
        formula = cells_to_formula(selected, free, state.cells_by_level)
        if budget is not None:
            budget.complete(
                sum(map(len, state.cells_by_level.values())) + state.ec_skipped_sector_cells
            )
        return state, truth, formula

    with budget_scope(budget):
        try:
            state, truth, formula = construct(build_collins_proj_set)
        except ProjectionNullification:
            from ..cad_algorithms.projection.lazard_complete import build_lazard_complete_projection

            if budget is not None:
                budget.restart()
            notes = (*notes, "restarted with full Lazard projection and native valuation lifting")
            state, truth, formula = construct(build_lazard_complete_projection)
    result = LazyQEResult(
        formula,
        variables,
        free,
        quantifiers,
        to_sympy(matrix),
        state.tower,
        {level: tuple(cells) for level, cells in state.cells_by_level.items()},
        truth,
        state.stats(quantifiers, start, notes=notes),
        backend="partial-lazard-free"
        if state.tower.metadata.get("projection") == "lazard-complete"
        else "partial-collins-free",
    )
    return result if return_result else result.formula
