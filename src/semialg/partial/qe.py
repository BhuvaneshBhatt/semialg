from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field
from time import perf_counter

import sympy as sp

from ..algebraic.comparison import sort_samples
from ..algebraic.roots import isolate_real_roots
from ..algebraic.samples import sample_to_expr
from ..cad.decomposition import _build_stack, _stack_roots_over_point
from ..cad.lifting.stack import CADCell
from ..cad.projection.collins import ProjectionTower, build_collins_proj_set
from ..formula import Formula, equational_constraints, formula_polynomials
from ..qe.complete import evaluate_formula_on_cell, norm_internal_order
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
    witness: Mapping[str, sp.Expr] | None = None
    counterexample: Mapping[str, sp.Expr] | None = None


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
    necessary equalities; discriminants/coefficients are deliberately not used
    as ECs because their vanishing is not necessary for a single equation.
    """
    vars_tuple = tuple(variables)
    known: dict[str, sp.Expr] = {}
    original_keys: set[str] = set()
    for expr in constraints:
        norm = _normalize_ec(expr)
        if norm is not None:
            key = sp.srepr(norm)
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
                    key = sp.srepr(norm)
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
        level: tuple(sorted(items, key=lambda e: (sp.count_ops(e), len(sp.sstr(e)), sp.sstr(e))))
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
            poly = sp.Poly(specialized, var, domain="EX")
        except (sp.PolynomialError, ValueError, TypeError):
            continue
        if poly.degree() <= 0:
            return tuple()
        return tuple(isolate_real_roots(poly))
    return None


def _stack_for_parent(
    state: _TraversalState, parent: CADCell | None, level: int, *, ec_sections_only: bool = False
) -> tuple[CADCell, ...]:
    tower = state.tower
    variables = state.variables
    ec_roots = _specialized_ec_roots(state, parent, level) if ec_sections_only else None
    if ec_roots is not None:
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
        assert parent is not None
        roots = _stack_roots_over_point(
            tower.level(level).polynomials,
            variables,
            parent.sample,
            variables[level - 1],
            use_lazard=str(tower.metadata.get("projection", "")) == "lazard",
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
            if value != 0 and value.is_zero is False:
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
    variables, free, quantified, notes = norm_internal_order(vars_, quantifiers, matrix, None)
    if free:
        raise ValueError(
            "lazy_resolve_formula requires a sentence; use reduce for formulas with free variables"
        )
    quantifiers = tuple((q.lower(), sym) for q, sym in quantifiers)
    qmap = _quantifier_map(quantifiers)
    polys = formula_polynomials(matrix) or [sp.Integer(1)]
    tower = build_collins_proj_set(polys, variables)
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
    witness: dict[str, sp.Expr] | None = None
    counterexample: dict[str, sp.Expr] | None = None

    def visit_truth_path(level: int, parent: CADCell | None) -> tuple[bool, CADCell | None]:
        """Evaluate quantified truth while lifting only logically needed cells."""

        nonlocal witness, counterexample
        if parent is not None:
            prefix_truth = _prefix_truth(state, parent)
            if prefix_truth is not None:
                return prefix_truth, parent
        if level > len(variables):
            assert parent is not None
            state.evaluated_leaf_cells += 1
            truth = evaluate_formula_on_cell(matrix, parent, variables)
            if truth and witness is None:
                witness = {sp.sstr(k): v for k, v in _sample_mapping(variables, parent).items()}
            if not truth and counterexample is None:
                counterexample = {
                    sp.sstr(k): v for k, v in _sample_mapping(variables, parent).items()
                }
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

    if not variables:
        truth = bool(sp.simplify(sp.sympify(matrix)))
        return LazyResolveResult(
            truth, "complete", "partial-collins-lazy", state.stats(quantifiers, start, notes=notes)
        )

    truth, leaf = visit_truth_path(1, None)
    if leaf is not None:
        mapping = {sp.sstr(k): v for k, v in _sample_mapping(variables, leaf).items()}
        if truth:
            witness = witness or mapping
        else:
            counterexample = counterexample or mapping
    return LazyResolveResult(
        truth_value=truth,
        status="complete",
        backend="partial-collins-lazy",
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
    if quantifiers is None or len(tuple(quantifiers)) == 0:
        variables = tuple(vars_)
        quantifiers = tuple(("exists", sym) for sym in variables)
    else:
        variables, free, quantified, notes = norm_internal_order(vars_, quantifiers, matrix, None)
        if free:
            # For instance finding, treat free variables as existential search
            # variables ahead of explicitly quantified ones.
            quantifiers = tuple(("exists", sym) for sym in free) + tuple(
                (q.lower(), sym) for q, sym in quantifiers
            )
        else:
            quantifiers = tuple((q.lower(), sym) for q, sym in quantifiers)
    qmap = _quantifier_map(quantifiers)
    polys = formula_polynomials(matrix) or [sp.Integer(1)]
    tower = build_collins_proj_set(polys, variables)
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
            assert parent is not None
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

    found = visit_instance_path(1, None) if variables else False
    instance = _sample_mapping(variables, found_cell) if found_cell is not None and found else None
    return LazyInstanceResult(
        instance=instance,
        status="complete",
        backend="partial-collins-lazy",
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
