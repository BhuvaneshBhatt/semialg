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
from ..formula import Formula, formula_polynomials
from ..qe.complete import evaluate_formula_on_cell, norm_internal_order


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
            notes=tuple(notes),
        )


def _quantifier_map(quantifiers: Sequence[tuple[str, sp.Symbol]]) -> dict[sp.Symbol, str]:
    out: dict[sp.Symbol, str] = {}
    for q, sym in quantifiers:
        q = q.lower()
        if q not in {"exists", "forall"}:
            raise ValueError(f"unsupported quantifier {q!r}")
        out[sym] = q
    return out


def _stack_for_parent(
    state: _TraversalState, parent: CADCell | None, level: int
) -> tuple[CADCell, ...]:
    tower = state.tower
    variables = state.variables
    if level == 1:
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
    return {sym: sample_to_expr(sample) for sym, sample in zip(variables, cell.sample, strict=True)}


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
    state = _TraversalState(tower=tower, variables=variables, matrix=matrix)
    witness: dict[str, sp.Expr] | None = None
    counterexample: dict[str, sp.Expr] | None = None

    def visit_truth_path(level: int, parent: CADCell | None) -> tuple[bool, CADCell | None]:
        nonlocal witness, counterexample
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
        stack = _stack_for_parent(state, parent, level)
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
    state = _TraversalState(tower=tower, variables=variables, matrix=matrix)
    found_cell: CADCell | None = None

    def visit_instance_path(level: int, parent: CADCell | None) -> bool:
        nonlocal found_cell
        if level > len(variables):
            assert parent is not None
            state.evaluated_leaf_cells += 1
            truth = evaluate_formula_on_cell(matrix, parent, variables)
            if truth:
                found_cell = parent
            return truth
        q = qmap.get(variables[level - 1], "exists")
        stack = _stack_for_parent(state, parent, level)
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
