from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field
from time import perf_counter

import sympy as sp

from ..algebraic.comparison import sort_samples
from ..algebraic.roots import isolate_real_roots
from ..algebraic.sample_points import choose_sector_sample
from ..algebraic.samples import AlgebraicRoot, Sample, sample_to_expr
from .lifting.sign_invariance import SignInvarianceCheck, verify_cad_sign_inv
from .lifting.stack import CADCell, sign_table
from .projection.collins import ProjectionTower, build_collins_proj_set


@dataclass(frozen=True)
class CADDiagnostics:
    cell_count_by_level: Mapping[int, int]
    proj_poly_count_by_level: Mapping[int, int]
    max_stack_size: int
    timing_by_stage: Mapping[str, float] = field(default_factory=dict)
    invariant_failures: tuple[str, ...] = ()


@dataclass(frozen=True)
class CompleteCAD:
    tower: ProjectionTower
    cells_by_level: dict[int, tuple[CADCell, ...]]
    complete: bool = True
    backend: str = "collins-complete"
    diagnostics: CADDiagnostics | None = None

    @property
    def cells(self) -> tuple[CADCell, ...]:
        return self.cells_by_level.get(len(self.tower.variables), tuple())

    def cell_count_by_level(self) -> dict[int, int]:
        return {level: len(cells) for level, cells in self.cells_by_level.items()}

    def proj_poly_count_by_level(self) -> dict[int, int]:
        return self.tower.poly_count_by_level()

    def max_stack_size(self) -> int:
        max_size = 0
        for level, cells in self.cells_by_level.items():
            if level == 1:
                max_size = max(max_size, len(cells))
                continue
            counts: dict[tuple[int, ...] | None, int] = {}
            for cell in cells:
                counts[cell.parent_index] = counts.get(cell.parent_index, 0) + 1
            if counts:
                max_size = max(max_size, max(counts.values()))
        return max_size

    def verify_sign_invariance(self) -> SignInvarianceCheck:
        return verify_cad_sign_inv(self.cells_by_level, self.tower)

    def failed_invariants(self) -> tuple[str, ...]:
        return self.verify_sign_invariance().failures


def _poly_key(poly: sp.Poly) -> str:
    return sp.sstr(sp.expand(poly.as_expr()))


def _stack_roots_over_point(
    polys: Sequence[sp.Poly],
    variables: Sequence[sp.Symbol],
    prefix: Sequence[Sample],
    var: sp.Symbol,
    *,
    use_lazard: bool = False,
) -> tuple[AlgebraicRoot, ...]:
    roots: list[AlgebraicRoot] = []
    substitutions = {variables[i]: sample_to_expr(prefix[i]) for i in range(len(prefix))}
    sample_prefix = [sample_to_expr(prefix[i]) for i in range(len(prefix))]
    for poly in polys:
        if use_lazard and prefix:
            from .lifting.lazard import lazard_evaluate

            expr = sp.expand(
                lazard_evaluate(poly.as_expr(), variables[: len(prefix)], sample_prefix).final_expr
            )
        else:
            expr = sp.expand(poly.as_expr().subs(substitutions))
        if expr == 0:
            continue
        try:
            univar = sp.Poly(expr, var, domain="EX")
        except Exception:
            continue
        if univar.degree() > 0:
            roots.extend(isolate_real_roots(univar))
    return tuple(root for root in sort_samples(tuple(roots)) if isinstance(root, AlgebraicRoot))


def _section_poly_for_root(
    level_polys: Sequence[sp.Poly],
    variables: Sequence[sp.Symbol],
    sample: Sequence[Sample],
    level: int,
) -> tuple[sp.Expr | None, str | None]:
    substitutions = {variables[i]: sample_to_expr(sample[i]) for i in range(level)}
    for poly in level_polys:
        if sp.simplify(poly.as_expr().subs(substitutions)) == 0:
            return poly.as_expr(), _poly_key(poly)
    return None, None


def _build_stack(
    parent: CADCell | None, roots: Sequence[Sample], level: int, tower: ProjectionTower
) -> tuple[CADCell, ...]:
    variables = tower.variables
    prefix = parent.sample if parent is not None else tuple()
    parent_index = parent.index if parent is not None else tuple()
    level_polys = tower.level(level).polynomials
    cells: list[CADCell] = []
    bounds: list[Sample | None] = [None, *roots, None]
    position = 0
    for pos in range(len(bounds) - 1):
        left = bounds[pos]
        right = bounds[pos + 1]
        sector_sample = choose_sector_sample(left, right)
        sample = (*prefix, sector_sample)
        cells.append(
            CADCell(
                level=level,
                index=(*parent_index, position),
                sample=sample,
                interval=(left, right),
                signs=sign_table(level_polys, sample),
                kind="sector",
                parent_index=parent.index if parent is not None else None,
                stack_position=position,
            )
        )
        position += 1
        if right is not None:
            sample = (*prefix, right)
            section_expr, section_key = _section_poly_for_root(
                level_polys, variables, sample, level
            )
            if section_expr is None:
                section_expr = (
                    getattr(right, "polynomial", None).as_expr()
                    if isinstance(right, AlgebraicRoot)
                    else sp.Integer(0)
                )
            cells.append(
                CADCell(
                    level=level,
                    index=(*parent_index, position),
                    sample=sample,
                    interval=(right, right),
                    section_polynomial=section_expr,
                    signs=sign_table(level_polys, sample),
                    kind="section",
                    parent_index=parent.index if parent is not None else None,
                    stack_position=position,
                    root_index=pos,
                    defining_polynomial_key=section_key,
                )
            )
            position += 1
    return tuple(cells)


def _make_diagnostics(
    tower: ProjectionTower,
    cells_by_level: dict[int, tuple[CADCell, ...]],
    timings: Mapping[str, float],
) -> CADDiagnostics:
    check = verify_cad_sign_inv(cells_by_level, tower)
    temp = CompleteCAD(tower=tower, cells_by_level=cells_by_level)
    return CADDiagnostics(
        cell_count_by_level=temp.cell_count_by_level(),
        proj_poly_count_by_level=tower.poly_count_by_level(),
        max_stack_size=temp.max_stack_size(),
        timing_by_stage=dict(timings),
        invariant_failures=check.failures,
    )


def decomp_from_proj_tower(tower: ProjectionTower, *, backend: str | None = None) -> CompleteCAD:
    """Lift a CAD from an already-built projection tower.

    This is the common lifting engine used by the complete Collins path and by
    safe reduced-CAD attempts. Correctness is not inferred from this function
    alone: callers using reduced towers must attach a separate proof/fallback
    certificate.
    """

    timings: dict[str, float] = {"projection": 0.0}
    variables = tower.variables
    use_lazard = str(tower.metadata.get("projection", "")) == "lazard"
    cells_by_level: dict[int, tuple[CADCell, ...]] = {}
    if not variables:
        diagnostics = _make_diagnostics(tower, {0: tuple()}, timings)
        return CompleteCAD(
            tower=tower,
            cells_by_level={0: tuple()},
            backend=backend or str(tower.metadata.get("projection", "custom")),
            diagnostics=diagnostics,
        )
    lift_start = perf_counter()
    roots: list[AlgebraicRoot] = []
    for poly in tower.level(1).polynomials:
        roots.extend(isolate_real_roots(poly))
    cells_by_level[1] = _build_stack(None, sort_samples(tuple(roots)), 1, tower)
    for level in range(2, len(variables) + 1):
        next_cells: list[CADCell] = []
        for parent in cells_by_level[level - 1]:
            stack_roots = _stack_roots_over_point(
                tower.level(level).polynomials,
                variables,
                parent.sample,
                variables[level - 1],
                use_lazard=use_lazard,
            )
            next_cells.extend(_build_stack(parent, stack_roots, level, tower))
        cells_by_level[level] = tuple(next_cells)
    timings["lifting"] = perf_counter() - lift_start
    diag_start = perf_counter()
    diagnostics = _make_diagnostics(tower, cells_by_level, timings)
    timings["diagnostics"] = perf_counter() - diag_start
    diagnostics = CADDiagnostics(
        cell_count_by_level=diagnostics.cell_count_by_level,
        proj_poly_count_by_level=tower.poly_count_by_level(),
        max_stack_size=diagnostics.max_stack_size,
        timing_by_stage=timings,
        invariant_failures=diagnostics.invariant_failures,
    )
    return CompleteCAD(
        tower=tower,
        cells_by_level=cells_by_level,
        backend=backend or str(tower.metadata.get("projection", "custom")),
        diagnostics=diagnostics,
    )


def decomp_collins_complete(
    polys: Sequence[sp.Expr | sp.Poly], variables: Sequence[sp.Symbol]
) -> CompleteCAD:
    start = perf_counter()
    tower = build_collins_proj_set(polys, variables)
    projection_time = perf_counter() - start
    cad = decomp_from_proj_tower(tower, backend="collins-complete")
    timings = dict(cad.diagnostics.timing_by_stage if cad.diagnostics is not None else {})
    timings["projection"] = projection_time
    diagnostics = CADDiagnostics(
        cell_count_by_level=cad.cell_count_by_level(),
        proj_poly_count_by_level=tower.poly_count_by_level(),
        max_stack_size=cad.max_stack_size(),
        timing_by_stage=timings,
        invariant_failures=cad.failed_invariants(),
    )
    return CompleteCAD(tower=tower, cells_by_level=cad.cells_by_level, diagnostics=diagnostics)
