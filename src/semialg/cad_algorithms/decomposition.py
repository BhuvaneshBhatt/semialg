from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field
from time import perf_counter

import sympy as sp

from .._immutable import freeze_mapping
from ..algebraic.comparison import sort_samples
from ..algebraic.roots import isolate_real_roots
from ..algebraic.sample_points import choose_sector_sample
from ..algebraic.samples import AlgebraicRoot, Sample, sample_to_expr
from ..algebraic.signs import sign_at_sample
from ..errors import VarietyCADUnsupported
from .lifting.groebner import GroebnerLiftOps, lift_groebner_variety
from .lifting.sign_invariance import SignInvarianceCheck, verify_cad_sign_inv
from .lifting.stack import CADCell, sign_table
from .performance_cache import COMPLETE_CADS, STATS
from .polynomial_utils import exact_univariate_poly, polynomial_family_key
from .polynomial_utils import polynomial_key as _poly_key
from .projection.collins import ProjectionTower, build_collins_proj_set
from .projection.groebner import GroebnerVarietyProjection


@dataclass(frozen=True)
class CADDiagnostics:
    cell_count_by_level: Mapping[int, int]
    proj_poly_count_by_level: Mapping[int, int]
    max_stack_size: int
    timing_by_stage: Mapping[str, float] = field(default_factory=dict)
    invariant_failures: tuple[str, ...] = ()

    def __post_init__(self) -> None:
        object.__setattr__(self, "cell_count_by_level", freeze_mapping(self.cell_count_by_level))
        object.__setattr__(
            self, "proj_poly_count_by_level", freeze_mapping(self.proj_poly_count_by_level)
        )
        object.__setattr__(self, "timing_by_stage", freeze_mapping(self.timing_by_stage))


@dataclass(frozen=True)
class CompleteCAD:
    tower: ProjectionTower
    cells_by_level: Mapping[int, tuple[CADCell, ...]]
    complete: bool = True
    backend: str = "collins-complete"
    diagnostics: CADDiagnostics | None = None

    def __post_init__(self) -> None:
        object.__setattr__(self, "cells_by_level", freeze_mapping(self.cells_by_level))

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
            univar = exact_univariate_poly(expr, var, algebraic_extension=False)
        except (sp.PolynomialError, TypeError, ValueError):
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
    """Return the original tower polynomial defining a section sample.

    Testing the unspecialized polynomial with ``sign_at_sample`` preserves its
    provenance even when the base coordinates are algebraic.  Direct SymPy
    substitution can otherwise leave an opaque ordered-root expression that
    fails to simplify to zero, causing the section to be labelled only by its
    specialized fiber polynomial.
    """

    fiber_var = variables[level - 1]
    for poly in level_polys:
        if fiber_var not in poly.as_expr().free_symbols:
            continue
        if sign_at_sample(poly, sample) == 0:
            return poly.as_expr(), _poly_key(poly)
    return None, None


def _build_stack(
    parent: CADCell | None, roots: Sequence[Sample], level: int, tower: ProjectionTower
) -> tuple[CADCell, ...]:
    """Lift one CAD stack over a parent cell from the ordered exact fiber roots."""
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


def decomp_groebner_variety(
    projection: GroebnerVarietyProjection,
    *,
    lift_ops: GroebnerLiftOps | None = None,
) -> CompleteCAD:
    """Lift a section-only CAD of a certified zero-dimensional equality variety.

    The result is complete *on the common equality variety*, not a partition of
    all of real space.  It is therefore intended for formula decomposition when
    those equalities are necessary conditions for truth.
    """

    cells_by_level, timings = lift_groebner_variety(projection, ops=lift_ops)
    tower = projection.tower
    diagnostics = CADDiagnostics(
        cell_count_by_level={level: len(cells) for level, cells in cells_by_level.items()},
        proj_poly_count_by_level=tower.poly_count_by_level(),
        max_stack_size=CompleteCAD(tower=tower, cells_by_level=cells_by_level).max_stack_size(),
        timing_by_stage=timings,
        invariant_failures=(),
    )
    return CompleteCAD(
        tower=tower,
        cells_by_level=cells_by_level,
        complete=False,
        backend="groebner-variety",
        diagnostics=diagnostics,
    )


def try_decomp_groebner_variety(
    projection: GroebnerVarietyProjection,
    *,
    lift_ops: GroebnerLiftOps | None = None,
) -> CompleteCAD | None:
    """Return a certified variety CAD, or ``None`` when exact lifting declines."""

    try:
        return decomp_groebner_variety(projection, lift_ops=lift_ops)
    except VarietyCADUnsupported:
        return None


def decomp_collins_complete(
    polys: Sequence[sp.Expr | sp.Poly], variables: Sequence[sp.Symbol]
) -> CompleteCAD:
    """Construct a complete Collins CAD with exact projection, lifting, and diagnostic metadata."""
    vars_tuple = tuple(variables)
    if not vars_tuple:
        tower = ProjectionTower(
            variables=tuple(),
            levels=tuple(),
            original_polynomials=tuple(),
            metadata={"projection": "collins-complete", "zero_variable": True},
        )
        cell = CADCell(level=0, index=tuple(), sample=tuple(), kind="sector")
        diagnostics = CADDiagnostics(
            cell_count_by_level={0: 1},
            proj_poly_count_by_level={},
            max_stack_size=1,
            timing_by_stage={},
            invariant_failures=(),
        )
        return CompleteCAD(
            tower=tower,
            cells_by_level={0: (cell,)},
            backend="collins-complete",
            diagnostics=diagnostics,
        )
    try:
        input_polys = tuple(
            poly
            if isinstance(poly, sp.Poly) and tuple(poly.gens) == vars_tuple
            else sp.Poly(poly.as_expr() if isinstance(poly, sp.Poly) else poly, *vars_tuple)
            for poly in polys
        )
        cache_key = (polynomial_family_key(input_polys), vars_tuple)
    except (sp.PolynomialError, TypeError, ValueError):
        cache_key = None
    if cache_key is not None:
        cached = COMPLETE_CADS.get(cache_key)
        if cached is not None:
            STATS.complete_cad_hits += 1
            return cached  # type: ignore[return-value]
        STATS.complete_cad_misses += 1

    start = perf_counter()
    tower = build_collins_proj_set(polys, vars_tuple)
    projection_time = perf_counter() - start
    cad = decomp_from_proj_tower(tower, backend="collins-complete")
    timings = dict(cad.diagnostics.timing_by_stage if cad.diagnostics is not None else {})
    timings["projection"] = projection_time
    diagnostics = CADDiagnostics(
        cell_count_by_level=cad.cell_count_by_level(),
        proj_poly_count_by_level=tower.poly_count_by_level(),
        max_stack_size=cad.max_stack_size(),
        timing_by_stage=timings,
        invariant_failures=(
            cad.diagnostics.invariant_failures if cad.diagnostics is not None else ()
        ),
    )
    result = CompleteCAD(tower=tower, cells_by_level=cad.cells_by_level, diagnostics=diagnostics)
    if cache_key is not None:
        COMPLETE_CADS.put(cache_key, result)
    return result
