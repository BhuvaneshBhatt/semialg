from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field
from typing import Literal

import sympy as sp

from ..algebraic.comparison import compare_samples
from ..formula import Formula, equational_constraints, formula_polynomials
from ..qe.complete import evaluate_formula_on_cell
from .decomposition import CompleteCAD, decomp_collins_complete, decomp_from_proj_tower
from .lifting.stack import CADCell
from .projection.lazard import build_lazard_proj_set
from .projection.mccallum import build_mccallum_proj_set
from .projection.reduced import ProjectionValidity, ReducedProjectionTower

ReducedBackend = Literal["mccallum", "lazard", "tticad"]


@dataclass(frozen=True)
class NullificationEvent:
    """A reduced-projection side-condition failure candidate.

    The event records that a polynomial used for lifting at ``level`` becomes
    identically zero over a lower-dimensional base cell. McCallum-style
    projection cannot use such a cell without extra delineating-polynomial
    repair or fallback.
    """

    level: int
    cell_index: tuple[int, ...]
    variable: sp.Symbol
    polynomial: sp.Expr
    reason: str = "polynomial nullifies over lower-dimensional cell"


@dataclass(frozen=True)
class ReducedSideReport:
    """Intrinsic reduced-CAD side-condition report.

    This is separate from the Collins-refinement certificate. A reduced CAD is
    accepted only when the intrinsic side conditions and the refinement-based
    invariant certificate both pass.
    """

    valid: bool
    checked_conditions: tuple[str, ...]
    failed_conditions: tuple[str, ...] = ()
    nullification_events: tuple[NullificationEvent, ...] = ()
    lazard_valuations: Mapping[tuple[int, tuple[int, ...], str], tuple[int, ...]] = field(
        default_factory=dict
    )


@dataclass(frozen=True)
class ReducedCertificate:
    """Post-lifting proof object for a reduced-CAD attempt.

    The certificate is deliberately conservative. A reduced CAD is accepted
    only when every cell of a full Collins refinement maps to a reduced cell and
    the requested invariant is constant on each reduced cell. For raw
    polynomial decomposition the invariant is the sign vector of the input
    polynomials. For formula decomposition the invariant is the truth value of
    the formula matrix.
    """

    valid: bool
    invariant: Literal["sign", "truth"]
    checked_refinement_cells: int
    checked_reduced_cells: int
    failures: tuple[str, ...] = ()
    grouped_values: Mapping[tuple[int, ...], tuple[object, ...]] = field(default_factory=dict)


@dataclass(frozen=True)
class SafeReducedCAD:
    """A reduced-CAD attempt with proof-carrying Collins fallback."""

    requested_backend: ReducedBackend
    effective_backend: str
    complete: bool
    used_fallback: bool
    reduced_projection: ReducedProjectionTower
    cad: CompleteCAD
    validity: ProjectionValidity
    side_conditions: ReducedSideReport | None = None
    certificate: ReducedCertificate | None = None
    fallback_cad: CompleteCAD | None = None

    @property
    def cells(self):
        return self.cad.cells


def _poly_key(poly: sp.Poly) -> str:
    return sp.sstr(sp.expand(poly.as_expr()))


def _same_sample(left, right) -> bool:
    return compare_samples(left, right) == 0


def _sample_in_cell(sample: Sequence[object], cell: CADCell) -> bool:
    if len(sample) < cell.level:
        return False
    for coord, level_cell in zip(sample, _cell_prefixes(cell), strict=True):
        left, right = level_cell.interval or (None, None)
        if level_cell.kind == "section":
            if left is None or not _same_sample(coord, left):
                return False
            continue
        if left is not None and compare_samples(coord, left) <= 0:
            return False
        if right is not None and compare_samples(coord, right) >= 0:
            return False
    return True


def _cell_prefixes(cell: CADCell) -> tuple[CADCell, ...]:
    # A final cell already carries only its own interval. The parent intervals
    # are reconstructed from the index externally in _reduced_cell_contains.
    return (cell,)


def _reduced_cell_contains(
    reduced_cells_by_level: Mapping[int, Sequence[CADCell]],
    reduced_cell: CADCell,
    sample: Sequence[object],
) -> bool:
    if len(sample) < reduced_cell.level:
        return False
    for level in range(1, reduced_cell.level + 1):
        prefix = reduced_cell.index[:level]
        try:
            level_cell = next(
                item for item in reduced_cells_by_level[level] if item.index == prefix
            )
        except StopIteration:
            return False
        coord = sample[level - 1]
        left, right = level_cell.interval or (None, None)
        if level_cell.kind == "section":
            if left is None or compare_samples(coord, left) != 0:
                return False
        else:
            if left is not None and compare_samples(coord, left) <= 0:
                return False
            if right is not None and compare_samples(coord, right) >= 0:
                return False
    return True


def _locate_reduced_cell(reduced: CompleteCAD, full_cell: CADCell) -> CADCell | None:
    candidates = reduced.cells_by_level.get(full_cell.level, tuple())
    for cell in candidates:
        if _reduced_cell_contains(reduced.cells_by_level, cell, full_cell.sample):
            return cell
    return None


def _certify_sign_invariance(
    *,
    reduced: CompleteCAD,
    refiner: CompleteCAD,
    input_polys: Sequence[sp.Expr | sp.Poly],
    variables: Sequence[sp.Symbol],
) -> ReducedCertificate:
    polys = tuple(
        poly if isinstance(poly, sp.Poly) else sp.Poly(sp.expand(poly), *variables)
        for poly in input_polys
    )
    keys = tuple(_poly_key(poly) for poly in polys)
    grouped: dict[tuple[int, ...], set[tuple[int | None, ...]]] = {}
    failures: list[str] = []
    for full_cell in refiner.cells:
        reduced_cell = _locate_reduced_cell(reduced, full_cell)
        if reduced_cell is None:
            failures.append(f"could not locate reduced cell for full cell {full_cell.index}")
            continue
        signs = tuple(full_cell.signs.get(key) for key in keys)
        grouped.setdefault(reduced_cell.index, set()).add(signs)
    for index, values in grouped.items():
        if len(values) > 1:
            failures.append(
                f"input sign vector is not invariant on reduced cell {index}: {sorted(values)!r}"
            )
    frozen = {index: tuple(values) for index, values in grouped.items()}
    return ReducedCertificate(
        valid=not failures,
        invariant="sign",
        checked_refinement_cells=len(refiner.cells),
        checked_reduced_cells=len(grouped),
        failures=tuple(failures),
        grouped_values=frozen,
    )


def certify_truth_inv(
    *,
    reduced: CompleteCAD,
    refiner: CompleteCAD,
    formula: Formula,
    variables: Sequence[sp.Symbol],
) -> ReducedCertificate:
    grouped: dict[tuple[int, ...], set[bool]] = {}
    failures: list[str] = []
    for full_cell in refiner.cells:
        reduced_cell = _locate_reduced_cell(reduced, full_cell)
        if reduced_cell is None:
            failures.append(f"could not locate reduced cell for full cell {full_cell.index}")
            continue
        grouped.setdefault(reduced_cell.index, set()).add(
            evaluate_formula_on_cell(formula, full_cell, variables)
        )
    for index, values in grouped.items():
        if len(values) > 1:
            failures.append(
                f"formula truth value is not invariant on reduced cell {index}: {sorted(values)!r}"
            )
    frozen = {index: tuple(values) for index, values in grouped.items()}
    return ReducedCertificate(
        valid=not failures,
        invariant="truth",
        checked_refinement_cells=len(refiner.cells),
        checked_reduced_cells=len(grouped),
        failures=tuple(failures),
        grouped_values=frozen,
    )


def proj_validity_from_cert(
    original: ProjectionValidity,
    certificate: ReducedCertificate,
    *,
    backend: ReducedBackend,
) -> ProjectionValidity:
    if certificate.valid:
        return ProjectionValidity(
            theory=backend,
            valid=True,
            complete_if_used=True,
            reason=(
                f"active {backend} reduced CAD certified by Collins-refinement "
                f"{certificate.invariant}-invariance check"
            ),
            checked_conditions=tuple(
                dict.fromkeys(
                    (
                        *original.checked_conditions,
                        "reduced lifting",
                        "Collins refinement certificate",
                    )
                )
            ),
            failed_conditions=(),
            fallback_backend=None,
            details={**dict(original.details), "certificate": certificate},
        )
    return ProjectionValidity(
        theory=backend,
        valid=False,
        complete_if_used=False,
        reason=f"active {backend} reduced CAD failed certification; using Collins fallback",
        checked_conditions=tuple(
            dict.fromkeys(
                (*original.checked_conditions, "reduced lifting", "Collins refinement certificate")
            )
        ),
        failed_conditions=tuple(
            dict.fromkeys((*original.failed_conditions, *certificate.failures))
        ),
        fallback_backend="collins-complete",
        details={**dict(original.details), "certificate": certificate},
    )


def build_reduced_proj(
    polys: Sequence[sp.Expr | sp.Poly],
    variables: Sequence[sp.Symbol],
    *,
    backend: ReducedBackend,
    equational_constraints: Sequence[sp.Expr],
    certify_reduced: bool,
) -> ReducedProjectionTower:
    if backend == "mccallum":
        return build_mccallum_proj_set(
            polys,
            variables,
            equational_constraints=equational_constraints,
            certify=certify_reduced,
        )
    if backend == "lazard":
        return build_lazard_proj_set(
            polys,
            variables,
            equational_constraints=equational_constraints,
            certify=certify_reduced,
        )
    raise ValueError(f"unknown reduced CAD backend: {backend!r}")


def _specialize_with_lazard(
    poly: sp.Poly, variables: Sequence[sp.Symbol], sample: Sequence[object], level: int
) -> tuple[sp.Expr, tuple[int, ...]]:
    from .lifting.lazard import lazard_evaluate

    assigned = variables[: level - 1]
    values = [
        sp.expand(
            getattr(item, "expr", None) or item.as_expr() if hasattr(item, "as_expr") else item
        )
        for item in ()
    ]
    values = []
    from ..algebraic.samples import sample_to_expr

    for item in sample[: level - 1]:
        values.append(sample_to_expr(item))
    result = lazard_evaluate(poly.as_expr(), assigned, values)
    return sp.expand(result.final_expr), result.valuation


def _specialized_expr(
    poly: sp.Poly,
    variables: Sequence[sp.Symbol],
    sample: Sequence[object],
    level: int,
    *,
    backend: ReducedBackend,
) -> tuple[sp.Expr, tuple[int, ...]]:
    from ..algebraic.samples import sample_to_expr

    if backend == "lazard":
        return _specialize_with_lazard(poly, variables, sample, level)
    substitutions = {variables[i]: sample_to_expr(sample[i]) for i in range(level - 1)}
    return sp.expand(poly.as_expr().subs(substitutions)), ()


def _scan_reduced_conditions(cad: CompleteCAD, *, backend: ReducedBackend) -> ReducedSideReport:
    variables = cad.tower.variables
    events: list[NullificationEvent] = []
    valuations: dict[tuple[int, tuple[int, ...], str], tuple[int, ...]] = {}
    checked = ["well-orientedness/nullification scan"]
    if backend == "lazard":
        checked.append("Lazard valuation-aware specialization")
    for level in range(2, len(variables) + 1):
        parent_cells = cad.cells_by_level.get(level - 1, ())
        level_polys = cad.tower.level(level).polynomials
        variable = variables[level - 1]
        for parent in parent_cells:
            for poly in level_polys:
                expr, valuation = _specialized_expr(
                    poly, variables, parent.sample, level, backend=backend
                )
                if backend == "lazard" and valuation:
                    valuations[(level, parent.index, _poly_key(poly))] = valuation
                if sp.expand(expr) == 0:
                    # Lazard cancellation can turn a nullified expression into a
                    # nonzero quotient. Only report an event if the final
                    # Lazard-evaluated expression is still zero.
                    events.append(
                        NullificationEvent(
                            level=level,
                            cell_index=parent.index,
                            variable=variable,
                            polynomial=poly.as_expr(),
                        )
                    )
    failed: tuple[str, ...] = tuple(
        f"nullification at level {event.level} over cell {event.cell_index}: {sp.sstr(event.polynomial)}"
        for event in events
    )
    return ReducedSideReport(
        valid=not events,
        checked_conditions=tuple(checked),
        failed_conditions=failed,
        nullification_events=tuple(events),
        lazard_valuations=valuations,
    )


def _validity_with_condition(
    original: ProjectionValidity, side: ReducedSideReport, *, backend: ReducedBackend
) -> ProjectionValidity:
    if side.valid:
        return ProjectionValidity(
            theory=backend,
            valid=original.valid,
            complete_if_used=original.complete_if_used,
            reason=original.reason,
            checked_conditions=tuple(
                dict.fromkeys((*original.checked_conditions, *side.checked_conditions))
            ),
            failed_conditions=original.failed_conditions,
            fallback_backend=original.fallback_backend,
            details={**dict(original.details), "side_conditions": side},
        )
    return ProjectionValidity(
        theory=backend,
        valid=False,
        complete_if_used=False,
        reason=f"{backend} side conditions failed; using Collins fallback",
        checked_conditions=tuple(
            dict.fromkeys((*original.checked_conditions, *side.checked_conditions))
        ),
        failed_conditions=tuple(
            dict.fromkeys((*original.failed_conditions, *side.failed_conditions))
        ),
        fallback_backend="collins-complete",
        details={**dict(original.details), "side_conditions": side},
    )


def _dedupe_polys(polys: Sequence[sp.Poly]) -> tuple[sp.Poly, ...]:
    seen: set[str] = set()
    out: list[sp.Poly] = []
    for poly in polys:
        key = _poly_key(poly)
        if key not in seen and not poly.is_zero and poly.total_degree() > 0:
            seen.add(key)
            out.append(poly)
    return tuple(out)


def delin_repair_polys(
    event: NullificationEvent, variables: Sequence[sp.Symbol]
) -> tuple[sp.Poly, ...]:
    """Return lower-level delineating polynomials for one nullification event."""

    level = event.level
    if level <= 1:
        return tuple()
    lower_gens = tuple(variables[: level - 1])
    var = event.variable
    expr = sp.expand(event.polynomial)
    repair: list[sp.Poly] = []
    try:
        as_univar = sp.Poly(expr, var)
        for coeff in as_univar.all_coeffs():
            coeff = sp.expand(coeff)
            if coeff != 0 and lower_gens:
                try:
                    repair.append(sp.Poly(coeff, *lower_gens))
                except Exception:
                    pass
        if as_univar.degree() > 1 and lower_gens:
            for candidate in (
                sp.discriminant(expr, var),
                sp.resultant(expr, sp.diff(expr, var), var),
            ):
                candidate = sp.expand(candidate)
                if candidate != 0:
                    try:
                        repair.append(sp.Poly(candidate, *lower_gens))
                    except Exception:
                        pass
    except Exception:
        pass
    return _dedupe_polys(repair)


def _augment_tower_repairs(
    projection: ReducedProjectionTower,
    side_conditions: ReducedSideReport,
    variables: Sequence[sp.Symbol],
) -> ReducedProjectionTower | None:
    """Add local delineating factors for failed side-condition cells."""

    if not side_conditions.nullification_events:
        return None
    from .projection.collins import ProjectionLevel, ProjectionPolynomial, ProjectionTower

    additions: dict[int, list[sp.Poly]] = {}
    for event in side_conditions.nullification_events:
        target_level = max(1, event.level - 1)
        additions.setdefault(target_level, []).extend(delin_repair_polys(event, variables))
    if not any(additions.values()):
        return None
    new_levels = []
    added_by_level: dict[int, int] = {}
    for level in projection.tower.levels:
        extra = _dedupe_polys(additions.get(level.level, ()))
        if not extra:
            new_levels.append(level)
            continue
        existing = list(level.polynomials)
        combined = _dedupe_polys(tuple(existing) + extra)
        entries = list(level.entries)
        existing_keys = {_poly_key(poly) for poly in existing}
        for poly in extra:
            if _poly_key(poly) not in existing_keys:
                entries.append(
                    ProjectionPolynomial(
                        poly=poly,
                        level=level.level,
                        source="delineating-repair",
                        parents=(),
                        operation_variable=level.variable,
                        expression=poly.as_expr(),
                    )
                )
        added_by_level[level.level] = len(combined) - len(existing)
        new_levels.append(ProjectionLevel(level.level, level.variable, combined, tuple(entries)))
    if not added_by_level:
        return None
    tower = ProjectionTower(
        variables=projection.tower.variables,
        levels=tuple(new_levels),
        original_polynomials=projection.tower.original_polynomials,
        metadata={
            **dict(projection.tower.metadata),
            "delineating_repair": True,
            "repair_added_by_level": added_by_level,
        },
    )
    validity = ProjectionValidity(
        theory=projection.requested_theory,
        valid=False,
        complete_if_used=False,
        reason="reduced projection locally augmented with delineating repair factors; certification pending",
        checked_conditions=tuple(
            dict.fromkeys(
                (*projection.validity.checked_conditions, "localized delineating-polynomial repair")
            )
        ),
        failed_conditions=("repair certificate pending",),
        fallback_backend="collins-complete",
        details={**dict(projection.validity.details), "repair_added_by_level": added_by_level},
    )
    return ReducedProjectionTower(projection.requested_theory, tower, validity)


def _attempt_reduced_certifi(
    *,
    reduced_projection: ReducedProjectionTower,
    polys: Sequence[sp.Expr | sp.Poly],
    variables: Sequence[sp.Symbol],
    backend: ReducedBackend,
    formula: Formula | None = None,
    collins_cad: CompleteCAD | None = None,
) -> tuple[CompleteCAD, ReducedSideReport, ReducedCertificate, ProjectionValidity]:
    reduced_cad = decomp_from_proj_tower(reduced_projection.tower, backend=f"{backend}-reduced")
    side_conditions = _scan_reduced_conditions(reduced_cad, backend=backend)
    collins_cad = collins_cad or decomp_collins_complete(polys, variables)
    if formula is None:
        certificate = _certify_sign_invariance(
            reduced=reduced_cad,
            refiner=collins_cad,
            input_polys=polys,
            variables=variables,
        )
    else:
        certificate = certify_truth_inv(
            reduced=reduced_cad,
            refiner=collins_cad,
            formula=formula,
            variables=variables,
        )
    certificate_validity = proj_validity_from_cert(
        reduced_projection.validity, certificate, backend=backend
    )
    validity = _validity_with_condition(certificate_validity, side_conditions, backend=backend)
    return reduced_cad, side_conditions, certificate, validity


def decompose_reduced_safe(
    polys: Sequence[sp.Expr | sp.Poly],
    variables: Sequence[sp.Symbol],
    *,
    backend: ReducedBackend,
    equational_constraints: Sequence[sp.Expr] = (),
    certify_reduced: bool = True,
) -> SafeReducedCAD:
    """Attempt active reduced CAD, then localized repair before fallback."""

    vars_tuple = tuple(variables)
    reduced_projection = build_reduced_proj(
        polys,
        vars_tuple,
        backend=backend,
        equational_constraints=equational_constraints,
        certify_reduced=certify_reduced,
    )
    collins_cad = decomp_collins_complete(polys, vars_tuple)
    reduced_cad, side_conditions, certificate, validity = _attempt_reduced_certifi(
        reduced_projection=reduced_projection,
        polys=polys,
        variables=vars_tuple,
        backend=backend,
        collins_cad=collins_cad,
    )
    if certificate.valid and side_conditions.valid:
        return SafeReducedCAD(
            requested_backend=backend,
            effective_backend=f"{backend}-reduced-certified",
            complete=True,
            used_fallback=False,
            reduced_projection=reduced_projection,
            cad=reduced_cad,
            validity=validity,
            side_conditions=side_conditions,
            certificate=certificate,
            fallback_cad=collins_cad,
        )

    repaired_projection = _augment_tower_repairs(reduced_projection, side_conditions, vars_tuple)
    if repaired_projection is not None:
        repaired_cad, repaired_side, repaired_cert, repaired_validity = _attempt_reduced_certifi(
            reduced_projection=repaired_projection,
            polys=polys,
            variables=vars_tuple,
            backend=backend,
            collins_cad=collins_cad,
        )
        if repaired_cert.valid and repaired_side.valid:
            repaired_validity = ProjectionValidity(
                theory=backend,
                valid=True,
                complete_if_used=True,
                reason=f"{backend} accepted after localized delineating-polynomial repair and certification",
                checked_conditions=tuple(
                    dict.fromkeys(
                        (*repaired_validity.checked_conditions, "localized repair accepted")
                    )
                ),
                failed_conditions=(),
                fallback_backend=None,
                details={
                    **dict(repaired_validity.details),
                    "original_side_conditions": side_conditions,
                },
            )
            return SafeReducedCAD(
                requested_backend=backend,
                effective_backend=f"{backend}-repaired-certified",
                complete=True,
                used_fallback=False,
                reduced_projection=repaired_projection,
                cad=repaired_cad,
                validity=repaired_validity,
                side_conditions=repaired_side,
                certificate=repaired_cert,
                fallback_cad=collins_cad,
            )

    return SafeReducedCAD(
        requested_backend=backend,
        effective_backend="collins-complete",
        complete=True,
        used_fallback=True,
        reduced_projection=reduced_projection,
        cad=collins_cad,
        validity=validity,
        side_conditions=side_conditions,
        certificate=certificate,
        fallback_cad=collins_cad,
    )


def decomp_form_reduced_safe(
    formula: Formula,
    variables: Sequence[sp.Symbol],
    *,
    backend: ReducedBackend,
    certify_reduced: bool = True,
) -> SafeReducedCAD:
    """Formula-oriented reduced CAD with truth-invariance and local repair."""

    vars_tuple = tuple(variables)
    polys = formula_polynomials(formula)
    reduced_projection = build_reduced_proj(
        polys,
        vars_tuple,
        backend=backend,
        equational_constraints=equational_constraints(formula),
        certify_reduced=certify_reduced,
    )
    collins_cad = decomp_collins_complete(polys, vars_tuple)
    reduced_cad, side_conditions, certificate, validity = _attempt_reduced_certifi(
        reduced_projection=reduced_projection,
        polys=polys,
        variables=vars_tuple,
        backend=backend,
        formula=formula,
        collins_cad=collins_cad,
    )
    if certificate.valid and side_conditions.valid:
        return SafeReducedCAD(
            requested_backend=backend,
            effective_backend=f"{backend}-reduced-certified",
            complete=True,
            used_fallback=False,
            reduced_projection=reduced_projection,
            cad=reduced_cad,
            validity=validity,
            side_conditions=side_conditions,
            certificate=certificate,
            fallback_cad=collins_cad,
        )

    repaired_projection = _augment_tower_repairs(reduced_projection, side_conditions, vars_tuple)
    if repaired_projection is not None:
        repaired_cad, repaired_side, repaired_cert, repaired_validity = _attempt_reduced_certifi(
            reduced_projection=repaired_projection,
            polys=polys,
            variables=vars_tuple,
            backend=backend,
            formula=formula,
            collins_cad=collins_cad,
        )
        if repaired_cert.valid and repaired_side.valid:
            repaired_validity = ProjectionValidity(
                theory=backend,
                valid=True,
                complete_if_used=True,
                reason=f"{backend} formula CAD accepted after localized delineating-polynomial repair and certification",
                checked_conditions=tuple(
                    dict.fromkeys(
                        (*repaired_validity.checked_conditions, "localized repair accepted")
                    )
                ),
                failed_conditions=(),
                fallback_backend=None,
                details={
                    **dict(repaired_validity.details),
                    "original_side_conditions": side_conditions,
                },
            )
            return SafeReducedCAD(
                requested_backend=backend,
                effective_backend=f"{backend}-repaired-certified",
                complete=True,
                used_fallback=False,
                reduced_projection=repaired_projection,
                cad=repaired_cad,
                validity=repaired_validity,
                side_conditions=repaired_side,
                certificate=repaired_cert,
                fallback_cad=collins_cad,
            )

    return SafeReducedCAD(
        requested_backend=backend,
        effective_backend="collins-complete",
        complete=True,
        used_fallback=True,
        reduced_projection=reduced_projection,
        cad=collins_cad,
        validity=validity,
        side_conditions=side_conditions,
        certificate=certificate,
        fallback_cad=collins_cad,
    )


__all__ = [
    "NullificationEvent",
    "ReducedBackend",
    "ReducedCertificate",
    "ReducedSideReport",
    "SafeReducedCAD",
    "decomp_form_reduced_safe",
    "decompose_reduced_safe",
]
