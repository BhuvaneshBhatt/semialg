from __future__ import annotations

import contextlib
from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field
from typing import Literal

import sympy as sp

from ..algebraic.samples import FiberRootContext
from ..algebraic.signs import sign_at_sample
from ..errors import ExactEvaluationFailure
from ..formula import Formula, equational_constraints, formula_polynomials, to_sympy
from ..qe.complete import evaluate_formula_on_cell
from .decomposition import CompleteCAD, decomp_collins_complete, decomp_from_proj_tower
from .invariant_scope import _missing_boundaries, _with_invariant
from .lifting.stack import CADCell
from .polynomial_utils import polynomial_key as _poly_key
from .polynomial_utils import univariate_leading_coefficient
from .projection.ec_hypotheses import ec_parent_admissible, replay_ec_projection
from .projection.lazard import build_lazard_proj_set
from .projection.mccallum import build_mccallum_proj_set
from .projection.reduced import ProjectionValidity, ReducedProjectionTower, _candidate_ec

ReducedBackend = Literal["mccallum", "lazard"]


@dataclass(frozen=True)
class NullificationEvent:
    """A reduced-projection side-condition failure candidate.

    The event records that a lifting polynomial becomes identically zero at a
    base sample. Extending that observation to the whole cell requires the
    relevant coefficient invariants. Reduced lifting needs repair or fallback
    when its nullification conditions cannot be established.
    """

    level: int
    cell_index: tuple[int, ...]
    variable: sp.Symbol
    polynomial: sp.Expr
    reason: str = "polynomial nullifies at the base sample"


@dataclass(frozen=True)
class ReducedSideReport:
    """Sample-level nullification and valuation diagnostics for reduced CAD.

    The safe driver combines these lifting checks with an intrinsic theorem
    certificate or a Collins-refinement certificate for the requested invariant.
    The sample checks alone do not establish degree invariance on a base cell.
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

    An intrinsic certificate records the reduced theorem's side conditions.
    Otherwise every cell of a complete Collins refinement must map to a reduced
    cell, with a constant requested invariant on each reduced cell. Polynomial
    decomposition uses the input sign vector; formula decomposition uses the
    matrix truth value. Refinement certificates report the number of checked
    reference cells; intrinsic certificates report zero.
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


def _locate_reduced_cell(reduced: CompleteCAD, full_cell: CADCell) -> CADCell | None:
    try:
        return reduced._point_locator.locate(full_cell.sample)
    except (ExactEvaluationFailure, ArithmeticError, ValueError, IndexError, NotImplementedError):
        # A failed exact boundary evaluation cannot establish refinement
        # containment. The certificate records a failure and requests fallback.
        return None


def _refinement_failures(reduced, refiner, variables, polynomials) -> tuple[str, ...]:
    """Require an ambient sign-invariant reference covering the requested boundaries.

    Agreement at sample points proves an invariant only when the reference CAD
    already separates every polynomial factor used by that invariant. A finite
    variety or a decomposition for an unrelated family cannot supply this proof.
    """

    variables = tuple(variables)
    if reduced.tower.variables != variables or refiner.tower.variables != variables:
        return ("refinement variable order does not match the requested invariant",)
    if not refiner.complete or not refiner.tower.metadata.get("complete", False):
        return ("reference must be a complete ambient sign-invariant CAD",)
    if refiner.diagnostics is not None and refiner.diagnostics.invariant_failures:
        return tuple(refiner.diagnostics.invariant_failures)
    missing = _missing_boundaries(refiner.tower, polynomials)
    if missing:
        return (f"reference omits polynomial boundaries: {missing!r}",)
    return ()


def _certify_sign_invariance(
    *,
    reduced: CompleteCAD,
    refiner: CompleteCAD,
    input_polys: Sequence[sp.Expr | sp.Poly],
    variables: Sequence[sp.Symbol],
) -> ReducedCertificate:
    reference_failures = _refinement_failures(reduced, refiner, variables, input_polys)
    if reference_failures:
        return ReducedCertificate(False, "sign", 0, 0, reference_failures)
    polys = tuple(
        poly if isinstance(poly, sp.Poly) else sp.Poly(sp.expand(poly), *variables)
        for poly in input_polys
    )
    keys = tuple(_poly_key(poly) for poly in polys)
    grouped: dict[tuple[int, ...], set[tuple[int, ...]]] = {}
    failures: list[str] = []
    for full_cell in refiner.cells:
        reduced_cell = _locate_reduced_cell(reduced, full_cell)
        if reduced_cell is None:
            failures.append(f"could not locate reduced cell for full cell {full_cell.index}")
            continue
        signs = tuple(
            full_cell.signs[key]
            if key in full_cell.signs
            else sign_at_sample(poly, full_cell.sample)
            for key, poly in zip(keys, polys, strict=True)
        )
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
    """Check truth invariance against a complete CAD for the formula polynomials."""

    reference_failures = _refinement_failures(
        reduced, refiner, variables, formula_polynomials(formula)
    )
    if reference_failures:
        return ReducedCertificate(False, "truth", 0, 0, reference_failures)
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
    """Merge a reduced-CAD refinement certificate into the projection-validity record."""
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
    if backend == "tticad":
        raise ValueError(
            "TTICAD uses the family-aware decompose_tticad_safe() driver, not "
            "decomp_form_reduced_safe()"
        )
    raise ValueError(f"unknown reduced CAD backend: {backend!r}")


def _scan_reduced_conditions(cad: CompleteCAD, *, backend: ReducedBackend) -> ReducedSideReport:
    """Check reduced-CAD base samples for nullification and Lazard valuations."""
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
            if cad.tower.metadata.get("ec_lifting") and not ec_parent_admissible(cad.tower, parent):
                continue
            plans = {plan.level: plan for plan in cad.tower.metadata.get("ec_levels", ())}
            checked_polys = level_polys
            if level in plans and plans[level].operator == "ec-sign":
                checked_polys = (sp.Poly(plans[level].equational_constraint, *variables[:level]),)
            for poly in checked_polys:
                if poly.degree(variable) <= 0:
                    continue
                context = FiberRootContext(poly, variables[:level], parent.sample, variable)
                nullified = context.is_nullified_at_parent()
                valuation = ()
                if backend == "lazard":
                    from .lifting.lazard import native_lazard_evaluate

                    valuation = (0,) * (level - 1)
                    if nullified:
                        valuation = native_lazard_evaluate(
                            poly, variables[:level], parent.sample
                        ).valuation
                        nullified = False
                if backend == "lazard" and valuation:
                    valuations[(level, parent.index, _poly_key(poly))] = valuation
                if nullified:
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
                with contextlib.suppress(
                    sp.PolynomialError, TypeError, ValueError, NotImplementedError, ArithmeticError
                ):
                    repair.append(sp.Poly(coeff, *lower_gens))
        if as_univar.degree() > 1 and lower_gens:
            from .polynomial_utils import derivative_resultant_expression

            for candidate in (
                sp.discriminant(expr, var),
                derivative_resultant_expression(as_univar, var),
            ):
                candidate = sp.expand(candidate)
                if candidate != 0:
                    with contextlib.suppress(
                        sp.PolynomialError,
                        TypeError,
                        ValueError,
                        NotImplementedError,
                        ArithmeticError,
                    ):
                        repair.append(sp.Poly(candidate, *lower_gens))
    except (sp.PolynomialError, TypeError, ValueError, NotImplementedError, ArithmeticError):
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


def _intrinsic_formula_certificate(
    *,
    reduced: CompleteCAD,
    formula: Formula,
    side_conditions: ReducedSideReport,
    backend: ReducedBackend,
) -> ReducedCertificate | None:
    """Certify formula truth-invariance from the reduced-projection theorem.

    Recursive lifting uses the recorded leading coefficients, discriminants
    and resultants to preserve degree and order on admissible EC sections.
    Positive Lazard valuations do not establish truth invariance of the
    original EC: its Taylor residue may be nonzero where the EC itself vanishes.
    Missing hypotheses require an ambient sign-invariant fallback.
    """

    constraints = tuple(equational_constraints(formula))
    reduced_levels = tuple(reduced.tower.metadata.get("reduced_levels", ()))
    variables = reduced.tower.variables
    if reduced.tower.metadata.get("ec_lifting"):
        if (
            not constraints
            or not side_conditions.valid
            or not replay_ec_projection(reduced.tower, formula)
        ):
            return None
        if any(
            any(order > 0 for order in valuation)
            for valuation in side_conditions.lazard_valuations.values()
        ):
            return None
        for plan in reduced.tower.metadata["ec_levels"]:
            if plan.level <= 1:
                continue
            for parent in reduced.cells_by_level[plan.level - 1]:
                if not ec_parent_admissible(reduced.tower, parent):
                    continue
                for coefficient in plan.degree_coefficients:
                    if (
                        sign_at_sample(
                            sp.Poly(coefficient, *variables[: plan.level - 1], domain="EX"),
                            parent.sample,
                        )
                        == 0
                    ):
                        return None
        return ReducedCertificate(
            valid=True,
            invariant="truth",
            checked_refinement_cells=0,
            checked_reduced_cells=len(reduced.cells),
            grouped_values={
                (-1,): (
                    "recursive EC degree/order induction",
                    "replayed EC implications and projection closure",
                )
            },
        )
    if not constraints or reduced_levels != (len(variables),) or not side_conditions.valid:
        return None
    if backend == "lazard" and any(
        any(order > 0 for order in valuation)
        for valuation in side_conditions.lazard_valuations.values()
    ):
        return None
    active = [
        poly
        for poly in reduced.tower.level(len(variables)).polynomials
        if poly.degree(variables[-1]) > 0
    ]
    designated = _candidate_ec(active, constraints)
    if designated is None:
        return None
    leading = univariate_leading_coefficient(designated, variables[-1])
    if leading.free_symbols:
        leading_poly = sp.Poly(leading, *variables[:-1], domain="EX")
        # Lower levels use full projection, so the projected leading
        # coefficient has invariant sign on each of these base cells.
        if any(
            sign_at_sample(leading_poly, parent.sample) == 0
            for parent in reduced.cells_by_level.get(len(variables) - 1, ())
        ):
            return None
    checked = (
        f"{backend} reduced-projection theorem",
        "formula-wide equational constraint",
        "single top-level reduced step",
        "degree-invariant designated equational constraint",
        "nullification-free reduced lifting",
    )
    return ReducedCertificate(
        valid=True,
        invariant="truth",
        checked_refinement_cells=0,
        checked_reduced_cells=len(reduced.cells),
        grouped_values={(-1,): checked},
    )


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
    if formula is not None:
        intrinsic = _intrinsic_formula_certificate(
            reduced=reduced_cad, formula=formula, side_conditions=side_conditions, backend=backend
        )
        if intrinsic is not None:
            validity = ProjectionValidity(
                theory=backend,
                valid=True,
                complete_if_used=True,
                reason=f"active {backend} reduced CAD certified by reduced-projection theorem",
                checked_conditions=tuple(
                    dict.fromkeys(
                        (
                            *reduced_projection.validity.checked_conditions,
                            "formula-wide equational constraint",
                            *side_conditions.checked_conditions,
                            f"{backend} reduced-projection theorem",
                        )
                    )
                ),
                failed_conditions=(),
                fallback_backend=None,
                details={**dict(reduced_projection.validity.details), "certificate": intrinsic},
            )
            return (
                _with_invariant(reduced_cad, "truth", formulas=(to_sympy(formula),)),
                side_conditions,
                intrinsic,
                validity,
            )
    if formula is not None and reduced_cad.tower.metadata.get("ec_lifting"):
        # Sample agreement with a reference CAD cannot establish delineability
        # of an EC-only stack whose degree/order induction has failed.
        certificate = ReducedCertificate(
            False,
            "truth",
            0,
            len(reduced_cad.cells),
            ("recursive EC degree/order hypotheses could not be established",),
        )
        validity = proj_validity_from_cert(
            reduced_projection.validity, certificate, backend=backend
        )
        return reduced_cad, side_conditions, certificate, validity
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
    if certificate.valid and side_conditions.valid:
        reduced_cad = _with_invariant(
            reduced_cad,
            certificate.invariant,
            formulas=() if formula is None else (to_sympy(formula),),
        )
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

    polys = tuple(polys)
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
        effective_backend=collins_cad.backend,
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
    from .projection.reduced import build_reduced_proj_tower

    reduced_projection = build_reduced_proj_tower(
        polys,
        vars_tuple,
        theory=backend,
        equational_constraints=equational_constraints(formula),
        certify=certify_reduced,
        recursive_ec=True,
    )
    collins_cad = None
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

    if collins_cad is None:
        collins_cad = decomp_collins_complete(polys, vars_tuple)

    return SafeReducedCAD(
        requested_backend=backend,
        effective_backend=collins_cad.backend,
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
