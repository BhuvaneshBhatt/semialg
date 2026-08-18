from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass

import sympy as sp

from ..cad.decomposition import CompleteCAD, decomp_collins_complete, decomp_from_proj_tower
from ..cad.projection.reduced import (
    ProjectionValidity,
    ReducedProjectionTower,
    build_reduced_proj_tower,
)
from ..cad.reduced import (
    ReducedCertificate,
    ReducedSideReport,
    _augment_tower_repairs,
    _locate_reduced_cell,
    _scan_reduced_conditions,
    _validity_with_condition,
    proj_validity_from_cert,
)
from ..formula import Formula, formula_polynomials
from ..qe.complete import evaluate_formula_on_cell
from ..tti import FormulaFamily, extract_formula_families


@dataclass(frozen=True)
class TTICADValidity:
    """Truth-table invariance status for a TTICAD attempt."""

    valid: bool
    reason: str
    family_count: int
    designated_constraints: tuple[sp.Expr, ...]
    checked_conditions: tuple[str, ...]
    failed_conditions: tuple[str, ...]
    family_truth_vectors: Mapping[tuple[int, ...], tuple[tuple[bool, ...], ...]] | None = None


@dataclass(frozen=True)
class SafeTTICAD:
    """TTICAD projection attempt with proof-carrying Collins fallback."""

    formula: Formula
    variables: tuple[sp.Symbol, ...]
    families: tuple[FormulaFamily, ...]
    reduced_projection: ReducedProjectionTower
    validity: TTICADValidity
    projection_validity: ProjectionValidity
    cad: CompleteCAD
    effective_backend: str = "collins-complete"
    complete: bool = True
    used_fallback: bool = True
    certificate: ReducedCertificate | None = None
    side_conditions: ReducedSideReport | None = None
    fallback_cad: CompleteCAD | None = None

    @property
    def family_count(self) -> int:
        return len(self.families)


def _designated_constraints(families: Sequence[FormulaFamily]) -> tuple[sp.Expr, ...]:
    out: list[sp.Expr] = []
    for family in families:
        if family.designated_ec is not None:
            expr = sp.expand(family.designated_ec)
            if expr not in out:
                out.append(expr)
    return tuple(out)


def certify_truth_table_inv(
    *,
    reduced: CompleteCAD,
    refiner: CompleteCAD,
    families: Sequence[FormulaFamily],
    variables: Sequence[sp.Symbol],
) -> ReducedCertificate:
    grouped: dict[tuple[int, ...], set[tuple[bool, ...]]] = {}
    failures: list[str] = []
    for full_cell in refiner.cells:
        reduced_cell = _locate_reduced_cell(reduced, full_cell)
        if reduced_cell is None:
            failures.append(f"could not locate reduced cell for full cell {full_cell.index}")
            continue
        vector = tuple(
            evaluate_formula_on_cell(family.formula, full_cell, variables) for family in families
        )
        grouped.setdefault(reduced_cell.index, set()).add(vector)
    for index, values in grouped.items():
        if len(values) > 1:
            failures.append(
                f"truth table is not invariant on reduced cell {index}: {sorted(values)!r}"
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


def decompose_tticad_safe(
    formula: Formula,
    variables: Sequence[sp.Symbol],
    *,
    certify_reduced: bool = True,
) -> SafeTTICAD:
    """Attempt active TTICAD and fall back unless truth-table invariance is certified.

    The reduced TTICAD tower is actively lifted. Completeness is accepted only
    if a full Collins refinement proves that every reduced cell has an invariant
    truth vector for all formula families. Otherwise the returned CAD is the
    Collins-complete fallback and the failure reasons are recorded.
    """

    vars_tuple = tuple(variables)
    families = extract_formula_families(formula)
    ecs = _designated_constraints(families)
    polys = formula_polynomials(formula)
    reduced = build_reduced_proj_tower(
        polys,
        vars_tuple,
        theory="tticad",
        equational_constraints=ecs,
        certify=certify_reduced,
    )
    reduced_cad = decomp_from_proj_tower(reduced.tower, backend="tticad-reduced")
    side_conditions = _scan_reduced_conditions(reduced_cad, backend="tticad")
    collins_cad = decomp_collins_complete(polys, vars_tuple)
    certificate = certify_truth_table_inv(
        reduced=reduced_cad,
        refiner=collins_cad,
        families=families,
        variables=vars_tuple,
    )
    cert_projection_validity = proj_validity_from_cert(
        reduced.validity, certificate, backend="tticad"
    )
    projection_validity = _validity_with_condition(
        cert_projection_validity, side_conditions, backend="tticad"
    )
    if certificate.valid and side_conditions.valid:
        validity = TTICADValidity(
            valid=True,
            reason="active TTICAD certified by Collins-refinement truth-table-invariance check",
            family_count=len(families),
            designated_constraints=ecs,
            checked_conditions=(
                "family extraction",
                "projection construction",
                *side_conditions.checked_conditions,
                "reduced lifting",
                "truth-table invariance",
            ),
            failed_conditions=(),
            family_truth_vectors=certificate.grouped_values,
        )
        return SafeTTICAD(
            formula=formula,
            variables=vars_tuple,
            families=families,
            reduced_projection=reduced,
            validity=validity,
            projection_validity=projection_validity,
            cad=reduced_cad,
            effective_backend="tticad-reduced-certified",
            used_fallback=False,
            certificate=certificate,
            side_conditions=side_conditions,
            fallback_cad=collins_cad,
        )

    repaired = _augment_tower_repairs(reduced, side_conditions, vars_tuple)
    if repaired is not None:
        repaired_cad = decomp_from_proj_tower(repaired.tower, backend="tticad-repaired")
        repaired_side = _scan_reduced_conditions(repaired_cad, backend="tticad")
        repaired_cert = certify_truth_table_inv(
            reduced=repaired_cad,
            refiner=collins_cad,
            families=families,
            variables=vars_tuple,
        )
        repaired_validity = proj_validity_from_cert(
            repaired.validity, repaired_cert, backend="tticad"
        )
        repaired_validity = _validity_with_condition(
            repaired_validity, repaired_side, backend="tticad"
        )
        if repaired_cert.valid and repaired_side.valid:
            validity = TTICADValidity(
                valid=True,
                reason="TTICAD accepted after localized delineating-polynomial repair and truth-table certification",
                family_count=len(families),
                designated_constraints=ecs,
                checked_conditions=(
                    "family extraction",
                    "projection construction",
                    *repaired_side.checked_conditions,
                    "localized repair",
                    "truth-table invariance",
                ),
                failed_conditions=(),
                family_truth_vectors=repaired_cert.grouped_values,
            )
            return SafeTTICAD(
                formula=formula,
                variables=vars_tuple,
                families=families,
                reduced_projection=repaired,
                validity=validity,
                projection_validity=repaired_validity,
                cad=repaired_cad,
                effective_backend="tticad-repaired-certified",
                used_fallback=False,
                certificate=repaired_cert,
                side_conditions=repaired_side,
                fallback_cad=collins_cad,
            )
    validity = TTICADValidity(
        valid=False,
        reason="active TTICAD failed truth-table-invariance certification; using Collins fallback",
        family_count=len(families),
        designated_constraints=ecs,
        checked_conditions=(
            "family extraction",
            "projection construction",
            "reduced lifting",
            "truth-table invariance",
        ),
        failed_conditions=tuple(
            dict.fromkeys((*side_conditions.failed_conditions, *certificate.failures))
        ),
        family_truth_vectors=certificate.grouped_values,
    )
    return SafeTTICAD(
        formula=formula,
        variables=vars_tuple,
        families=families,
        reduced_projection=reduced,
        validity=validity,
        projection_validity=projection_validity,
        cad=collins_cad,
        effective_backend="collins-complete",
        used_fallback=True,
        certificate=certificate,
        side_conditions=side_conditions,
        fallback_cad=collins_cad,
    )


__all__ = ["SafeTTICAD", "TTICADValidity", "decompose_tticad_safe"]
