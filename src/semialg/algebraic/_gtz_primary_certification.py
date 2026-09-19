"""Schema validation and independent replay for GTZ primary decomposition."""

from __future__ import annotations

import sympy as sp

from ._gtz_primary_types import (
    FractionFieldZeroDimensionalCertificate,
    GTZContractedComponentCertificate,
    GTZNodeCertificate,
    GTZPrimaryDecompositionCertificate,
)
from .equality_ideal import EqualityIdealContext
from .gtz import (
    IndependentLocalizationCertificate,
    SaturationStabilizationCertificate,
    verify_independent_localization_certificate,
    verify_localization_contraction_certificate,
    verify_saturation_stabilization_certificate,
)
from .gtz_primary import (
    _canonical_qq_basis,
    _cleanup_components,
    _component_from_zero_dim,
    _ff_ideal_equal,
    _ideal_intersection,
    _intersection_all,
    _qq_ideal_equal,
    _remove_redundant,
    verify_fraction_field_zero_dimensional_certificate,
)
from .gtz_zero_dim import (
    ZeroDimensionalPrimaryCertificate,
    verify_zero_dimensional_primary_certificate,
)
from .hilbert import ideal_degree


def _validate_gtz_node_schema(node: GTZNodeCertificate) -> bool:
    """Cheap structural checks before any Groebner or ideal arithmetic."""
    if not isinstance(node, GTZNodeCertificate):
        return False
    vars_ = tuple(node.variables)
    if len(set(vars_)) != len(vars_) or any(not isinstance(v, sp.Symbol) for v in vars_):
        return False
    if node.kind == "unit":
        return (
            node.zero_dimensional_certificate is None
            and node.localization_certificate is None
            and node.local_primary_certificate is None
            and not node.contracted_components
            and not node.hull_generators
            and node.split_certificate is None
            and node.residual is None
        )
    if node.kind == "zero_dimensional":
        cert = node.zero_dimensional_certificate
        return (
            isinstance(cert, ZeroDimensionalPrimaryCertificate)
            and tuple(cert.variables) == vars_
            and node.localization_certificate is None
            and node.local_primary_certificate is None
            and not node.contracted_components
            and not node.hull_generators
            and node.split_certificate is None
            and node.residual is None
        )
    if node.kind not in {"localized", "localized_split"}:
        return False
    localization = node.localization_certificate
    local = node.local_primary_certificate
    if not isinstance(localization, IndependentLocalizationCertificate):
        return False
    if not isinstance(local, FractionFieldZeroDimensionalCertificate):
        return False
    if tuple(localization.variables) != vars_:
        return False
    if tuple(local.parameters) != tuple(localization.independent_variables):
        return False
    if tuple(local.variables) != tuple(localization.dependent_variables):
        return False
    if len(node.contracted_components) != len(local.components):
        return False
    if any(
        not isinstance(item, GTZContractedComponentCertificate)
        for item in node.contracted_components
    ):
        return False
    if node.zero_dimensional_certificate is not None:
        return False
    if node.kind == "localized":
        return node.split_certificate is None and node.residual is None
    return (
        isinstance(node.split_certificate, SaturationStabilizationCertificate)
        and isinstance(node.residual, GTZNodeCertificate)
        and tuple(node.split_certificate.variables) == vars_
        and tuple(node.residual.variables) == vars_
    )


def _validate_gtz_certificate_schema(certificate: GTZPrimaryDecompositionCertificate) -> bool:
    if not isinstance(certificate, GTZPrimaryDecompositionCertificate):
        return False
    vars_ = tuple(certificate.variables)
    return (
        isinstance(certificate.root, GTZNodeCertificate)
        and tuple(certificate.root.variables) == vars_
        and _validate_gtz_node_schema(certificate.root)
        and all(
            hasattr(component, "equations") and hasattr(component, "radical")
            for component in certificate.final_components
        )
    )


def _verify_node(node: GTZNodeCertificate) -> tuple[bool, tuple[object, ...]]:
    from ..algebraic_decomposition import PrimaryComponent

    if not _validate_gtz_node_schema(node):
        return False, tuple()
    try:
        vars_ = tuple(node.variables)
        source = _canonical_qq_basis(node.source_generators, vars_)
        if source != tuple(node.source_generators):
            return False, tuple()
        context = EqualityIdealContext.build(source, vars_)
        if int(context.dimension) != int(node.dimension):
            return False, tuple()
        if node.kind == "unit":
            return context.inconsistent, tuple()
        if node.kind == "zero_dimensional":
            certificate = node.zero_dimensional_certificate
            if certificate is None or not verify_zero_dimensional_primary_certificate(certificate):
                return False, tuple()
            if not _qq_ideal_equal(certificate.source_generators, source, vars_):
                return False, tuple()
            components = tuple(
                _component_from_zero_dim(component, vars_) for component in certificate.components
            )
            return True, components
        if node.kind not in {"localized", "localized_split"}:
            return False, tuple()
        localization = node.localization_certificate
        local = node.local_primary_certificate
        if localization is None or local is None:
            return False, tuple()
        if not verify_independent_localization_certificate(localization):
            return False, tuple()
        if not _qq_ideal_equal(localization.source_generators, source, vars_):
            return False, tuple()
        if tuple(local.parameters) != tuple(localization.independent_variables):
            return False, tuple()
        if tuple(local.variables) != tuple(localization.dependent_variables):
            return False, tuple()
        if not _ff_ideal_equal(
            local.source_generators,
            localization.localized_basis,
            local.variables,
            local.parameters,
        ):
            return False, tuple()
        if not verify_fraction_field_zero_dimensional_certificate(local):
            return False, tuple()
        if len(node.contracted_components) != len(local.components):
            return False, tuple()
        top = []
        for evidence, local_component in zip(
            node.contracted_components, local.components, strict=True
        ):
            if not verify_localization_contraction_certificate(evidence.primary_contraction):
                return False, tuple()
            if not verify_localization_contraction_certificate(evidence.radical_contraction):
                return False, tuple()
            if not _ff_ideal_equal(
                evidence.primary_contraction.localized_generators,
                local_component.generators,
                local.variables,
                local.parameters,
            ):
                return False, tuple()
            if not _ff_ideal_equal(
                evidence.radical_contraction.localized_generators,
                local_component.radical,
                local.variables,
                local.parameters,
            ):
                return False, tuple()
            equations = _canonical_qq_basis(
                evidence.primary_contraction.contracted_generators, vars_
            )
            radical = _canonical_qq_basis(evidence.radical_contraction.contracted_generators, vars_)
            if equations != tuple(evidence.equations) or radical != tuple(evidence.radical):
                return False, tuple()
            qctx = EqualityIdealContext.build(equations, vars_)
            pctx = EqualityIdealContext.build(radical, vars_)
            if qctx.dimension != node.dimension or pctx.dimension != node.dimension:
                return False, tuple()
            if evidence.dimension != node.dimension:
                return False, tuple()
            if ideal_degree(equations, vars_) != evidence.degree:
                return False, tuple()
            top.append(
                PrimaryComponent(
                    equations,
                    radical,
                    evidence.dimension,
                    evidence.degree,
                    True,
                    "gtz_localized_primary_contraction",
                )
            )
        hull = _intersection_all([component.equations for component in top], vars_)
        if not _qq_ideal_equal(hull, node.hull_generators, vars_):
            return False, tuple()
        if node.kind == "localized":
            if node.split_certificate is not None or node.residual is not None:
                return False, tuple()
            if not _qq_ideal_equal(hull, source, vars_):
                return False, tuple()
            return True, tuple(top)
        split = node.split_certificate
        residual = node.residual
        if split is None or residual is None:
            return False, tuple()
        if not verify_saturation_stabilization_certificate(split):
            return False, tuple()
        if not _qq_ideal_equal(split.source_generators, source, vars_):
            return False, tuple()
        if not _qq_ideal_equal(split.saturated_generators, hull, vars_):
            return False, tuple()
        if tuple(residual.variables) != vars_:
            return False, tuple()
        if not _qq_ideal_equal(residual.source_generators, split.companion_generators, vars_):
            return False, tuple()
        if residual.dimension > node.dimension and residual.kind != "unit":
            return False, tuple()
        if _qq_ideal_equal(residual.source_generators, source, vars_):
            return False, tuple()
        ok, residual_components = _verify_node(residual)
        if not ok:
            return False, tuple()
        reconstructed = _ideal_intersection(hull, split.companion_generators, vars_)
        if not _qq_ideal_equal(reconstructed, source, vars_):
            return False, tuple()
        return True, tuple(top) + tuple(residual_components)
    except (
        ArithmeticError,
        ValueError,
        TypeError,
        ZeroDivisionError,
        sp.PolynomialError,
        sp.polys.polyerrors.CoercionFailed,
    ):
        return False, tuple()


def verify_gtz_primary_decomposition_certificate(
    certificate: GTZPrimaryDecompositionCertificate,
) -> bool:
    """Replay the complete recursive GTZ proof tree independently."""
    if not _validate_gtz_certificate_schema(certificate):
        return False
    try:
        vars_ = tuple(certificate.variables)
        source = _canonical_qq_basis(certificate.source_generators, vars_)
        if source != tuple(certificate.source_generators):
            return False
        if tuple(certificate.root.variables) != vars_:
            return False
        if not _qq_ideal_equal(certificate.root.source_generators, source, vars_):
            return False
        ok, raw = _verify_node(certificate.root)
        if not ok:
            return False
        cleaned = _cleanup_components(raw, vars_)
        if len(cleaned) != len(certificate.final_components):
            return False
        for left, right in zip(cleaned, certificate.final_components, strict=True):
            if not _qq_ideal_equal(left.equations, right.equations, vars_):
                return False
            if not _qq_ideal_equal(left.radical, right.radical, vars_):
                return False
            if left.dimension != right.dimension or left.degree != right.degree:
                return False
            if not right.primary:
                return False
        reconstructed = _intersection_all(
            [component.equations for component in certificate.final_components], vars_
        )
        if not _qq_ideal_equal(reconstructed, source, vars_):
            return False
        radicals = [
            tuple(_canonical_qq_basis(component.radical, vars_))
            for component in certificate.final_components
        ]
        irredundant = len(set(radicals)) == len(radicals)
        if irredundant != certificate.irredundant:
            return False
        # Exact ideal redundancy check as a final independent cleanup proof.
        if tuple(_remove_redundant(certificate.final_components, vars_)) != tuple(
            certificate.final_components
        ):
            return False
        return True
    except (
        ArithmeticError,
        ValueError,
        TypeError,
        ZeroDivisionError,
        sp.PolynomialError,
        sp.polys.polyerrors.CoercionFailed,
    ):
        return False
