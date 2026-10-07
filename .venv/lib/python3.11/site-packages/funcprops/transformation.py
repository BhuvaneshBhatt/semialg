"""Transformation geometry and change-of-variables certificates."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import sympy as sp
from sympy import S
from sympy.assumptions import Q, ask

from ._semialg import SemialgUnavailable
from .exactness import certify_identity_equal, certify_zero
from .exceptions import SYMBOLIC_QUERY_ERRORS
from .real_analysis.common import (
    normalize_domain_condition,
    normalize_expressions,
    normalize_variables,
)
from .registry_derivatives import JacobianCertificate, registry_jacobian
from .results import Sign, TruthValue


def _domain_condition(domain, variables):
    try:
        return normalize_domain_condition(domain, variables)
    except (TypeError, ValueError):
        return None


def _truth_from_nonzero(expr, assumptions=S.true):
    expr = sp.sympify(expr)
    zero = certify_zero(expr, assumptions=assumptions)
    if zero is True:
        source = "identically-zero" if assumptions is S.true else "assumption-zero"
        return TruthValue.FALSE, source
    if zero is False:
        source = "constant-nonzero" if not expr.free_symbols else "assumption-nonzero"
        return TruthValue.TRUE, source
    return TruthValue.UNKNOWN, "unresolved"


@dataclass(frozen=True)
class DeterminantNonvanishingCertificate:
    determinant: Any
    variables: tuple[sp.Symbol, ...]
    domain: Any
    classification: TruthValue
    source: str
    witness: Any = None
    assumptions: Any = S.true

    @property
    def primary(self):
        return self.classification

    def verify(self):
        det = sp.sympify(self.determinant)
        if self.classification is TruthValue.TRUE:
            if self.source == "constant-nonzero":
                return not det.free_symbols and det != 0
            if self.source == "assumption-nonzero":
                try:
                    return ask(Q.nonzero(det), self.assumptions) is True
                except SYMBOLIC_QUERY_ERRORS:
                    return False
            if self.source == "semialgebraic-no-zero":
                condition = _domain_condition(self.domain, self.variables)
                try:
                    from ._semialg import is_satisfiable

                    query = (
                        sp.Eq(det, 0)
                        if condition is S.true
                        else sp.And(condition, sp.Eq(det, 0))
                    )
                    return not is_satisfiable(query, self.variables)
                except (SemialgUnavailable, *SYMBOLIC_QUERY_ERRORS):
                    return False
            return False
        if self.classification is TruthValue.FALSE:
            if self.source in {"identically-zero", "assumption-zero"}:
                return certify_zero(det, assumptions=self.assumptions) is True
            if self.source == "semialgebraic-zero-witness":
                condition = _domain_condition(self.domain, self.variables)
                if condition is None:
                    return False
                try:
                    from ._semialg import is_satisfiable

                    query = (
                        sp.Eq(det, 0)
                        if condition is S.true
                        else sp.And(condition, sp.Eq(det, 0))
                    )
                    return is_satisfiable(query, self.variables) is True
                except (SemialgUnavailable, *SYMBOLIC_QUERY_ERRORS):
                    return False
            return False
        return self.source in {"unresolved", "nonsquare-map", "domain-unresolved"}


@dataclass(frozen=True)
class LocalInvertibilityCertificate:
    expressions: tuple[Any, ...]
    variables: tuple[sp.Symbol, ...]
    jacobian: JacobianCertificate
    determinant: Any
    nonvanishing: DeterminantNonvanishingCertificate
    classification: TruthValue
    source: str = "inverse-function-theorem"

    @property
    def primary(self):
        return self.classification

    def verify(self):
        if not self.jacobian.verify():
            return False
        if len(self.expressions) != len(self.variables):
            return (
                self.classification is TruthValue.UNKNOWN
                and self.source == "nonsquare-map"
            )
        if (
            certify_identity_equal(self.jacobian.jacobian.det(), self.determinant)
            is not True
        ):
            return False
        if self.classification is TruthValue.TRUE:
            return (
                self.nonvanishing.classification is TruthValue.TRUE
                and self.nonvanishing.verify()
            )
        if self.classification is TruthValue.FALSE:
            return (
                self.nonvanishing.classification is TruthValue.FALSE
                and self.nonvanishing.verify()
            )
        return self.nonvanishing.classification is TruthValue.UNKNOWN


@dataclass(frozen=True)
class ChangeOfVariablesCertificate:
    expressions: tuple[Any, ...]
    variables: tuple[sp.Symbol, ...]
    jacobian: JacobianCertificate
    determinant: Any
    absolute_determinant: Any
    local_invertibility: LocalInvertibilityCertificate
    classification: TruthValue
    orientation: Sign
    domain: Any = S.Reals
    source: str = "local-change-of-variables"

    @property
    def primary(self):
        return self.classification

    @property
    def volume_factor(self):
        return self.absolute_determinant

    def verify(self):
        if not self.jacobian.verify() or not self.local_invertibility.verify():
            return False
        if (
            certify_identity_equal(self.absolute_determinant, sp.Abs(self.determinant))
            is not True
        ):
            return False
        if self.classification is TruthValue.TRUE:
            return self.local_invertibility.classification is TruthValue.TRUE
        if self.classification is TruthValue.FALSE:
            return self.local_invertibility.classification is TruthValue.FALSE
        if self.classification is TruthValue.UNKNOWN:
            return self.local_invertibility.classification is TruthValue.UNKNOWN
        return False


@dataclass(frozen=True)
class TransformationGeometry:
    expressions: tuple[Any, ...]
    variables: tuple[sp.Symbol, ...]
    domain: Any
    jacobian: JacobianCertificate
    determinant: Any
    determinant_nonvanishing: DeterminantNonvanishingCertificate
    local_invertibility: LocalInvertibilityCertificate
    change_of_variables: ChangeOfVariablesCertificate


def jacobian(expressions, variables, *, return_result=False):
    """Return the Jacobian matrix; request its certificate with ``return_result=True``."""
    certificate = registry_jacobian(
        normalize_expressions(expressions), normalize_variables(variables)
    )
    return certificate if return_result else certificate.jacobian


def _det_nonzero_from_jacobian(
    expressions, variables, jacobian, *, domain=S.Reals, assumptions=S.true
):
    if len(expressions) != len(variables):
        return DeterminantNonvanishingCertificate(
            S.NaN,
            variables,
            domain,
            TruthValue.UNKNOWN,
            "nonsquare-map",
            assumptions=assumptions,
        )
    det = sp.factor(sp.trigsimp(jacobian.jacobian.det()))
    verdict, source = _truth_from_nonzero(det, assumptions)
    if verdict is not TruthValue.UNKNOWN:
        return DeterminantNonvanishingCertificate(
            det, variables, domain, verdict, source, assumptions=assumptions
        )

    condition = _domain_condition(domain, variables)
    if condition is None:
        return DeterminantNonvanishingCertificate(
            det,
            variables,
            domain,
            TruthValue.UNKNOWN,
            "domain-unresolved",
            assumptions=assumptions,
        )

    # A parameter is not a transformation coordinate.  Semialgebraic zero
    # search must not choose parameter values as if they were points in the
    # requested domain.
    unresolved_parameters = det.free_symbols.difference(variables)
    if unresolved_parameters:
        return DeterminantNonvanishingCertificate(
            det,
            variables,
            domain,
            TruthValue.UNKNOWN,
            "parameter-dependent-determinant",
            assumptions=assumptions,
        )
    try:
        from ._semialg import is_satisfiable

        query = (
            sp.Eq(det, 0) if condition is S.true else sp.And(condition, sp.Eq(det, 0))
        )
        sat = is_satisfiable(query, variables)
        if sat is False:
            return DeterminantNonvanishingCertificate(
                det,
                variables,
                domain,
                TruthValue.TRUE,
                "semialgebraic-no-zero",
                assumptions=assumptions,
            )
        if sat is True:
            return DeterminantNonvanishingCertificate(
                det,
                variables,
                domain,
                TruthValue.FALSE,
                "semialgebraic-zero-witness",
                witness=query,
                assumptions=assumptions,
            )
    except (SemialgUnavailable, *SYMBOLIC_QUERY_ERRORS):
        pass
    return DeterminantNonvanishingCertificate(
        det,
        variables,
        domain,
        TruthValue.UNKNOWN,
        "unresolved",
        assumptions=assumptions,
    )


def determinant_nonvanishing(
    expressions, variables, *, domain=S.Reals, assumptions=S.true, return_result=False
):
    """Return whether the Jacobian determinant is nonzero throughout the domain."""
    expressions = normalize_expressions(expressions)
    variables = normalize_variables(variables)
    jacobian = registry_jacobian(expressions, variables)
    certificate = _det_nonzero_from_jacobian(
        expressions, variables, jacobian, domain=domain, assumptions=assumptions
    )
    return certificate if return_result else certificate.classification


def _local_invertibility_from_facts(expressions, variables, jacobian, nonvanishing):
    if len(expressions) != len(variables):
        return LocalInvertibilityCertificate(
            expressions,
            variables,
            jacobian,
            S.NaN,
            nonvanishing,
            TruthValue.UNKNOWN,
            "nonsquare-map",
        )
    determinant = sp.factor(sp.trigsimp(jacobian.jacobian.det()))
    return LocalInvertibilityCertificate(
        expressions,
        variables,
        jacobian,
        determinant,
        nonvanishing,
        nonvanishing.classification,
    )


def local_invertibility(
    expressions, variables, *, domain=S.Reals, assumptions=S.true, return_result=False
):
    """Return whether the inverse function theorem certifies local invertibility."""
    expressions = normalize_expressions(expressions)
    variables = normalize_variables(variables)
    jacobian = registry_jacobian(expressions, variables)
    nonvanishing = _det_nonzero_from_jacobian(
        expressions, variables, jacobian, domain=domain, assumptions=assumptions
    )
    certificate = _local_invertibility_from_facts(
        expressions, variables, jacobian, nonvanishing
    )
    return certificate if return_result else certificate.classification


def _orientation(det, assumptions=S.true):
    try:
        if ask(Q.positive(det), assumptions) is True:
            return Sign.POSITIVE
        if ask(Q.negative(det), assumptions) is True:
            return Sign.NEGATIVE
    except SYMBOLIC_QUERY_ERRORS:
        pass
    return Sign.UNKNOWN


def _change_of_variables_from_facts(
    expressions, variables, jacobian, local_invertibility, *, domain, assumptions
):
    determinant = local_invertibility.determinant
    absolute_determinant = sp.Abs(determinant) if determinant is not S.NaN else S.NaN
    return ChangeOfVariablesCertificate(
        expressions,
        variables,
        jacobian,
        determinant,
        absolute_determinant,
        local_invertibility,
        local_invertibility.classification,
        _orientation(determinant, assumptions)
        if determinant is not S.NaN
        else Sign.UNKNOWN,
        domain,
    )


def change_of_variables(expressions, variables, *, domain=S.Reals, assumptions=S.true):
    """Certify the local differential change-of-variables factor ``Abs(det J)``."""
    expressions = normalize_expressions(expressions)
    variables = normalize_variables(variables)
    jacobian = registry_jacobian(expressions, variables)
    nonvanishing = _det_nonzero_from_jacobian(
        expressions, variables, jacobian, domain=domain, assumptions=assumptions
    )
    local = _local_invertibility_from_facts(
        expressions, variables, jacobian, nonvanishing
    )
    return _change_of_variables_from_facts(
        expressions,
        variables,
        jacobian,
        local,
        domain=domain,
        assumptions=assumptions,
    )


def transformation_geometry(
    expressions, variables, *, domain=S.Reals, assumptions=S.true
):
    """Return the combined transformation-geometry summary."""
    expressions = normalize_expressions(expressions)
    variables = normalize_variables(variables)
    jacobian = registry_jacobian(expressions, variables)
    nonvanishing = _det_nonzero_from_jacobian(
        expressions, variables, jacobian, domain=domain, assumptions=assumptions
    )
    local = _local_invertibility_from_facts(
        expressions, variables, jacobian, nonvanishing
    )
    change_of_variables = _change_of_variables_from_facts(
        expressions,
        variables,
        jacobian,
        local,
        domain=domain,
        assumptions=assumptions,
    )
    return TransformationGeometry(
        expressions,
        variables,
        domain,
        jacobian,
        local.determinant,
        nonvanishing,
        local,
        change_of_variables,
    )
