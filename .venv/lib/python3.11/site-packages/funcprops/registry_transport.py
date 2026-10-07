"""Certified property transport for native registry theorems.

The transport layer maps localized theorem subjects into user analysis variables.
It supports scalar paths as well as genuinely multivariate subject maps while
keeping theorem instantiation separate from property-specific composition logic.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import sympy as sp
from sympy import S

from .conditional import restrict_domain
from .registry import FunctionSlot, InstantiatedPropertyRule
from .results import (
    Analyticity,
    Continuity,
    Convexity,
    Differentiability,
    Monotonicity,
    Sign,
    TruthValue,
)


@dataclass(frozen=True)
class AnalysisSubject:
    """Variables and domain on which a transported property is queried."""

    variables: tuple[sp.Symbol, ...]
    domain: Any = S.Reals

    def __post_init__(self):
        variables = tuple(sp.sympify(v) for v in self.variables)
        if not variables or not all(isinstance(v, sp.Symbol) for v in variables):
            raise TypeError(
                "AnalysisSubject variables must be a nonempty tuple of Symbols"
            )
        if len(set(variables)) != len(variables):
            raise ValueError("AnalysisSubject variables cannot contain duplicates")
        object.__setattr__(self, "variables", variables)


@dataclass(frozen=True)
class SubjectMap:
    """Map registry subject slots into one or more analysis variables.

    The full Jacobian is stored for multivariable analysis, while ``tangent``
    provides the one-variable view used for scalar paths.
    """

    source_slots: tuple[FunctionSlot, ...]
    expressions: tuple[Any, ...]
    analysis_variables: tuple[sp.Symbol, ...]
    jacobian: sp.ImmutableMatrix

    @classmethod
    def from_instantiated(cls, instantiated: InstantiatedPropertyRule, variables):
        if isinstance(variables, (tuple, list, sp.Tuple)):
            vv = tuple(sp.sympify(v) for v in variables)
        else:
            vv = (sp.sympify(variables),)
        if not vv or not all(isinstance(v, sp.Symbol) for v in vv):
            raise TypeError("analysis variables must be Symbols")
        if len(set(vv)) != len(vv):
            raise ValueError("analysis variables cannot contain duplicates")
        expressions = tuple(sp.sympify(value) for value in instantiated.subject)
        jacobian = sp.ImmutableMatrix(
            [[sp.diff(value, var) for var in vv] for value in expressions]
        )
        return cls(tuple(instantiated.rule.subject), expressions, vv, jacobian)

    @property
    def is_scalar_path(self) -> bool:
        return len(self.analysis_variables) == 1

    @property
    def tangent(self) -> tuple[Any, ...]:
        if not self.is_scalar_path:
            raise ValueError("tangent is defined only for a scalar analysis path")
        return tuple(self.jacobian[:, 0])

    @property
    def varying_positions(self) -> tuple[int, ...]:
        return tuple(
            slot.position
            for row, slot in enumerate(self.source_slots)
            if any(self.jacobian[row, col] != 0 for col in range(self.jacobian.cols))
        )


@dataclass(frozen=True)
class TransportedPropertyRule:
    """An instantiated registry theorem transported to analysis variables."""

    instantiated: InstantiatedPropertyRule
    analysis_variables: tuple[sp.Symbol, ...]
    subject_map: SubjectMap
    conditions: Any
    result: Any
    complete: bool
    source: str = "native-registry-path-transport"
    evidence: tuple[Any, ...] = ()


def _compose_monotonicity(outer, inner):
    inc = {
        Monotonicity.STRICTLY_INCREASING,
        Monotonicity.NONDECREASING,
    }
    dec = {
        Monotonicity.STRICTLY_DECREASING,
        Monotonicity.NONINCREASING,
    }
    if outer is Monotonicity.CONSTANT or inner is Monotonicity.CONSTANT:
        return Monotonicity.CONSTANT
    if outer in inc and inner in inc:
        if (
            outer is Monotonicity.STRICTLY_INCREASING
            and inner is Monotonicity.STRICTLY_INCREASING
        ):
            return Monotonicity.STRICTLY_INCREASING
        return Monotonicity.NONDECREASING
    if outer in dec and inner in dec:
        if (
            outer is Monotonicity.STRICTLY_DECREASING
            and inner is Monotonicity.STRICTLY_DECREASING
        ):
            return Monotonicity.STRICTLY_INCREASING
        return Monotonicity.NONDECREASING
    if outer in inc and inner in dec:
        if (
            outer is Monotonicity.STRICTLY_INCREASING
            and inner is Monotonicity.STRICTLY_DECREASING
        ):
            return Monotonicity.STRICTLY_DECREASING
        return Monotonicity.NONINCREASING
    if outer in dec and inner in inc:
        if (
            outer is Monotonicity.STRICTLY_DECREASING
            and inner is Monotonicity.STRICTLY_INCREASING
        ):
            return Monotonicity.STRICTLY_DECREASING
        return Monotonicity.NONINCREASING
    return Monotonicity.UNKNOWN


def _monotonicity_from_derivative_sign(value):
    return {
        Sign.POSITIVE: Monotonicity.STRICTLY_INCREASING,
        Sign.NONNEGATIVE: Monotonicity.NONDECREASING,
        Sign.ZERO: Monotonicity.CONSTANT,
        Sign.NONPOSITIVE: Monotonicity.NONINCREASING,
        Sign.NEGATIVE: Monotonicity.STRICTLY_DECREASING,
    }.get(value, Monotonicity.UNKNOWN)


def _preimage_of_zero_result(result, subject_map: SubjectMap, var, domain):
    """Pull an outer zero description back along a scalar subject path."""
    expressions = subject_map.expressions
    dom = domain if isinstance(domain, sp.Set) else S.Reals

    if len(expressions) == 1:
        inner = expressions[0]
        if isinstance(result, sp.FiniteSet):
            pieces = []
            for root in result:
                try:
                    pieces.append(sp.solveset(sp.Eq(inner, root), var, domain=dom))
                except (NotImplementedError, TypeError, ValueError):
                    return None
            return sp.Union(*pieces) if pieces else S.EmptySet
        if result is S.EmptySet:
            return S.EmptySet
        if isinstance(result, sp.Equality):
            # Registry equalities are stated in the theorem subject.  Pull them
            # back through the instantiated subject map before solving in the
            # analysis variable; solving the outer equality directly can ignore
            # the inner map entirely.
            slots = subject_map.source_slots
            if len(slots) != 1:
                return None
            slot_symbol = slots[0]
            pulled = result.xreplace({slot_symbol: inner})
            try:
                return sp.solveset(pulled, var, domain=dom)
            except (NotImplementedError, TypeError, ValueError):
                return None
        return None

    if result is S.EmptySet:
        return S.EmptySet
    if not isinstance(result, sp.FiniteSet):
        return None

    pieces = []
    for point in result:
        coords = tuple(point) if isinstance(point, (tuple, sp.Tuple)) else None
        if coords is None or len(coords) != len(expressions):
            return None
        candidate = dom
        for inner, target in zip(expressions, coords, strict=True):
            try:
                solution = sp.solveset(sp.Eq(inner, target), var, domain=dom)
            except (NotImplementedError, TypeError, ValueError):
                return None
            candidate = sp.Intersection(candidate, solution)
        pieces.append(candidate)
    return sp.Union(*pieces) if pieces else S.EmptySet


def _inner_regularity(prop, expressions, var, restricted):
    evidence = []
    for inner in expressions:
        if prop == "continuity":
            from .singularities import continuous

            value = continuous(inner, var, domain=restricted)
            ok = value is Continuity.CONTINUOUS
        elif prop == "differentiability":
            from .real_analysis.differentiability import differentiability

            value = differentiability(inner, var, domain=restricted)
            ok = value in {
                Differentiability.DIFFERENTIABLE,
                Differentiability.CK,
                Differentiability.SMOOTH,
            }
        elif prop == "smoothness":
            from .real_analysis.smoothness import smoothness

            value = smoothness(inner, var, domain=restricted)
            ok = value is TruthValue.TRUE
        else:
            from .singularities import analytic

            value = analytic(inner, var, domain=restricted)
            ok = value is Analyticity.ANALYTIC
        evidence.append((inner, value))
        if not ok:
            return False, tuple(evidence)
    return True, tuple(evidence)


def _transport_scalar_rule(
    instantiated: InstantiatedPropertyRule,
    analysis_var,
    *,
    domain=True,
    registry=None,
) -> TransportedPropertyRule:
    """Transport a native property theorem along one scalar analysis path.

    One or more subject slots may be followed along a single analysis variable.
    Pointwise facts substitute directly; regularity requires every subject
    expression to have the required regularity; finite joint zero sets are
    pulled back by simultaneous preimage; and multiargument monotonicity may be
    certified from registered partial derivatives.
    """
    var = sp.sympify(analysis_var)
    if not isinstance(var, sp.Symbol):
        raise TypeError("analysis_var must be a Symbol")
    subject_map = SubjectMap.from_instantiated(instantiated, var)
    prop = instantiated.rule.property.lower()
    cond = instantiated.conditions
    expressions = subject_map.expressions

    if all(inner == var for inner in expressions) and len(expressions) == 1:
        return TransportedPropertyRule(
            instantiated,
            (var,),
            subject_map,
            cond,
            instantiated.result,
            True,
            evidence=("identity-subject-map",),
        )

    if prop == "sign":
        return TransportedPropertyRule(
            instantiated,
            (var,),
            subject_map,
            cond,
            instantiated.result,
            True,
            evidence=(("pointwise-path-composition", expressions),),
        )

    if prop in {"continuity", "differentiability", "smoothness", "analyticity"}:
        base_domain = S.Reals if domain is True else domain
        restricted = restrict_domain(base_domain, cond, var)
        ok, inner_evidence = _inner_regularity(prop, expressions, var, restricted)
        return TransportedPropertyRule(
            instantiated,
            (var,),
            subject_map,
            cond,
            instantiated.result if ok else None,
            ok,
            evidence=(("inner-regularity", inner_evidence),),
        )

    if prop in {"zeros", "real_zeros", "zero_set"}:
        pulled = _preimage_of_zero_result(instantiated.result, subject_map, var, domain)
        return TransportedPropertyRule(
            instantiated,
            (var,),
            subject_map,
            cond,
            pulled,
            pulled is not None,
            evidence=(("path-preimage", instantiated.result),),
        )

    if prop == "monotonicity":
        if len(expressions) == 1:
            from .real_analysis.monotonicity import monotonicity

            base_domain = S.Reals if domain is True else domain
            restricted = restrict_domain(base_domain, cond, var)
            inner_direction = monotonicity(expressions[0], var, domain=restricted)
            result = _compose_monotonicity(instantiated.result, inner_direction)
            complete = result is not Monotonicity.UNKNOWN
            return TransportedPropertyRule(
                instantiated,
                (var,),
                subject_map,
                cond,
                result,
                complete,
                evidence=(("inner-monotonicity", inner_direction),),
            )

        from .real_analysis.sign import sign
        from .registry_derivatives import registry_derivative

        cert = registry_derivative(instantiated.expression, var, registry=registry)
        if not cert.complete:
            return TransportedPropertyRule(
                instantiated,
                (var,),
                subject_map,
                sp.And(cond, cert.conditions),
                Monotonicity.UNKNOWN,
                False,
                evidence=(("chain-rule", cert),),
            )
        combined_condition = sp.And(cond, cert.conditions)
        base_domain = S.Reals if domain is True else domain
        restricted = restrict_domain(base_domain, combined_condition, var)
        derivative_sign = sign(cert.derivative, var, domain=restricted)
        result = _monotonicity_from_derivative_sign(derivative_sign)
        return TransportedPropertyRule(
            instantiated,
            (var,),
            subject_map,
            combined_condition,
            result,
            result is not Monotonicity.UNKNOWN,
            evidence=(
                ("chain-rule", cert),
                ("path-derivative-sign", derivative_sign),
            ),
        )

    return TransportedPropertyRule(
        instantiated,
        (var,),
        subject_map,
        cond,
        instantiated.result,
        False,
        evidence=(("unsupported-property-transport", prop),),
    )


def _constant_matrix_convexity(matrix):
    matrix = sp.ImmutableMatrix(matrix)
    if matrix.rows != matrix.cols or matrix != matrix.T:
        return Convexity.UNKNOWN
    if matrix == sp.zeros(matrix.rows):
        return Convexity.AFFINE
    try:
        if matrix.is_positive_definite is True:
            return Convexity.STRONGLY_CONVEX
        if matrix.is_negative_definite is True:
            return Convexity.STRONGLY_CONCAVE
        if matrix.is_positive_semidefinite is True:
            return Convexity.CONVEX
        if matrix.is_negative_semidefinite is True:
            return Convexity.CONCAVE
    except (TypeError, ValueError, NotImplementedError):
        return Convexity.UNKNOWN
    return Convexity.UNKNOWN


def _multivariate_inner_regularity(prop, expressions, variables):
    """Conservative structural regularity check for multivariate inner maps.

    The checker certifies the inexpensive exact class (polynomial/rational-without
    denominator zeros, exp/trig and algebraic compositions whose derivatives
    are symbolically finite).  Anything unclear remains incomplete.
    """
    evidence = []
    for inner in expressions:
        _, denominator = sp.fraction(sp.together(inner))
        if denominator != 1 and denominator.has(*variables):
            evidence.append((inner, "variable-denominator"))
            return False, tuple(evidence)
        try:
            derivatives = [sp.diff(inner, v) for v in variables]
        except (TypeError, ValueError, NotImplementedError):
            return False, tuple(evidence)
        bad = any(d.has(sp.Derivative) for d in derivatives)
        evidence.append((inner, tuple(derivatives)))
        if bad:
            return False, tuple(evidence)
        if prop in {"smoothness", "analyticity"} and any(
            derivative.has(sp.Derivative) for derivative in derivatives
        ):
            return False, tuple(evidence)
    return True, tuple(evidence)


def _transport_multivar_rule(instantiated, variables, *, domain=True, registry=None):
    variables = tuple(sp.sympify(v) for v in variables)
    subject_map = SubjectMap.from_instantiated(instantiated, variables)
    prop = instantiated.rule.property.lower()
    cond = instantiated.conditions

    if prop == "sign":
        return TransportedPropertyRule(
            instantiated,
            variables,
            subject_map,
            cond,
            instantiated.result,
            True,
            source="native-registry-multivariate-transport",
            evidence=(("pointwise-multivariate-composition", subject_map.expressions),),
        )

    if prop in {"continuity", "differentiability", "smoothness", "analyticity"}:
        ok, evidence = _multivariate_inner_regularity(
            prop, subject_map.expressions, variables
        )
        return TransportedPropertyRule(
            instantiated,
            variables,
            subject_map,
            cond,
            instantiated.result if ok else None,
            ok,
            source="native-registry-multivariate-transport",
            evidence=(("inner-multivariate-regularity", evidence),),
        )

    if prop == "convexity":
        from .registry_derivatives import registry_chain_rule

        cert = registry_chain_rule(
            instantiated.expression, variables, registry=registry, second_order=True
        )
        if not cert.complete or cert.hessian is None:
            return TransportedPropertyRule(
                instantiated,
                variables,
                subject_map,
                sp.And(cond, cert.conditions),
                Convexity.UNKNOWN,
                False,
                source="native-registry-multivariate-transport",
                evidence=(("multivariate-chain-rule", cert),),
            )
        if any(entry.has(*variables) for entry in cert.hessian):
            result = Convexity.UNKNOWN
        else:
            result = _constant_matrix_convexity(cert.hessian)
        return TransportedPropertyRule(
            instantiated,
            variables,
            subject_map,
            sp.And(cond, cert.conditions),
            result,
            result is not Convexity.UNKNOWN,
            source="native-registry-multivariate-transport",
            evidence=(("multivariate-chain-rule", cert),),
        )

    return TransportedPropertyRule(
        instantiated,
        variables,
        subject_map,
        cond,
        instantiated.result,
        False,
        source="native-registry-multivariate-transport",
        evidence=(("unsupported-multivariate-property-transport", prop),),
    )


def transport_instantiated_rule(
    instantiated: InstantiatedPropertyRule,
    analysis_var,
    *,
    domain=True,
    registry=None,
) -> TransportedPropertyRule:
    """Transport a native theorem to one or more analysis variables.

    A single symbol uses scalar-path transport. A tuple/list of symbols invokes
    genuinely multivariate transport.
    """
    if isinstance(analysis_var, (tuple, list, sp.Tuple)):
        variables = tuple(analysis_var)
        if len(variables) == 1:
            return _transport_scalar_rule(
                instantiated, variables[0], domain=domain, registry=registry
            )
        return _transport_multivar_rule(
            instantiated, variables, domain=domain, registry=registry
        )
    return _transport_scalar_rule(
        instantiated, analysis_var, domain=domain, registry=registry
    )
