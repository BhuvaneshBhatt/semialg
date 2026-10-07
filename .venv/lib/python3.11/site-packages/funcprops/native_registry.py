from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import sympy as sp

from .registry import FunctionSlot, PropertyRule
from .registry_derivatives import DerivativeRule, IdentityRule
from .results import Analyticity


@dataclass(frozen=True)
class NativeFamilyRegistration:
    """Exact native-registry facts for one function family.

    Registration data stays explicit so CI can count native coverage without
    executing property analyzers. Only facts whose slot semantics are
    unambigous are promoted here.
    """

    family: str
    head: Any
    property_rules: tuple[PropertyRule, ...] = ()
    derivative_rules: tuple[DerivativeRule, ...] = ()
    identity_rules: tuple[IdentityRule, ...] = ()


def _unary_entire(family: str, head, derivative, *, identities=()):
    z = FunctionSlot(0, "z")
    return NativeFamilyRegistration(
        family=family,
        head=head,
        property_rules=(
            PropertyRule(
                head,
                "analyticity",
                z,
                result=Analyticity.ANALYTIC,
                proof_kind="exact-special-function-theorem",
                source=f"native:{family}:entire",
            ),
        ),
        derivative_rules=(
            DerivativeRule(
                head,
                with_respect_to=z,
                result=derivative(z),
                source=f"native:{family}:derivative",
            ),
        ),
        identity_rules=tuple(identities(z) if callable(identities) else identities),
    )


def _derivative_family(family, head, slot, result, *, source=None, **rule_options):
    """Declare a family with one exact argument-derivative theorem."""
    return NativeFamilyRegistration(
        family,
        head,
        derivative_rules=(
            DerivativeRule(
                head,
                with_respect_to=slot,
                result=result,
                source=source or f"native:{family}:argument-derivative",
                **rule_options,
            ),
        ),
    )


def native_family_registrations() -> tuple[NativeFamilyRegistration, ...]:
    """Return the exact native special-function theorem catalogue."""

    registrations = [
        _unary_entire(
            "erf",
            sp.erf,
            lambda z: 2 * sp.exp(-(z**2)) / sp.sqrt(sp.pi),
        ),
        _unary_entire(
            "erfc",
            sp.erfc,
            lambda z: -2 * sp.exp(-(z**2)) / sp.sqrt(sp.pi),
            identities=lambda z: (
                IdentityRule(
                    sp.erfc,
                    1 - sp.erf(z),
                    source="native:erfc:erf-reduction",
                ),
            ),
        ),
        _unary_entire(
            "erfi",
            sp.erfi,
            lambda z: 2 * sp.exp(z**2) / sp.sqrt(sp.pi),
        ),
        _unary_entire(
            "airyai",
            sp.airyai,
            lambda z: sp.airyaiprime(z),
        ),
        _unary_entire(
            "airybi",
            sp.airybi,
            lambda z: sp.airybiprime(z),
        ),
        _unary_entire(
            "airyaiprime",
            sp.airyaiprime,
            lambda z: z * sp.airyai(z),
        ),
        _unary_entire(
            "airybiprime",
            sp.airybiprime,
            lambda z: z * sp.airybi(z),
        ),
    ]

    nu = FunctionSlot(0, "nu")
    z = FunctionSlot(1, "z")
    registrations.append(
        _derivative_family(
            "besselj",
            sp.besselj,
            z,
            (sp.besselj(nu - 1, z) - sp.besselj(nu + 1, z)) / 2,
        )
    )
    registrations.append(
        _derivative_family(
            "besseli",
            sp.besseli,
            z,
            (sp.besseli(nu - 1, z) + sp.besseli(nu + 1, z)) / 2,
        )
    )

    z0 = FunctionSlot(0, "z")
    registrations.append(
        _derivative_family(
            "gamma",
            sp.gamma,
            z0,
            sp.gamma(z0) * sp.polygamma(0, z0),
            source="native:gamma:derivative",
        )
    )

    z = FunctionSlot(0, "z")
    k = FunctionSlot(1, "k")
    registrations.append(
        NativeFamilyRegistration(
            "LambertW",
            sp.LambertW,
            derivative_rules=(
                DerivativeRule(
                    sp.LambertW,
                    with_respect_to=z,
                    result=sp.LambertW(z) / (z * (1 + sp.LambertW(z))),
                    conditions=sp.And(sp.Ne(z, 0), sp.Ne(sp.LambertW(z), -1)),
                    source="native:LambertW:principal-derivative",
                    arity=1,
                ),
                DerivativeRule(
                    sp.LambertW,
                    with_respect_to=z,
                    result=sp.LambertW(z, k) / (z * (1 + sp.LambertW(z, k))),
                    conditions=sp.And(sp.Ne(z, 0), sp.Ne(sp.LambertW(z, k), -1)),
                    source="native:LambertW:branch-derivative",
                    arity=2,
                ),
            ),
        )
    )

    # Orthogonal-polynomial argument derivatives.  Degree/shape arguments are
    # ordinary function slots; only the polynomial variable is differentiated.
    n = FunctionSlot(0, "n")
    z = FunctionSlot(1, "z")
    registrations.extend(
        (
            NativeFamilyRegistration(
                "hermite",
                sp.hermite,
                derivative_rules=(
                    DerivativeRule(
                        sp.hermite,
                        with_respect_to=z,
                        result=2 * n * sp.hermite(n - 1, z),
                        source="native:hermite:argument-derivative",
                    ),
                ),
            ),
            NativeFamilyRegistration(
                "laguerre",
                sp.laguerre,
                derivative_rules=(
                    DerivativeRule(
                        sp.laguerre,
                        with_respect_to=z,
                        result=-sp.assoc_laguerre(n - 1, 1, z),
                        source="native:laguerre:argument-derivative",
                    ),
                ),
            ),
            NativeFamilyRegistration(
                "chebyshevt",
                sp.chebyshevt,
                derivative_rules=(
                    DerivativeRule(
                        sp.chebyshevt,
                        with_respect_to=z,
                        result=n * sp.chebyshevu(n - 1, z),
                        source="native:chebyshevt:argument-derivative",
                    ),
                ),
            ),
        )
    )

    n = FunctionSlot(0, "n")
    a = FunctionSlot(1, "a")
    z = FunctionSlot(2, "z")
    registrations.append(
        NativeFamilyRegistration(
            "gegenbauer",
            sp.gegenbauer,
            derivative_rules=(
                DerivativeRule(
                    sp.gegenbauer,
                    with_respect_to=z,
                    result=2 * a * sp.gegenbauer(n - 1, a + 1, z),
                    source="native:gegenbauer:argument-derivative",
                ),
            ),
        )
    )

    n = FunctionSlot(0, "n")
    a = FunctionSlot(1, "a")
    b = FunctionSlot(2, "b")
    z = FunctionSlot(3, "z")
    registrations.append(
        NativeFamilyRegistration(
            "jacobi",
            sp.jacobi,
            derivative_rules=(
                DerivativeRule(
                    sp.jacobi,
                    with_respect_to=z,
                    result=(n + a + b + 1) * sp.jacobi(n - 1, a + 1, b + 1, z) / 2,
                    source="native:jacobi:argument-derivative",
                ),
            ),
        )
    )

    return tuple(registrations)


def install_native_rules(registry) -> None:
    """Install all exact native rules into ``registry``."""

    for family in native_family_registrations():
        for rule in family.property_rules:
            registry.register_property_rule(rule)
        for rule in family.derivative_rules:
            registry.register_derivative_rule(rule)
        for rule in family.identity_rules:
            registry.register_identity_rule(rule)


def native_registry_snapshot() -> dict[str, Any]:
    """Machine-readable native coverage used by CI's monotonic guard."""

    families = native_family_registrations()
    property_counts: dict[str, int] = {}
    derivative_positions: dict[str, list[int]] = {}
    identity_counts: dict[str, int] = {}
    for fam in families:
        for rule in fam.property_rules:
            property_counts[rule.property] = property_counts.get(rule.property, 0) + 1
        if fam.derivative_rules:
            derivative_positions[fam.family] = sorted(
                {r.with_respect_to.position for r in fam.derivative_rules}
            )
        if fam.identity_rules:
            identity_counts[fam.family] = len(fam.identity_rules)
    return {
        "family_count": len(families),
        "families": sorted(f.family for f in families),
        "property_rule_counts": dict(sorted(property_counts.items())),
        "derivative_positions": dict(sorted(derivative_positions.items())),
        "identity_rule_counts": dict(sorted(identity_counts.items())),
    }
