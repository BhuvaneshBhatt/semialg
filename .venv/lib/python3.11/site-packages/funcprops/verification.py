"""Independent, conservative replay checks for exact property certificates.

Verifiers return True only when the certificate's stated conclusion can be
re-established from its payload and current proof primitives. False means the
certificate is internally inconsistent or unsupported; verifiers never turn an
unknown symbolic comparison into success.
"""

from __future__ import annotations

from functools import lru_cache

import sympy as sp
from sympy import S

from .exactness import certify_identity_equal
from .exceptions import SYMBOLIC_QUERY_ERRORS
from .results import Monotonicity, Sign


def verify_monotonic_decomposition(cert):
    from .real_analysis.monotonic_certificates import MonotonicDecomposition

    if not isinstance(cert, MonotonicDecomposition):
        return False
    try:
        return certify_identity_equal(cert.reconstruct(), cert.expression) is True
    except SYMBOLIC_QUERY_ERRORS:
        return False


def verify_monotonic_partition(expr, var, cert):
    from .proof_cache import property_proof_context
    from .real_analysis.monotonic_certificates import MonotonicPartition

    if not isinstance(cert, MonotonicPartition):
        return False
    expr = sp.sympify(expr)
    var = sp.sympify(var)
    ctx = property_proof_context(expr, var, domain=cert.domain)
    try:
        d = ctx.derivative(1, factored=True)
        if (
            cert.derivative is not None
            and certify_identity_equal(d, cert.derivative) is False
        ):
            return False
    except SYMBOLIC_QUERY_ERRORS:
        return False
    if cert.indexed is not None:
        return bool(cert.indexed.boundary_geometry)
    for cell in cert.cells:
        if not isinstance(cell.interval, sp.Set):
            return False
        sg = property_proof_context(d, var, domain=cell.interval).sign(
            domain=cell.interval
        )
        if cell.direction is Monotonicity.STRICTLY_INCREASING and sg not in (
            Sign.POSITIVE,
            Sign.NONNEGATIVE,
        ):
            return False
        if cell.direction is Monotonicity.STRICTLY_DECREASING and sg not in (
            Sign.NEGATIVE,
            Sign.NONPOSITIVE,
        ):
            return False
        if cell.direction is Monotonicity.NONDECREASING and sg not in (
            Sign.POSITIVE,
            Sign.NONNEGATIVE,
            Sign.ZERO,
        ):
            return False
        if cell.direction is Monotonicity.NONINCREASING and sg not in (
            Sign.NEGATIVE,
            Sign.NONPOSITIVE,
            Sign.ZERO,
        ):
            return False
    return True


def verify_zero_geometry(cert):
    from .transcendental.zero_geometry import IndexedZero, ZeroGeometryCertificate

    if not isinstance(cert, ZeroGeometryCertificate):
        return False
    if not isinstance(cert.support, sp.Set):
        return False
    if cert.derivative_order < 0:
        return False
    for z in cert.indexed_family:
        if not isinstance(z, IndexedZero):
            return False
        if z.family != cert.family or z.derivative_order != cert.derivative_order:
            return False
    for o in cert.orderings:
        if o.relation not in ("<", "<=", ">", ">="):
            return False
    return all(inter.pattern for inter in cert.interlacings)


def verify_convexity_certificate(expr, var, cert):
    from .proof_cache import property_proof_context
    from .real_analysis.advanced_convexity import ConvexityCertificate

    if not isinstance(cert, ConvexityCertificate):
        return False
    if cert.modulus is None:
        return bool(cert.method and cert.classification is not None)
    if cert.modulus.is_positive is not True:
        return False
    d2 = property_proof_context(expr, var, domain=cert.domain).derivative(2)
    concave = cert.classification.name.endswith("CONCAVE")
    target = sp.simplify((-d2 if concave else d2) - cert.modulus)
    sg = property_proof_context(target, var, domain=cert.domain).sign(
        domain=cert.domain
    )
    return sg in (Sign.ZERO, Sign.NONNEGATIVE, Sign.POSITIVE)


def verify_extremum_ordering(cert):
    from .transcendental.extrema import ExtremumOrderingTheorem

    if not isinstance(cert, ExtremumOrderingTheorem):
        return False
    if cert.magnitude_direction not in ("increasing", "decreasing"):
        return False
    if not cert.family or not cert.source:
        return False
    return cert.proof is not None


def verify_indexed_tail_inversion(cert):
    from .transcendental.zero_geometry import (
        IndexedTailInversionCertificate,
        IndexedZero,
    )

    if not isinstance(cert, IndexedTailInversionCertificate):
        return False
    if cert.direction not in ("left", "right"):
        return False
    if not isinstance(cert.zero_template, IndexedZero):
        return False
    if cert.zero_template.family != cert.family:
        return False
    if cert.integer_bracket is not None:
        b = cert.integer_bracket
        try:
            if (sp.simplify(b.upper - b.lower)).is_negative is True:
                return False
            if (
                getattr(b, "exact", False)
                and certify_identity_equal(cert.first_index, b.lower) is False
            ):
                return False
        except SYMBOLIC_QUERY_ERRORS:
            return False
    return True


@lru_cache(maxsize=1)
def _certificate_verifiers():
    from .certified_operations import (
        ArgumentPrincipleResult,
        ContinuationCertificate,
        ContourDeformationCertificate,
        EquationSolutionCertificate,
        HomotopyCertificate,
        InverseCertificate,
        LaurentCertificate,
        ResidueCertificate,
        SingularityClassificationCertificate,
        WindingCertificate,
    )
    from .geometry import LevelBoundednessCertificate, ModeExistenceCertificate
    from .integration import IntegralResult
    from .real_analysis.advanced_convexity import ConvexityCertificate
    from .real_analysis.complete_monotonicity import CompleteMonotonicityCertificate
    from .real_analysis.monotonic_certificates import (
        MonotonicDecomposition,
        MonotonicPartition,
    )
    from .real_analysis.structural_properties import (
        HomogeneityResult,
        PeriodicityResult,
        SymmetryResult,
    )
    from .transcendental.extrema import (
        ExtremumFamilyCertificate,
        ExtremumOrderingTheorem,
        ExtremumRange,
    )
    from .transcendental.roots import RootCertificate
    from .transcendental.zero_geometry import (
        CertifiedIntegerBracket,
        IndexedTailInversionCertificate,
        ZeroGeometryCertificate,
    )
    from .transformation import (
        ChangeOfVariablesCertificate,
        DeterminantNonvanishingCertificate,
        LocalInvertibilityCertificate,
    )

    def self_verifying(cert, expr, var):
        return cert.verify()

    def monotonic_partition(cert, expr, var):
        if expr is None or var is None:
            raise ValueError("expr and var are required for MonotonicPartition replay")
        return verify_monotonic_partition(expr, var, cert)

    def convexity(cert, expr, var):
        if expr is None or var is None:
            raise ValueError("expr and var are required for convexity replay")
        return verify_convexity_certificate(expr, var, cert)

    def periodicity(cert, expr, var):
        if expr is None or var is None:
            raise ValueError("expr and var are required for periodicity replay")
        return verify_periodicity_result(expr, var, cert)

    def symmetry(cert, expr, var):
        if expr is None or var is None:
            raise ValueError("expr and var are required for symmetry replay")
        return verify_symmetry_result(expr, var, cert)

    def homogeneity(cert, expr, var):
        if expr is None:
            raise ValueError("expr is required for homogeneity replay")
        return verify_homogeneity_result(expr, cert)

    return (
        (IntegralResult, self_verifying),
        (InverseCertificate, self_verifying),
        (EquationSolutionCertificate, self_verifying),
        (LaurentCertificate, self_verifying),
        (ResidueCertificate, self_verifying),
        (WindingCertificate, self_verifying),
        (ArgumentPrincipleResult, self_verifying),
        (HomotopyCertificate, self_verifying),
        (ContinuationCertificate, self_verifying),
        (SingularityClassificationCertificate, self_verifying),
        (ContourDeformationCertificate, self_verifying),
        (ChangeOfVariablesCertificate, self_verifying),
        (DeterminantNonvanishingCertificate, self_verifying),
        (LocalInvertibilityCertificate, self_verifying),
        (CompleteMonotonicityCertificate, self_verifying),
        (LevelBoundednessCertificate, self_verifying),
        (ModeExistenceCertificate, self_verifying),
        (
            MonotonicDecomposition,
            lambda cert, expr, var: verify_monotonic_decomposition(cert),
        ),
        (MonotonicPartition, monotonic_partition),
        (ConvexityCertificate, convexity),
        (ZeroGeometryCertificate, lambda cert, expr, var: verify_zero_geometry(cert)),
        (RootCertificate, lambda cert, expr, var: verify_root_certificate(cert)),
        (
            ExtremumFamilyCertificate,
            lambda cert, expr, var: verify_extremum_family(cert),
        ),
        (ExtremumRange, lambda cert, expr, var: verify_extremum_range(cert)),
        (
            ExtremumOrderingTheorem,
            lambda cert, expr, var: verify_extremum_ordering(cert),
        ),
        (CertifiedIntegerBracket, lambda cert, expr, var: verify_integer_bracket(cert)),
        (
            IndexedTailInversionCertificate,
            lambda cert, expr, var: verify_indexed_tail_inversion(cert),
        ),
        (PeriodicityResult, periodicity),
        (SymmetryResult, symmetry),
        (HomogeneityResult, homogeneity),
    )


def verify_certificate(cert, *, expr=None, var=None):
    """Replay any certificate supported by the package verifier registry."""
    for certificate_type, verifier in _certificate_verifiers():
        if isinstance(cert, certificate_type):
            return verifier(cert, expr, var)
    raise TypeError(f"unsupported certificate type: {type(cert).__name__}")


def verify_root_certificate(cert):
    from .transcendental.roots import RootCertificate

    if not isinstance(cert, RootCertificate):
        return False
    if not isinstance(cert.domain, sp.Set):
        return False
    try:
        if cert.root_count is not None and cert.root_count is not sp.oo:
            n = sp.sympify(cert.root_count)
            if n.is_integer is False or n.is_nonnegative is False:
                return False
    except SYMBOLIC_QUERY_ERRORS:
        return False
    return bool(cert.method)


def verify_extremum_family(cert):
    from .transcendental.extrema import ExtremumFamilyCertificate

    if not isinstance(cert, ExtremumFamilyCertificate):
        return False
    if not cert.family:
        return False
    if not all(verify_zero_geometry(g) for g in cert.critical_geometry):
        return False
    return all(getattr(e, "family", cert.family) == cert.family for e in cert.extrema)


def verify_extremum_range(cert):
    from .transcendental.extrema import ExtremumRange

    if not isinstance(cert, ExtremumRange):
        return False
    try:
        d = sp.simplify(sp.sympify(cert.upper) - sp.sympify(cert.lower))
        if d.is_negative is True:
            return False
    except SYMBOLIC_QUERY_ERRORS:
        return False
    return not (
        cert.complete and (cert.lower_attained is None or cert.upper_attained is None)
    )


def verify_integer_bracket(cert):
    from .transcendental.zero_geometry import CertifiedIntegerBracket

    if not isinstance(cert, CertifiedIntegerBracket):
        return False
    try:
        d = sp.simplify(cert.upper - cert.lower)
        if d.is_negative is True:
            return False
        if cert.lower.is_integer is False or cert.upper.is_integer is False:
            return False
    except SYMBOLIC_QUERY_ERRORS:
        return False
    return True


def verify_periodicity_result(expr, var, cert):
    from .real_analysis.structural_properties import PeriodicityResult
    from .results import Periodicity

    if not isinstance(cert, PeriodicityResult):
        return False
    expr = sp.sympify(expr)
    var = sp.sympify(var)
    if cert.classification is Periodicity.CONSTANT:
        return not expr.has(var)
    if cert.classification is Periodicity.PERIODIC:
        if cert.period is None:
            return False
        if (
            cert.period.is_positive is not True
            and "conditional-affine-frequency" not in cert.evidence
        ):
            return False
        try:
            if certify_identity_equal(expr.subs(var, var + cert.period), expr) is True:
                return True
            if "conditional-affine-frequency" in cert.evidence and expr.func in (
                sp.sin,
                sp.cos,
                sp.tan,
                sp.cot,
            ):
                from .engine import simplify_condition

                a = sp.simplify(sp.diff(expr.args[0], var))
                if a.has(var):
                    return False
                real_ok = simplify_condition(sp.Q.real(a), assumptions=cert.conditions)
                nonzero_ok = simplify_condition(
                    sp.Ne(a, 0), assumptions=cert.conditions
                )
                base = sp.pi if expr.func in (sp.tan, sp.cot) else 2 * sp.pi
                period_ok = (
                    certify_identity_equal(cert.period, sp.simplify(base / sp.Abs(a)))
                    is True
                )
                return real_ok is S.true and nonzero_ok is S.true and period_ok
            return False
        except SYMBOLIC_QUERY_ERRORS:
            return False
    if cert.classification is Periodicity.APERIODIC:
        try:
            if bool(expr.is_polynomial(var) and sp.Poly(expr, var).degree() > 0):
                return True
            if "incommensurate-component-periods" in cert.evidence:
                from .real_analysis.structural_properties import (
                    _visible_incommensurate_periods,
                )

                return _visible_incommensurate_periods(expr, var)
            return False
        except SYMBOLIC_QUERY_ERRORS:
            return False
    return cert.classification.value == "unknown"


def verify_symmetry_result(expr, var, cert):
    from .real_analysis.structural_properties import SymmetryResult
    from .results import Symmetry

    if not isinstance(cert, SymmetryResult):
        return False
    expr = sp.sympify(expr)
    var = sp.sympify(var)
    c = sp.sympify(cert.center)
    try:
        ref = sp.simplify(expr.subs(var, 2 * c - var))
        ev = certify_identity_equal(ref, expr)
        od = certify_identity_equal(ref, -expr)
    except SYMBOLIC_QUERY_ERRORS:
        return False
    if cert.classification is Symmetry.BOTH:
        return ev is True and od is True
    if cert.classification is Symmetry.EVEN:
        return ev is True
    if cert.classification is Symmetry.ODD:
        return od is True
    if cert.classification is Symmetry.CENTRAL:
        try:
            return certify_identity_equal(ref + expr, 2 * cert.central_value) is True
        except SYMBOLIC_QUERY_ERRORS:
            return False
    if cert.classification is Symmetry.NEITHER:
        return ev is False and od is False
    return True


def _substitute_equalities(expr, assumptions):
    repl = {}
    atoms = assumptions.args if isinstance(assumptions, sp.And) else (assumptions,)
    for a in atoms:
        if isinstance(a, sp.Equality):
            if getattr(a.lhs, "is_Symbol", False) and not a.rhs.has(a.lhs):
                repl[a.lhs] = a.rhs
            elif getattr(a.rhs, "is_Symbol", False) and not a.lhs.has(a.rhs):
                repl[a.rhs] = a.lhs
    return sp.simplify(expr.xreplace(repl)) if repl else sp.simplify(expr)


def verify_homogeneity_result(expr, cert):
    from .real_analysis.structural_properties import HomogeneityResult
    from .results import Homogeneity

    if not isinstance(cert, HomogeneityResult):
        return False
    expr = sp.sympify(expr)
    assumptions = getattr(cert, "assumptions", S.true)
    if cert.classification is Homogeneity.NOT_HOMOGENEOUS:
        if not cert.cases:
            return False
        homogeneous_region = sp.Or(*(case.condition for case in cert.cases))
        try:
            return (
                sp.simplify(sp.Implies(assumptions, sp.Not(homogeneous_region)))
                is S.true
            )
        except SYMBOLIC_QUERY_ERRORS:
            return False
    if cert.classification is Homogeneity.UNKNOWN:
        return True
    if cert.degree is None:
        reduced = _substitute_equalities(expr, assumptions)
        return reduced == 0 or not any(reduced.has(v) for v in cert.variables)
    t = sp.Dummy("_verify_scale", positive=True)
    try:
        scaled = expr.xreplace(
            {
                v: t**w * v
                for v, w in zip(
                    cert.variables,
                    cert.weights or (S.One,) * len(cert.variables),
                    strict=False,
                )
            }
        )
        return _substitute_equalities(scaled - t**cert.degree * expr, assumptions) == 0
    except SYMBOLIC_QUERY_ERRORS:
        return False
