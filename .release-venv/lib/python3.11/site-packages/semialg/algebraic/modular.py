"""Certified modular accelerators for exact algebraic arithmetic.

The modular computation is never trusted by itself.  Reconstructed candidates are
converted back to exact ``QQ`` coefficients and independently verified before
a result is marked certified.  Factorization uses exact multiplication,
Groebner bases use exact Buchberger/two-way membership checks, and resultants
use exact fixed-degree Sylvester determinants on a degree-complete grid.
"""

from __future__ import annotations

from collections.abc import Sequence
from math import lcm

import sympy as sp
from sympy.polys.domains import ZZ
from sympy.polys.factortools import dup_zz_factor

from ._modular_types import (
    ModularFactorizationResult,
    ModularFractionFieldGroebnerCertificate,
    ModularFractionFieldGroebnerResult,
    ModularGroebnerCertificate,
    ModularGroebnerResult,
    ModularResultantCertificate,
    ModularResultantResult,
    ModularSubresultantCertificate,
    ModularSubresultantResult,
)


def _mul_low(left: Sequence[sp.Rational], right: Sequence[sp.Rational]):
    out = [sp.Rational(0)] * (len(left) + len(right) - 1)
    for i, a in enumerate(left):
        for j, b in enumerate(right):
            out[i + j] += a * b
    while out and out[-1] == 0:
        out.pop()
    return tuple(out)


def _pow_low(poly: Sequence[sp.Rational], exponent: int):
    out = (sp.Rational(1),)
    for _ in range(exponent):
        out = _mul_low(out, poly)
    return out


def factor_univariate_qq_modular(
    coeffs: Sequence[object],
) -> ModularFactorizationResult:
    """Factor a nonzero polynomial in ``QQ[z]`` via modular ``ZZ`` factoring.

    Denominators are cleared exactly, SymPy's modular integer factorizer is
    used for the primitive integer polynomial, and the result is normalized to
    monic rational factors.  Exact reconstruction over ``QQ`` is mandatory.
    """
    qcoeffs = [sp.Rational(c) for c in coeffs]
    while qcoeffs and qcoeffs[-1] == 0:
        qcoeffs.pop()
    if not qcoeffs:
        raise ValueError("cannot factor the zero polynomial")
    if len(qcoeffs) == 1:
        return ModularFactorizationResult(qcoeffs[0], tuple(), True)

    den = 1
    for c in qcoeffs:
        den = lcm(den, int(c.q))
    integers = [int(c * den) for c in qcoeffs]
    content = 0
    for value in integers:
        content = int(sp.igcd(content, value))
    if content == 0:
        raise ValueError("cannot factor the zero polynomial")
    primitive_low = [value // content for value in integers]
    if primitive_low[-1] < 0:
        primitive_low = [-value for value in primitive_low]
        content = -content

    zz_unit, zz_factors = dup_zz_factor(list(reversed(primitive_low)), ZZ)
    scalar = sp.Rational(content * int(zz_unit), den)
    factors: list[tuple[tuple[sp.Rational, ...], int]] = []
    unit = scalar
    for factor_high, multiplicity in zz_factors:
        factor_low = [sp.Rational(v) for v in reversed(factor_high)]
        lc = factor_low[-1]
        monic = tuple(sp.cancel(v / lc) for v in factor_low)
        unit *= lc**multiplicity
        factors.append((monic, int(multiplicity)))

    reconstructed = (sp.Rational(unit),)
    for factor, multiplicity in factors:
        reconstructed = _mul_low(reconstructed, _pow_low(factor, multiplicity))
    target = tuple(qcoeffs)
    certified = reconstructed == target
    if not certified:
        raise ArithmeticError("modular factor reconstruction failed exact verification")
    return ModularFactorizationResult(sp.Rational(unit), tuple(factors), True)


def _prime_stream():
    # Deterministic 30-bit and smaller primes keep finite-field arithmetic fast
    # while making CRT growth substantial after only a few lucky images.
    seeds = (32003, 32009, 32027, 32029, 32051, 32057, 32059, 32063)
    for p in seeds:
        yield p
    p = 32069
    while True:
        p = int(sp.nextprime(p))
        yield p


def _qq_generators(generators, variables):
    vars_tuple = tuple(variables)
    exprs = tuple(sp.expand(sp.sympify(g)) for g in generators if sp.sympify(g) != 0)
    if not exprs:
        return tuple()
    try:
        return tuple(sp.Poly(g, *vars_tuple, domain=sp.QQ).as_expr() for g in exprs)
    except (sp.PolynomialError, sp.polys.polyerrors.CoercionFailed, ValueError, TypeError):
        return None


def _basis_image(exprs, variables, order, prime):
    integer_exprs = []
    for expr in exprs:
        poly = sp.Poly(expr, *variables, domain=sp.QQ)
        _, integer_poly = poly.clear_denoms(convert=True)
        integer_exprs.append(integer_poly.as_expr())
    try:
        gb = sp.groebner(integer_exprs, *variables, order=order, modulus=prime)
    except (sp.PolynomialError, sp.polys.polyerrors.CoercionFailed, ValueError, TypeError):
        return None
    image = {}
    for poly in gb.polys:
        lm = tuple(int(v) for v in poly.LM(order=gb.order).exponents)
        coeffs = {}
        for monom, coeff in poly.terms(order=gb.order):
            coeffs[tuple(int(v) for v in monom)] = int(coeff) % prime
        image[lm] = coeffs
    signature = tuple(sorted((lm, tuple(sorted(coeffs))) for lm, coeffs in image.items()))
    return signature, image


def _crt_pair(a, modulus, b, prime):
    step = ((int(b) - int(a)) * pow(int(modulus), -1, int(prime))) % int(prime)
    return int(a) + int(modulus) * step


def _rational_reconstruct(residue: int, modulus: int):
    """Wang rational reconstruction, returning ``None`` before uniqueness."""
    residue %= modulus
    r0, r1 = int(modulus), int(residue)
    s0, s1 = 0, 1
    # Integer form of |r|, |s| < sqrt(m/2): 2*r^2 < m.
    while r1 and 2 * r1 * r1 >= modulus:
        q = r0 // r1
        r0, r1 = r1, r0 - q * r1
        s0, s1 = s1, s0 - q * s1
    if not s1 or 2 * s1 * s1 >= modulus:
        return None
    if s1 < 0:
        r1, s1 = -r1, -s1
    if sp.igcd(r1, s1) != 1:
        return None
    q = sp.Rational(r1, s1)
    if int(q.p) * pow(int(q.q), -1, modulus) % modulus != residue:
        return None
    return q


def _reconstruct_basis(signature, residues, modulus, variables):
    polys = []
    for lm, support in signature:
        expr = sp.Integer(0)
        for monom in support:
            coeff = _rational_reconstruct(residues[lm][monom], modulus)
            if coeff is None:
                return None
            term = coeff
            for variable, exponent in zip(variables, monom, strict=False):
                term *= variable**exponent
            expr += term
        polys.append(sp.expand(expr))
    return tuple(polys)


def _monomial_exponents(variable_count: int, max_degree: int):
    out = []

    def rec(prefix, remaining, index):
        if index == variable_count - 1:
            out.append(tuple(prefix + [remaining]))
            return
        for value in range(remaining + 1):
            rec(prefix + [value], remaining - value, index + 1)

    if variable_count == 0:
        return (tuple(),)
    for degree in range(max_degree + 1):
        rec([], degree, 0)
    return tuple(out)


def _poly_dict(expr, variables):
    poly = sp.Poly(expr, *variables, domain=sp.QQ)
    return {tuple(int(v) for v in monom): sp.Rational(coeff) for monom, coeff in poly.terms()}


def _shift_monomial_dict(coeffs, shift):
    return {
        tuple(a + b for a, b in zip(monom, shift, strict=False)): coeff
        for monom, coeff in coeffs.items()
    }


def _membership_representation(target, generators, variables, *, extra_degree=3):
    """Find an exact Macaulay representation of ``target`` in ``generators``.

    This is a verifier, not a completeness procedure.  Failure to find a
    representation within the small exact search window simply rejects the
    modular acceleration and lets the caller use its ordinary exact path.
    """
    target_poly = sp.Poly(target, *variables, domain=sp.QQ)
    source_polys = tuple(sp.Poly(g, *variables, domain=sp.QQ) for g in generators)
    if target_poly.is_zero:
        return tuple(sp.Integer(0) for _ in generators)
    source_degree = max((p.total_degree() for p in source_polys), default=0)
    start = max(target_poly.total_degree(), source_degree)
    target_dict = _poly_dict(target, variables)
    source_dicts = tuple(_poly_dict(g, variables) for g in generators)
    n = len(variables)
    for degree_cap in range(start, start + extra_degree + 1):
        columns = []
        descriptors = []
        row_support = set(target_dict)
        for j, source in enumerate(source_polys):
            allowance = degree_cap - source.total_degree()
            if allowance < 0:
                continue
            for shift in _monomial_exponents(n, allowance):
                col = _shift_monomial_dict(source_dicts[j], shift)
                columns.append(col)
                descriptors.append((j, shift))
                row_support.update(col)
        if not columns:
            continue
        rows = tuple(sorted(row_support))
        A = sp.zeros(len(rows), len(columns))
        b = sp.zeros(len(rows), 1)
        row_index = {monom: i for i, monom in enumerate(rows)}
        for monom, coeff in target_dict.items():
            b[row_index[monom], 0] = coeff
        for col_index, col in enumerate(columns):
            for monom, coeff in col.items():
                A[row_index[monom], col_index] = coeff
        try:
            solution, params = A.gauss_jordan_solve(b)
        except ValueError:
            continue
        if params.rows:
            substitutions = {params[i, 0]: sp.Integer(0) for i in range(params.rows)}
            solution = solution.xreplace(substitutions)
        multipliers = [sp.Integer(0) for _ in generators]
        for value, (j, shift) in zip(solution, descriptors, strict=False):
            if value == 0:
                continue
            monomial = sp.Integer(1)
            for variable, exponent in zip(variables, shift, strict=False):
                monomial *= variable**exponent
            multipliers[j] += sp.Rational(value) * monomial
        check = sp.expand(
            target - sum(h * g for h, g in zip(multipliers, generators, strict=False))
        )
        if check == 0:
            return tuple(sp.expand(h) for h in multipliers)
    return None


def _canonical_exact_basis(polys, variables, order):
    gb = sp.groebner(polys, *variables, order=order, domain=sp.QQ)
    return gb, tuple(sp.expand(poly.as_expr()) for poly in gb.polys)


def _certify_reconstructed_groebner(source, candidate, variables, order, primes, modulus):
    # Exact Buchberger/reduction work is cheap when candidate is already the
    # reduced basis; importantly, this does not compute a basis of the source.
    gb, canonical = _canonical_exact_basis(candidate, variables, order)
    normalized_candidate = tuple(sp.expand(p) for p in candidate)
    if set(canonical) != set(normalized_candidate) or len(canonical) != len(normalized_candidate):
        return None
    for generator in source:
        _, remainder = gb.reduce(generator)
        if sp.expand(remainder) != 0:
            return None
    representations = []
    for polynomial in canonical:
        rep = _membership_representation(polynomial, source, variables)
        if rep is None:
            return None
        representations.append(rep)
    certificate = ModularGroebnerCertificate(
        tuple(variables),
        order,
        tuple(source),
        canonical,
        tuple(representations),
        tuple(primes),
        int(modulus),
        True,
    )
    return ModularGroebnerResult(canonical, certificate)


def modular_groebner_basis_qq(
    generators: Sequence[object],
    variables: Sequence[sp.Symbol],
    *,
    order: str = "grevlex",
    min_lucky_primes: int = 2,
    max_primes: int = 12,
) -> ModularGroebnerResult | None:
    """Propose a ``QQ`` Groebner basis modularly and certify it exactly.

    Lucky finite-field images with identical reduced-basis support are merged
    coefficientwise by CRT.  Rational reconstruction is attempted as the CRT
    modulus grows.  A proposal is accepted only after exact Groebner checking,
    exact reduction of every source generator, and exact Macaulay membership
    certificates expressing every reconstructed basis polynomial in the source
    ideal.  ``None`` means only that acceleration was inconclusive.
    """
    vars_tuple = tuple(variables)
    if not vars_tuple:
        return None
    source = _qq_generators(generators, vars_tuple)
    if source is None or not source:
        return None

    signature = None
    residues = None
    modulus = 1
    lucky_primes = []
    attempted = 0
    for prime in _prime_stream():
        if attempted >= max_primes:
            break
        attempted += 1
        image_data = _basis_image(source, vars_tuple, order, prime)
        if image_data is None:
            continue
        image_signature, image = image_data
        if signature is None or image_signature != signature:
            # Prefer the larger leading/support pattern when an early prime is
            # unlucky.  Equal-sized conflicts restart conservatively.
            if signature is None or len(image_signature) > len(signature):
                signature = image_signature
                residues = {lm: dict(coeffs) for lm, coeffs in image.items()}
                modulus = prime
                lucky_primes = [prime]
            elif len(image_signature) == len(signature):
                signature = image_signature
                residues = {lm: dict(coeffs) for lm, coeffs in image.items()}
                modulus = prime
                lucky_primes = [prime]
            continue
        if residues is None:
            return None
        for lm, support in signature:
            for monom in support:
                residues[lm][monom] = _crt_pair(
                    residues[lm][monom], modulus, image[lm][monom], prime
                )
        modulus *= prime
        lucky_primes.append(prime)
        if len(lucky_primes) < min_lucky_primes:
            continue
        candidate = _reconstruct_basis(signature, residues, modulus, vars_tuple)
        if candidate is None:
            continue
        certified = _certify_reconstructed_groebner(
            source, candidate, vars_tuple, order, lucky_primes, modulus
        )
        if certified is not None:
            return certified
    return None


def _qq_to_mod(value, prime: int) -> int:
    q = sp.Rational(value)
    den = int(q.q) % prime
    if den == 0:
        raise ZeroDivisionError
    return (int(q.p) % prime) * pow(den, -1, prime) % prime


def _frac_element_mod(expr, qq_field, ff_field, prime: int):
    value = qq_field.convert(sp.cancel(expr))

    def map_poly(poly):
        out = sp.Integer(0)
        for monom, coeff in poly.to_dict().items():
            term = sp.Integer(_qq_to_mod(coeff, prime))
            for parameter, exponent in zip(qq_field.symbols, monom, strict=True):
                term *= parameter ** int(exponent)
            out += term
        return out

    numerator = map_poly(value.numer)
    denominator = map_poly(value.denom)
    if denominator == 0:
        raise ZeroDivisionError
    return ff_field.convert(numerator) / ff_field.convert(denominator)


def _fraction_field_source_image(source, variables, parameters, prime):
    qq_field = sp.QQ.frac_field(*parameters)
    ff_field = sp.GF(prime).frac_field(*parameters)
    images = []
    for expr in source:
        poly = sp.Poly(expr, *variables, domain=qq_field)
        image = sp.Integer(0)
        for monom, coeff in poly.terms():
            mapped = _frac_element_mod(qq_field.to_sympy(coeff), qq_field, ff_field, prime)
            term = ff_field.to_sympy(mapped)
            for variable, exponent in zip(variables, monom, strict=True):
                term *= variable ** int(exponent)
            image += term
        images.append(sp.cancel(image))
    return tuple(images), ff_field


def _normalized_ff_coeff_signature(expr, parameters, field):
    value = field.convert(sp.cancel(expr))
    numerator = value.numer
    denominator = value.denom
    terms = denominator.terms()
    if not terms:
        raise ZeroDivisionError
    lead_coeff = terms[0][1]
    inv = 1 / lead_coeff
    numerator = numerator * inv
    denominator = denominator * inv
    num = {
        tuple(int(e) for e in mon): int(coeff) % field.domain.mod
        for mon, coeff in numerator.terms()
    }
    den = {
        tuple(int(e) for e in mon): int(coeff) % field.domain.mod
        for mon, coeff in denominator.terms()
    }
    signature = (tuple(sorted(num)), tuple(sorted(den)))
    return signature, num, den


def _fraction_field_basis_image(source, variables, parameters, order, prime):
    try:
        images, field = _fraction_field_source_image(source, variables, parameters, prime)
        gb = sp.groebner(images, *variables, order=order, domain=field)
    except (
        ArithmeticError,
        ValueError,
        TypeError,
        ZeroDivisionError,
        sp.PolynomialError,
        sp.polys.polyerrors.CoercionFailed,
    ):
        return None
    data = {}
    signatures = []
    for poly in gb.polys:
        lm = tuple(int(v) for v in poly.LM(order=gb.order).exponents)
        coeff_data = {}
        support_signature = []
        for monom, coeff in poly.terms(order=gb.order):
            mon = tuple(int(v) for v in monom)
            sig, num, den = _normalized_ff_coeff_signature(field.to_sympy(coeff), parameters, field)
            coeff_data[mon] = (num, den)
            support_signature.append((mon, sig))
        data[lm] = coeff_data
        signatures.append((lm, tuple(support_signature)))
    return tuple(signatures), data


def _merge_poly_residues(old, new, modulus, prime):
    out = dict(old)
    for monom in out:
        out[monom] = _crt_pair(out[monom], modulus, new[monom], prime)
    return out


def _reconstruct_parameter_poly(residues, modulus, parameters):
    expr = sp.Integer(0)
    for monom, residue in residues.items():
        coeff = _rational_reconstruct(residue, modulus)
        if coeff is None:
            return None
        term = coeff
        for parameter, exponent in zip(parameters, monom, strict=True):
            term *= parameter**exponent
        expr += term
    return sp.expand(expr)


def _reconstruct_fraction_basis(signature, residues, modulus, variables, parameters):
    polys = []
    for lm, support_signature in signature:
        expr = sp.Integer(0)
        for monom, _coeff_signature in support_signature:
            num_res, den_res = residues[lm][monom]
            numerator = _reconstruct_parameter_poly(num_res, modulus, parameters)
            denominator = _reconstruct_parameter_poly(den_res, modulus, parameters)
            if numerator is None or denominator is None or denominator == 0:
                return None
            coeff = sp.cancel(numerator / denominator)
            term = coeff
            for variable, exponent in zip(variables, monom, strict=True):
                term *= variable**exponent
            expr += term
        polys.append(sp.cancel(expr))
    return tuple(polys)


def _ff_membership_representation(target, generators, variables, parameters, *, extra_degree=3):
    domain = sp.QQ.frac_field(*parameters)
    target_poly = sp.Poly(target, *variables, domain=domain)
    source_polys = tuple(sp.Poly(g, *variables, domain=domain) for g in generators)
    if target_poly.is_zero:
        return tuple(sp.Integer(0) for _ in generators)
    source_degree = max((p.total_degree() for p in source_polys), default=0)
    start = max(target_poly.total_degree(), source_degree)
    n = len(variables)

    def poly_dict(poly):
        return {
            tuple(int(v) for v in mon): sp.cancel(domain.to_sympy(c)) for mon, c in poly.terms()
        }

    target_dict = poly_dict(target_poly)
    source_dicts = tuple(poly_dict(p) for p in source_polys)
    for degree_cap in range(start, start + extra_degree + 1):
        columns, descriptors = [], []
        row_support = set(target_dict)
        for j, source_poly in enumerate(source_polys):
            allowance = degree_cap - source_poly.total_degree()
            if allowance < 0:
                continue
            for shift in _monomial_exponents(n, allowance):
                col = _shift_monomial_dict(source_dicts[j], shift)
                columns.append(col)
                descriptors.append((j, shift))
                row_support.update(col)
        if not columns:
            continue
        rows = tuple(sorted(row_support))
        row_index = {m: i for i, m in enumerate(rows)}
        matrix = sp.zeros(len(rows), len(columns))
        rhs = sp.zeros(len(rows), 1)
        for mon, coeff in target_dict.items():
            rhs[row_index[mon], 0] = coeff
        for col_index, col in enumerate(columns):
            for mon, coeff in col.items():
                matrix[row_index[mon], col_index] = coeff
        try:
            solution, params = matrix.gauss_jordan_solve(rhs)
        except ValueError:
            continue
        if params.rows:
            solution = solution.xreplace({params[i, 0]: sp.Integer(0) for i in range(params.rows)})
        multipliers = [sp.Integer(0) for _ in generators]
        for value, (j, shift) in zip(solution, descriptors, strict=False):
            if value == 0:
                continue
            monomial = sp.prod(v**e for v, e in zip(variables, shift, strict=True))
            multipliers[j] += sp.cancel(value) * monomial
        if (
            sp.cancel(target - sum(h * g for h, g in zip(multipliers, generators, strict=True)))
            == 0
        ):
            return tuple(sp.cancel(h) for h in multipliers)
    return None


def _certify_fraction_field_groebner(
    source, candidate, variables, parameters, order, primes, modulus
):
    domain = sp.QQ.frac_field(*parameters)
    try:
        gb = sp.groebner(candidate, *variables, order=order, domain=domain)
    except (
        ArithmeticError,
        ValueError,
        TypeError,
        sp.PolynomialError,
        sp.polys.polyerrors.CoercionFailed,
    ):
        return None
    canonical = tuple(sp.cancel(p.as_expr()) for p in gb.polys)
    normalized = tuple(sp.cancel(p) for p in candidate)
    if canonical != normalized:
        return None
    for generator in source:
        _, remainder = gb.reduce(generator)
        if sp.cancel(remainder) != 0:
            return None
    representations = []
    for polynomial in canonical:
        rep = _ff_membership_representation(polynomial, source, variables, parameters)
        if rep is None:
            return None
        representations.append(rep)
    cert = ModularFractionFieldGroebnerCertificate(
        tuple(variables),
        tuple(parameters),
        order,
        tuple(source),
        canonical,
        tuple(representations),
        tuple(primes),
        int(modulus),
        True,
    )
    return ModularFractionFieldGroebnerResult(canonical, cert)


def _fraction_signature_rank(signature) -> tuple[int, int]:
    coefficient_support = 0
    for _lm, polynomial_support in signature:
        for _monom, (numerator_support, denominator_support) in polynomial_support:
            coefficient_support += len(numerator_support) + len(denominator_support)
    return len(signature), coefficient_support


def modular_groebner_basis_fraction_field(
    generators: Sequence[object],
    variables: Sequence[sp.Symbol],
    parameters: Sequence[sp.Symbol],
    *,
    order: str = "grevlex",
    min_lucky_primes: int = 2,
    max_primes: int = 12,
) -> ModularFractionFieldGroebnerResult | None:
    """Certified modular Groebner proposal over ``QQ(parameters)``.

    Finite-field rational-function images are normalized coefficientwise,
    merged by CRT, reconstructed over QQ, and independently verified in the
    exact fraction field.  Unlucky primes or unstable supports are discarded.
    """
    vars_ = tuple(variables)
    params = tuple(parameters)
    if not vars_ or not params:
        return None
    domain = sp.QQ.frac_field(*params)
    try:
        source = tuple(
            sp.cancel(sp.Poly(g, *vars_, domain=domain).as_expr())
            for g in generators
            if sp.cancel(g) != 0
        )
    except (ValueError, TypeError, sp.PolynomialError, sp.polys.polyerrors.CoercionFailed):
        return None
    if not source:
        return None
    signature = None
    residues = None
    modulus = 1
    lucky = []
    attempted = 0
    for prime in _prime_stream():
        if attempted >= max_primes:
            break
        attempted += 1
        image_data = _fraction_field_basis_image(source, vars_, params, order, prime)
        if image_data is None:
            continue
        image_signature, image = image_data
        if signature != image_signature:
            if signature is not None and _fraction_signature_rank(
                image_signature
            ) < _fraction_signature_rank(signature):
                continue
            signature = image_signature
            residues = {
                lm: {monom: (dict(num), dict(den)) for monom, (num, den) in coeffs.items()}
                for lm, coeffs in image.items()
            }
            modulus = prime
            lucky = [prime]
            continue
        if residues is None:
            return None
        for lm, support_signature in signature:
            for monom, _ in support_signature:
                old_num, old_den = residues[lm][monom]
                new_num, new_den = image[lm][monom]
                residues[lm][monom] = (
                    _merge_poly_residues(old_num, new_num, modulus, prime),
                    _merge_poly_residues(old_den, new_den, modulus, prime),
                )
        modulus *= prime
        lucky.append(prime)
        if len(lucky) < min_lucky_primes:
            continue
        candidate = _reconstruct_fraction_basis(signature, residues, modulus, vars_, params)
        if candidate is None:
            continue
        certified = _certify_fraction_field_groebner(
            source, candidate, vars_, params, order, lucky, modulus
        )
        if certified is not None:
            return certified
    return None


__all__ = [
    "ModularFactorizationResult",
    "factor_univariate_qq_modular",
    "ModularGroebnerCertificate",
    "ModularGroebnerResult",
    "modular_groebner_basis_qq",
    "verify_modular_groebner_certificate",
    "ModularFractionFieldGroebnerCertificate",
    "ModularFractionFieldGroebnerResult",
    "modular_groebner_basis_fraction_field",
    "verify_modular_fraction_field_groebner_certificate",
    "ModularResultantCertificate",
    "ModularResultantResult",
    "modular_resultant_qq",
    "verify_modular_resultant_certificate",
    "ModularSubresultantCertificate",
    "ModularSubresultantResult",
    "modular_subresultants_qq",
    "verify_modular_subresultant_certificate",
]

# Certificate replay is outside the computational hot paths.
from ._modular_certification import (  # noqa: E402, F401
    verify_modular_fraction_field_groebner_certificate,
    verify_modular_groebner_certificate,
    verify_modular_resultant_certificate,
    verify_modular_subresultant_certificate,
)

# Resultant/subresultant realization lives in its dedicated backend.
from ._modular_resultants import modular_resultant_qq, modular_subresultants_qq  # noqa: E402, F401
