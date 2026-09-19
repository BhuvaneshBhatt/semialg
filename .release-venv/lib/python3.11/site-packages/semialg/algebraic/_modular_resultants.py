"""Certified modular resultant and subresultant realization backend."""

from __future__ import annotations

from itertools import product

import sympy as sp

from ._modular_types import (
    ModularResultantCertificate,
    ModularResultantResult,
    ModularSubresultantCertificate,
    ModularSubresultantResult,
)
from .modular import _crt_pair, _prime_stream, _rational_reconstruct


def _input_expr(value):
    return value.as_expr() if isinstance(value, sp.Poly) else sp.sympify(value)


def _all_polynomial_symbols(first, second, variable):
    symbols = (set(_input_expr(first).free_symbols) | set(_input_expr(second).free_symbols)) - {
        variable
    }
    return tuple(sorted(symbols, key=lambda sym: (sym.name, sp.srepr(sym))))


def _qq_poly_expr(expr, variables):
    try:
        return sp.Poly(sp.expand(_input_expr(expr)), *variables, domain=sp.QQ)
    except (sp.PolynomialError, sp.polys.polyerrors.CoercionFailed, ValueError, TypeError):
        return None


def _modularize_qq_expr(expr, variables, prime):
    poly = _qq_poly_expr(expr, variables)
    if poly is None:
        return None
    out = sp.Integer(0)
    for monom, coeff in poly.terms():
        q = sp.Rational(coeff)
        try:
            value = int(q.p) * pow(int(q.q), -1, int(prime)) % int(prime)
        except ValueError:
            return None
        term = sp.Integer(value)
        for symbol, exponent in zip(variables, monom, strict=False):
            term *= symbol**exponent
        out += term
    return sp.expand(out)


def _poly_support_image(expr, variables, prime):
    try:
        poly = sp.Poly(expr, *variables, modulus=prime)
    except (sp.PolynomialError, sp.polys.polyerrors.CoercionFailed, ValueError, TypeError):
        return None
    return {tuple(int(v) for v in monom): int(coeff) % prime for monom, coeff in poly.terms()}


def _crt_support(residues, image, modulus, prime):
    if set(residues) != set(image):
        return None
    return {monom: _crt_pair(residues[monom], modulus, image[monom], prime) for monom in residues}


def _reconstruct_expr_from_support(residues, modulus, variables):
    expr = sp.Integer(0)
    for monom, residue in residues.items():
        coeff = _rational_reconstruct(residue, modulus)
        if coeff is None:
            return None
        term = coeff
        for symbol, exponent in zip(variables, monom, strict=False):
            term *= symbol**exponent
        expr += term
    return sp.expand(expr)


def _parameter_degree_bounds(first, second, variable, parameters):
    all_vars = (variable, *parameters)
    f = sp.Poly(first, *all_vars, domain=sp.QQ)
    g = sp.Poly(second, *all_vars, domain=sp.QQ)
    m = f.degree(variable)
    n = g.degree(variable)
    bounds = []
    for parameter in parameters:
        bounds.append(n * f.degree(parameter) + m * g.degree(parameter))
    return tuple(max(0, int(bound)) for bound in bounds)


def _fixed_degree_coefficients(expr, variable, degree, substitutions):
    specialized = sp.expand(_input_expr(expr).subs(substitutions))
    poly = sp.Poly(specialized, variable, domain=sp.QQ)
    return tuple(sp.Rational(poly.nth(i)) for i in range(degree, -1, -1))


def _sylvester_det_at(first, second, variable, m, n, substitutions):
    f = _fixed_degree_coefficients(first, variable, m, substitutions)
    g = _fixed_degree_coefficients(second, variable, n, substitutions)
    size = m + n
    matrix = sp.zeros(size, size)
    for row in range(n):
        for j, coeff in enumerate(f):
            matrix[row, row + j] = coeff
    for row in range(m):
        for j, coeff in enumerate(g):
            matrix[n + row, row + j] = coeff
    return sp.factor(matrix.det(method="domain-ge"))


def _resultant_grid(parameters, bounds):
    if not parameters:
        return (tuple(),)
    return tuple(product(*(range(bound + 1) for bound in bounds)))


def _resultant_image(first, second, variable, parameters, prime):
    all_vars = (variable, *parameters)
    fexpr = _modularize_qq_expr(first, all_vars, prime)
    gexpr = _modularize_qq_expr(second, all_vars, prime)
    if fexpr is None or gexpr is None:
        return None
    domain = sp.GF(prime).poly_ring(*parameters) if parameters else sp.GF(prime)
    try:
        f = sp.Poly(fexpr, variable, domain=domain)
        g = sp.Poly(gexpr, variable, domain=domain)
        result = f.resultant(g)
        expr = result.as_expr() if isinstance(result, sp.Poly) else sp.sympify(result)
    except (
        ArithmeticError,
        ValueError,
        TypeError,
        NotImplementedError,
        sp.PolynomialError,
        sp.polys.polyerrors.CoercionFailed,
    ):
        return None
    support = (
        _poly_support_image(expr, parameters or (sp.Dummy("_c"),), prime)
        if parameters
        else {tuple(): int(expr) % prime}
    )
    if support is None:
        return None
    return (f.degree(), g.degree(), tuple(sorted(support))), support


def modular_resultant_qq(
    first, second, variable, *, min_lucky_primes=2, max_primes=12, max_verification_points=4096
):
    """Reconstruct a polynomial resultant modularly and certify it exactly."""
    parameters = _all_polynomial_symbols(first, second, variable)
    all_vars = (variable, *parameters)
    if _qq_poly_expr(first, all_vars) is None or _qq_poly_expr(second, all_vars) is None:
        return None
    bounds = _parameter_degree_bounds(first, second, variable, parameters)
    grid = _resultant_grid(parameters, bounds)
    if len(grid) > max_verification_points:
        return None
    signature = None
    residues = None
    modulus = 1
    primes = []
    attempted = 0
    for prime in _prime_stream():
        if attempted >= max_primes:
            break
        attempted += 1
        image = _resultant_image(first, second, variable, parameters, prime)
        if image is None:
            continue
        image_signature, support = image
        if signature is None or image_signature != signature:
            signature = image_signature
            residues = dict(support)
            modulus = prime
            primes = [prime]
            continue
        merged = _crt_support(residues, support, modulus, prime)
        if merged is None:
            signature = image_signature
            residues = dict(support)
            modulus = prime
            primes = [prime]
            continue
        residues = merged
        modulus *= prime
        primes.append(prime)
        if len(primes) < min_lucky_primes:
            continue
        candidate = (
            _reconstruct_expr_from_support(residues, modulus, parameters)
            if parameters
            else sp.Rational(_rational_reconstruct(residues[tuple()], modulus))
            if _rational_reconstruct(residues[tuple()], modulus) is not None
            else None
        )
        if candidate is None:
            continue
        cert = ModularResultantCertificate(
            variable,
            parameters,
            sp.expand(_input_expr(first)),
            sp.expand(_input_expr(second)),
            sp.expand(candidate),
            bounds,
            grid,
            tuple(primes),
            int(modulus),
            True,
        )
        from ._modular_certification import verify_modular_resultant_certificate

        if verify_modular_resultant_certificate(cert):
            return ModularResultantResult(sp.expand(candidate), cert)
    return None


def _subresultant_image(first, second, variable, parameters, prime):
    all_vars = (variable, *parameters)
    fexpr = _modularize_qq_expr(first, all_vars, prime)
    gexpr = _modularize_qq_expr(second, all_vars, prime)
    if fexpr is None or gexpr is None:
        return None
    domain = sp.GF(prime).poly_ring(*parameters) if parameters else sp.GF(prime)
    try:
        f = sp.Poly(fexpr, variable, domain=domain)
        g = sp.Poly(gexpr, variable, domain=domain)
        if f.degree() < g.degree():
            f, g = g, f
        seq = tuple(sp.subresultants(f, g))
    except (
        ArithmeticError,
        ValueError,
        TypeError,
        NotImplementedError,
        sp.PolynomialError,
        sp.polys.polyerrors.CoercionFailed,
    ):
        return None
    vars_all = (variable, *parameters)
    supports = []
    signature = []
    for item in seq:
        expr = item.as_expr() if isinstance(item, sp.Poly) else sp.sympify(item)
        support = _poly_support_image(expr, vars_all, prime)
        if support is None:
            return None
        supports.append(support)
        degree = sp.Poly(expr, variable, domain=domain).degree()
        signature.append((degree, tuple(sorted(support))))
    return tuple(signature), tuple(supports)


def modular_subresultants_qq(
    first, second, variable, *, min_lucky_primes=2, max_primes=12, max_verification_points=4096
):
    """Reconstruct an exact subresultant sequence from finite-field images."""
    parameters = _all_polynomial_symbols(first, second, variable)
    all_vars = (variable, *parameters)
    if _qq_poly_expr(first, all_vars) is None or _qq_poly_expr(second, all_vars) is None:
        return None
    resultant = modular_resultant_qq(
        first,
        second,
        variable,
        min_lucky_primes=min_lucky_primes,
        max_primes=max_primes,
        max_verification_points=max_verification_points,
    )
    if resultant is None:
        return None
    signature = None
    residues = None
    modulus = 1
    primes = []
    attempted = 0
    vars_all = (variable, *parameters)
    for prime in _prime_stream():
        if attempted >= max_primes:
            break
        attempted += 1
        image = _subresultant_image(first, second, variable, parameters, prime)
        if image is None:
            continue
        image_signature, supports = image
        if signature is None or image_signature != signature:
            signature = image_signature
            residues = [dict(s) for s in supports]
            modulus = prime
            primes = [prime]
            continue
        merged_sequence = []
        ok = True
        for old, new in zip(residues, supports, strict=False):
            merged = _crt_support(old, new, modulus, prime)
            if merged is None:
                ok = False
                break
            merged_sequence.append(merged)
        if not ok:
            signature = image_signature
            residues = [dict(s) for s in supports]
            modulus = prime
            primes = [prime]
            continue
        residues = merged_sequence
        modulus *= prime
        primes.append(prime)
        if len(primes) < min_lucky_primes:
            continue
        candidates = tuple(
            _reconstruct_expr_from_support(support, modulus, vars_all) for support in residues
        )
        if any(candidate is None for candidate in candidates):
            continue
        cert = ModularSubresultantCertificate(
            variable,
            parameters,
            sp.expand(_input_expr(first)),
            sp.expand(_input_expr(second)),
            tuple(sp.expand(c) for c in candidates),
            resultant.certificate,
            tuple(primes),
            int(modulus),
            True,
        )
        from ._modular_certification import verify_modular_subresultant_certificate

        if verify_modular_subresultant_certificate(cert):
            return ModularSubresultantResult(
                tuple(sp.expand(c) for c in candidates), resultant.resultant, cert
            )
    return None
