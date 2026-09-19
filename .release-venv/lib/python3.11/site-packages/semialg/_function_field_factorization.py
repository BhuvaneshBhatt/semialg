"""Exact univariate polynomial arithmetic and factorization over function fields."""

from __future__ import annotations

from collections.abc import Sequence
from itertools import count, permutations

import sympy as sp

from .algebraic_function_fields import (
    CertifiedFactorization,
    FunctionFieldError,
    MonogenicElement,
    MonogenicFunctionField,
    RationalFunctionField,
)


def _poly_trim(p, field):
    p = list(p)
    while p and field.is_zero(p[-1]):
        p.pop()
    return p


def _poly_degree(p, field):
    return len(_poly_trim(p, field)) - 1


def _poly_add(a, b, field):
    n = max(len(a), len(b))
    out = []
    for i in range(n):
        out.append((a[i] if i < len(a) else field.zero) + (b[i] if i < len(b) else field.zero))
    return _poly_trim(out, field)


def _poly_sub(a, b, field):
    n = max(len(a), len(b))
    out = []
    for i in range(n):
        out.append((a[i] if i < len(a) else field.zero) - (b[i] if i < len(b) else field.zero))
    return _poly_trim(out, field)


def _poly_mul(a, b, field):
    if not a or not b:
        return []
    out = [field.zero] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            out[i + j] = out[i + j] + x * y
    return _poly_trim(out, field)


def _poly_divmod(a, b, field):
    a = _poly_trim(a, field)
    b = _poly_trim(b, field)
    if not b:
        raise ZeroDivisionError
    q = [field.zero] * max(0, len(a) - len(b) + 1)
    r = list(a)
    db = len(b) - 1
    inv = field.inv(b[-1])
    while r and len(r) - 1 >= db:
        k = len(r) - 1 - db
        c = r[-1] * inv
        q[k] = c
        for j in range(db + 1):
            r[k + j] = r[k + j] - c * b[j]
        r = _poly_trim(r, field)
    return _poly_trim(q, field), r


def _poly_monic(p, field):
    p = _poly_trim(p, field)
    if not p:
        return []
    inv = field.inv(p[-1])
    return [c * inv for c in p]


def _poly_gcd(a, b, field):
    a = _poly_trim(a, field)
    b = _poly_trim(b, field)
    while b:
        _q, r = _poly_divmod(a, b, field)
        a, b = b, r
    return _poly_monic(a, field)


def _poly_derivative(p, field):
    return _poly_trim([field.convert(i) * p[i] for i in range(1, len(p))], field)


def _poly_is_squarefree(p, field) -> bool:
    p = _poly_trim(p, field)
    if _poly_degree(p, field) <= 0:
        return True
    derivative = _poly_derivative(p, field)
    return _poly_degree(_poly_gcd(p, derivative, field), field) == 0


def _squarefree_decomposition_with_unit(coeffs, field):
    coeffs = _poly_trim([field.convert(c) for c in coeffs], field)
    degree = _poly_degree(coeffs, field)
    if degree < 0:
        raise FunctionFieldError("cannot decompose the zero polynomial")
    if degree == 0:
        return coeffs[0], ()
    unit = coeffs[-1]
    monic = _poly_monic(coeffs, field)
    return unit, _poly_squarefree_decomposition(monic, field)


def _poly_squarefree_decomposition(p, field):
    """Return the characteristic-zero squarefree decomposition of monic ``p``.

    Factors are returned as ``(monic_factor, multiplicity)`` pairs.  This is
    Yun's gcd-based decomposition specialized to a characteristic-zero field,
    so no p-th-root branch is needed.
    """
    p = _poly_trim(p, field)
    degree = _poly_degree(p, field)
    if degree <= 0:
        return ()
    p = _poly_monic(p, field)
    derivative = _poly_derivative(p, field)
    repeated = _poly_gcd(p, derivative, field)
    squarefree_quotient, remainder = _poly_divmod(p, repeated, field)
    if _poly_degree(remainder, field) >= 0:
        raise FunctionFieldError("inexact squarefree decomposition")

    parts = []
    multiplicity = 1
    w = squarefree_quotient
    c = repeated
    while _poly_degree(w, field) > 0:
        y = _poly_gcd(w, c, field)
        z, remainder = _poly_divmod(w, y, field)
        if _poly_degree(remainder, field) >= 0:
            raise FunctionFieldError("inexact squarefree multiplicity extraction")
        if _poly_degree(z, field) > 0:
            parts.append((tuple(_poly_monic(z, field)), multiplicity))
        w = y
        c, remainder = _poly_divmod(c, y, field)
        if _poly_degree(remainder, field) >= 0:
            raise FunctionFieldError("inexact squarefree gcd quotient")
        multiplicity += 1

    if _poly_degree(c, field) > 0:
        raise FunctionFieldError("squarefree decomposition left a residual factor")
    return tuple(parts)


def _trager_shift_is_good(coeffs, ext: MonogenicFunctionField, shift: int):
    """Return the exact norm when the Trager discriminant is nonzero."""
    shifted = _poly_shift_alpha(coeffs, ext, shift)
    norm = _norm_polynomial(shifted, ext)
    if _poly_is_squarefree(norm, ext.base):
        return norm
    return None


def _first_good_trager_shift(coeffs, ext: MonogenicFunctionField):
    """Choose an integer outside the exact finite bad-shift discriminant set."""
    # Trager's characteristic-zero norm argument proves that for squarefree f
    # the discriminant D(s) of Norm(f(z-s*alpha)) is not identically zero.
    # Consequently this unbounded deterministic search is provably terminating;
    # each membership test D(s) != 0 is exact via squarefreeness of the norm.
    for shift in count():
        norm = _trager_shift_is_good(coeffs, ext, shift)
        if norm is not None:
            return shift, norm
    raise AssertionError("unreachable: a nonzero univariate polynomial has finitely many roots")


def _poly_shift_alpha(coeffs, field: MonogenicFunctionField, shift: int):
    # p(z - shift*alpha), Horner in polynomial ring L[z].
    out = []
    linear = [field.neg(field.convert(shift) * field.alpha), field.one]
    for c in reversed(coeffs):
        out = _poly_mul(out, linear, field)
        if out:
            out[0] = out[0] + c
        else:
            out = [c]
    return _poly_trim(out, field)


def _perm_sign(p):
    inv = sum(p[i] > p[j] for i in range(len(p)) for j in range(i + 1, len(p)))
    return -1 if inv % 2 else 1


def _det_poly(matrix, field):
    n = len(matrix)
    total = []
    for perm in permutations(range(n)):
        term = [field.one]
        for i, j in enumerate(perm):
            term = _poly_mul(term, matrix[i][j], field)
        if _perm_sign(perm) < 0:
            term = [-c for c in term]
        total = _poly_add(total, term, field)
    return total


def _norm_polynomial(coeffs, ext: MonogenicFunctionField):
    # determinant of multiplication by p(z) on L/K as K[z]-module.
    d = ext.degree
    cols = []
    for j in range(d):
        basis = [ext.base.zero] * d
        basis[j] = ext.base.one
        bj = MonogenicElement(ext, tuple(basis))
        col = []
        for c in coeffs:
            prod = ext.mul(c, bj)
            for row in range(d):
                while len(col) <= row:
                    col.append([])
                # handled below
        cols.append(bj)
    matrix = [[[] for _ in range(d)] for _ in range(d)]
    for j in range(d):
        basis = [ext.base.zero] * d
        basis[j] = ext.base.one
        bj = MonogenicElement(ext, tuple(basis))
        for power, c in enumerate(coeffs):
            prod = ext.mul(c, bj)
            for i, bc in enumerate(prod.coeffs):
                while len(matrix[i][j]) <= power:
                    matrix[i][j].append(ext.base.zero)
                matrix[i][j][power] = bc
    return _det_poly(matrix, ext.base)


def _poly_reconstruct_factorization(unit, factors, field):
    out = [field.convert(unit)]
    for factor, multiplicity in factors:
        for _ in range(multiplicity):
            out = _poly_mul(out, list(factor), field)
    return _poly_trim(out, field)


def certified_factor_univariate(field, coeffs: Sequence[object]) -> CertifiedFactorization:
    """Factor and independently verify exact reconstruction in ``field[z]``."""
    target = _poly_trim([field.convert(c) for c in coeffs], field)
    if not target:
        raise FunctionFieldError("cannot factor the zero polynomial")
    unit, factors = field.factor_univariate(target)
    reconstructed = _poly_reconstruct_factorization(unit, factors, field)
    verified = len(reconstructed) == len(target) and all(
        field.is_zero(a - b) for a, b in zip(reconstructed, target, strict=True)
    )
    if not verified:
        raise FunctionFieldError("factorization failed exact reconstruction")
    method = (
        "modular_qq+exact_reconstruction"
        if isinstance(field, RationalFunctionField)
        and all(not sp.sympify(c).free_symbols.intersection(field.parameters) for c in target)
        else "trager_tower+exact_reconstruction"
        if isinstance(field, MonogenicFunctionField)
        else "rational_function_factor+exact_reconstruction"
    )
    return CertifiedFactorization(unit, tuple(factors), True, method)
