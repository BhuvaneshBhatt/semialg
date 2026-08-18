from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass
from itertools import product
from math import isqrt

import sympy as sp
from sympy import Eq

from .output_normalization import canon_int_result
from .thue import solve_binary_bounded


@dataclass(frozen=True)
class SpecialFamDesc:
    family_name: str
    variables: tuple[sp.Symbol, ...]
    metadata: dict


def _conjuncts(expr: sp.Expr) -> list[sp.Expr]:
    return list(expr.args) if isinstance(expr, sp.And) else [expr]


def _split_eq(expr: sp.Expr):
    atoms = _conjuncts(expr)
    eqs = [a for a in atoms if isinstance(a, Eq)]
    others = [a for a in atoms if not isinstance(a, Eq)]
    return eqs, others


def detect_sum_fam2(expr: sp.Expr, variables: Sequence[sp.Symbol]) -> SpecialFamDesc | None:
    variables = tuple(variables)
    eqs, others = _split_eq(expr)
    if len(eqs) != 1 or others:
        return None
    diff = sp.expand(eqs[0].lhs - eqs[0].rhs)
    poly = sp.Poly(diff, *variables)
    target = -poly.coeff_monomial(1)
    coeffs = []
    for v in variables:
        c = poly.coeff_monomial(v**2)
        if c not in (0, 1):
            return None
        coeffs.append(int(c))
    if sum(coeffs) == 0:
        return None
    if any(sum(mon) not in (0, 2) for mon, coeff in poly.terms() if coeff != 0):
        return None
    cross_terms = [
        mon for mon, coeff in poly.terms() if coeff != 0 and sum(mon) == 2 and max(mon) < 2
    ]
    if cross_terms:
        return None
    if target.is_integer:
        return SpecialFamDesc(
            "sum_of_squares",
            variables,
            {"target": int(target), "coefficients": tuple(coeffs), "dimension": len(variables)},
        )
    return None


def three_square_obstr(n: int) -> bool:
    while n % 4 == 0 and n > 0:
        n //= 4
    return n % 8 == 7


def solve_sum_of_squares_fam(
    expr: sp.Expr, variables: Sequence[sp.Symbol], *, search_limit: int = 500
):
    desc = detect_sum_fam2(expr, variables)
    if desc is None:
        return None
    n = int(desc.metadata["target"])
    vars_t = desc.variables
    coeffs = desc.metadata["coefficients"]
    if n < 0:
        return canon_int_result(
            vars_t,
            formula=sp.false,
            solutions=[],
            method="sum_of_squares_family",
            complete=True,
            provenance=["special_family"],
        )
    if len(vars_t) == 3 and all(c == 1 for c in coeffs) and three_square_obstr(n):
        return canon_int_result(
            vars_t,
            formula=sp.false,
            solutions=[],
            method="sum_of_squares_legendre_obstruction",
            complete=True,
            provenance=["special_family"],
            metadata=desc.metadata,
        )
    bound = isqrt(n) if n >= 0 else 0
    if (2 * bound + 1) ** len(vars_t) > search_limit**2:
        return None
    pts = []
    for vals in product(range(-bound, bound + 1), repeat=len(vars_t)):
        if sum(c * (v * v) for c, v in zip(coeffs, vals, strict=True)) == n:
            pts.append(vals)
    return canon_int_result(
        vars_t,
        solutions=pts,
        method="sum_of_squares_family",
        complete=True,
        provenance=["special_family"],
        metadata=desc.metadata,
    )


def detect_pythag_fam(expr: sp.Expr, variables: Sequence[sp.Symbol]) -> SpecialFamDesc | None:
    variables = tuple(variables)
    if len(variables) != 3:
        return None
    eqs, others = _split_eq(expr)
    if len(eqs) != 1 or others:
        return None
    x, y, z = variables
    diff = sp.expand(eqs[0].lhs - eqs[0].rhs)
    if (
        sp.expand(diff - (x**2 + y**2 - z**2)) == 0
        or sp.expand(diff - (-(x**2) - y**2 + z**2)) == 0
    ):
        return SpecialFamDesc("pythagorean_triples", variables, {})
    return None


def pythag_triple_form(variables: Sequence[sp.Symbol]) -> sp.Expr:
    x, y, z = tuple(variables)
    m, n, k = sp.symbols("m n k", integer=True)
    f1 = sp.And(sp.Eq(x, k * (m**2 - n**2)), sp.Eq(y, k * (2 * m * n)), sp.Eq(z, k * (m**2 + n**2)))
    f2 = sp.And(sp.Eq(y, k * (m**2 - n**2)), sp.Eq(x, k * (2 * m * n)), sp.Eq(z, k * (m**2 + n**2)))
    return sp.Or(f1, f2)


def solve_pythag_triples_fam(
    expr: sp.Expr, variables: Sequence[sp.Symbol], *, parameter_bound: int = 10
):
    desc = detect_pythag_fam(expr, variables)
    if desc is None:
        return None
    x, y, z = desc.variables
    pts = []
    for m in range(-parameter_bound, parameter_bound + 1):
        for n in range(-parameter_bound, parameter_bound + 1):
            for k in range(-parameter_bound, parameter_bound + 1):
                if k == 0:
                    continue
                a = k * (m * m - n * n)
                b = k * (2 * m * n)
                c = k * (m * m + n * n)
                pts.extend([(a, b, c), (b, a, c)])
    return canon_int_result(
        (x, y, z),
        formula=pythag_triple_form((x, y, z)),
        solutions=pts,
        method="pythagorean_triples_family",
        complete=False,
        provenance=["special_family"],
        metadata={"parameter_bound": parameter_bound},
    )


def detect_diag_fam(expr: sp.Expr, variables: Sequence[sp.Symbol]) -> SpecialFamDesc | None:
    variables = tuple(variables)
    eqs, others = _split_eq(expr)
    if len(eqs) != 1 or others:
        return None
    diff = sp.expand(eqs[0].lhs - eqs[0].rhs)
    poly = sp.Poly(diff, *variables)
    nonconst_terms = [(mon, coeff) for mon, coeff in poly.terms() if sum(mon) > 0 and coeff != 0]
    if not nonconst_terms:
        return None
    exponents = set()
    coeffs = []
    for mon, coeff in nonconst_terms:
        if coeff not in (1, -1):
            return None
        nz = [e for e in mon if e != 0]
        if len(nz) != 1:
            return None
        exponents.add(nz[0])
        coeffs.append(int(coeff))
    if len(exponents) != 1:
        return None
    k = list(exponents)[0]
    target = -poly.coeff_monomial(1)
    if target.is_integer:
        return SpecialFamDesc(
            "diagonal_sum_of_powers",
            variables,
            {"power": int(k), "target": int(target), "coefficients": tuple(coeffs)},
        )
    return None


def solve_diag_fam(expr: sp.Expr, variables: Sequence[sp.Symbol], *, search_bound: int = 20):
    desc = detect_diag_fam(expr, variables)
    if desc is None:
        return None
    k = desc.metadata["power"]
    target = desc.metadata["target"]
    coeffs = desc.metadata["coefficients"]
    vars_t = tuple(variables)
    if k % 2 == 0 and target < -sum(0 for _ in coeffs):
        return canon_int_result(
            vars_t,
            formula=sp.false,
            solutions=[],
            method="diagonal_sum_of_powers_even_obstruction",
            complete=True,
            provenance=["special_family"],
            metadata=desc.metadata,
        )
    pts = []
    for vals in product(range(-search_bound, search_bound + 1), repeat=len(vars_t)):
        if sum(c * (v**k) for c, v in zip(coeffs, vals, strict=True)) == target:
            pts.append(vals)
    if not pts:
        return None
    return canon_int_result(
        vars_t,
        solutions=pts,
        method="diagonal_sum_of_powers_family",
        complete=False,
        provenance=["special_family"],
        metadata=desc.metadata,
    )


def detect_two_cubes_family(expr: sp.Expr, variables: Sequence[sp.Symbol]) -> SpecialFamDesc | None:
    variables = tuple(variables)
    if len(variables) != 2:
        return None
    desc = detect_diag_fam(expr, variables)
    if (
        desc is not None
        and desc.metadata.get("power") == 3
        and all(c in (1, -1) for c in desc.metadata.get("coefficients", ()))
    ):
        return SpecialFamDesc(
            "two_cubes",
            variables,
            {"target": desc.metadata["target"], "coefficients": desc.metadata["coefficients"]},
        )
    return None


def solve_two_cubes_family(
    expr: sp.Expr, variables: Sequence[sp.Symbol], *, search_bound: int = 100
):
    desc = detect_two_cubes_family(expr, variables)
    if desc is None:
        return None
    x, y = tuple(variables)
    target = desc.metadata["target"]
    c1, c2 = desc.metadata["coefficients"]
    pts = []
    cubes = {}
    for i in range(-search_bound, search_bound + 1):
        cubes.setdefault(i**3, []).append(i)
    for a in range(-search_bound, search_bound + 1):
        needed = target - c1 * (a**3)
        if c2 == -1:
            needed = -needed
        for b in cubes.get(needed, []):
            pts.append((a, b))
    return canon_int_result(
        (x, y),
        solutions=pts,
        method="two_cubes_family",
        complete=False,
        provenance=["special_family"],
        metadata=desc.metadata,
    )


def detect_pell_family(expr: sp.Expr, variables: Sequence[sp.Symbol]) -> SpecialFamDesc | None:
    variables = tuple(variables)
    if len(variables) != 2:
        return None
    eqs, others = _split_eq(expr)
    if len(eqs) != 1 or others:
        return None
    x, y = variables
    diff = sp.expand(eqs[0].lhs - eqs[0].rhs)
    # x^2 - D y^2 - N == 0
    poly = sp.Poly(diff, x, y)
    if poly.coeff_monomial(x**2) != 1:
        return None
    if poly.coeff_monomial(x * y) != 0 or poly.coeff_monomial(y**2) >= 0:
        return None
    D = -poly.coeff_monomial(y**2)
    N = -poly.coeff_monomial(1)
    if D.is_integer and N.is_integer and int(D) > 0 and int(D) != int(sp.isqrt(int(D))) ** 2:
        return SpecialFamDesc("pell_equation", variables, {"D": int(D), "N": int(N)})
    return None


def solve_pell_family(expr: sp.Expr, variables: Sequence[sp.Symbol], *, search_bound: int = 200):
    desc = detect_pell_family(expr, variables)
    if desc is None:
        return None
    x, y = tuple(variables)
    D = desc.metadata["D"]
    N = desc.metadata["N"]
    pts = []
    for b in range(-search_bound, search_bound + 1):
        rhs = D * b * b + N
        if rhs >= 0:
            a = int(sp.isqrt(rhs))
            if a * a == rhs:
                pts.extend([(a, b), (-a, b)])
    return canon_int_result(
        (x, y),
        solutions=pts,
        method="pell_family_search",
        complete=False,
        provenance=["special_family", "pell_family"],
        metadata=desc.metadata,
    )


def solve_binary_scan(expr: sp.Expr, variables: Sequence[sp.Symbol], *, x_bound: int = 200):
    variables = tuple(variables)
    if len(variables) != 2:
        return None
    eqs, others = _split_eq(expr)
    if len(eqs) != 1:
        return None
    x, y = variables
    diff = sp.expand(eqs[0].lhs - eqs[0].rhs)
    poly = sp.Poly(diff, x, y)
    if poly.total_degree() < 3 or others:
        return None
    pts = set()
    # Direct vertical scan
    for xv in range(-x_bound, x_bound + 1):
        univ = sp.Poly(sp.expand(diff.subs(x, xv)), y)
        try:
            roots = univ.all_roots()
        except Exception:
            continue
        for r in roots:
            sr = sp.simplify(r)
            if sr.is_integer is True:
                pts.add((xv, int(sr)))
    # Rational-slope scan for homogeneous / near-homogeneous behaviour
    t = sp.Symbol("_t")
    try:
        slope_poly = sp.Poly(sp.expand(diff.subs({x: t, y: 1})), t)
        rat_roots = []
        for rr in slope_poly.ground_roots().keys():
            srr = sp.simplify(rr)
            if srr.is_rational:
                rat_roots.append(sp.Rational(srr))
        for rr in rat_roots:
            p, q = int(rr.p), int(rr.q)
            for k in range(-x_bound, x_bound + 1):
                a, b = p * k, q * k
                if (
                    abs(a) <= x_bound
                    and abs(b) <= x_bound
                    and sp.expand(diff.subs({x: a, y: b})) == 0
                ):
                    pts.add((a, b))
    except Exception:
        pass
    if not pts:
        return None
    return canon_int_result(
        (x, y),
        solutions=sorted(pts),
        method="binary_form_asymptotic_scan",
        complete=False,
        provenance=["special_family", "asymptotic_scan"],
        metadata={"x_bound": x_bound},
    )


def solve_int_fams(expr: sp.Expr, variables: Sequence[sp.Symbol]):
    # First give the dedicated Thue-family layer a chance.
    try:
        from .thue import solve_binary_lll

        thue = solve_binary_lll(expr, variables, search_bound=200)
    except Exception:
        thue = None
    if thue is None:
        try:
            thue = solve_binary_bounded(expr, variables, search_bound=200)
        except Exception:
            thue = None
    if thue is not None:
        return thue

    for solver in (
        lambda e, vs: solve_sum_of_squares_fam(e, vs),
        lambda e, vs: solve_pythag_triples_fam(e, vs),
        lambda e, vs: solve_diag_fam(e, vs),
        lambda e, vs: solve_two_cubes_family(e, vs),
        lambda e, vs: solve_pell_family(e, vs),
        lambda e, vs: solve_binary_scan(e, vs),
    ):
        try:
            out = solver(expr, variables)
        except Exception:
            out = None
        if out is not None:
            return out
    return None


__all__ = [
    "SpecialFamDesc",
    "detect_sum_fam2",
    "solve_sum_of_squares_fam",
    "detect_pythag_fam",
    "pythag_triple_form",
    "solve_pythag_triples_fam",
    "detect_diag_fam",
    "solve_diag_fam",
    "detect_two_cubes_family",
    "solve_two_cubes_family",
    "detect_pell_family",
    "solve_pell_family",
    "solve_binary_scan",
    "solve_int_fams",
]
