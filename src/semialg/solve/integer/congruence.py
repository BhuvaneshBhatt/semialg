from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass, field
from itertools import product
from math import gcd

import sympy as sp
from sympy import Eq, Ne
from sympy.matrices.normalforms import smith_normal_decomp


@dataclass
class ModularSolveResult:
    modulus: int
    variables: tuple[sp.Symbol, ...]
    points: list[tuple[int, ...]] = field(default_factory=list)
    formula: sp.Expr = sp.false
    method: str = "enumeration"
    complete: bool = True
    metadata: dict = field(default_factory=dict)


def _coerce_modulus(modulus: int) -> int:
    modulus = int(modulus)
    if modulus <= 0:
        raise ValueError("Modulus must be a positive integer")
    return modulus


def _is_prime_modulus(modulus: int) -> bool:
    return bool(sp.ntheory.primetest.isprime(_coerce_modulus(modulus)))


def _factor_prime_powers(modulus: int) -> list[int]:
    modulus = _coerce_modulus(modulus)
    fac = sp.factorint(modulus)
    return [int(p**e) for p, e in sorted(fac.items())]


def _crt_pair(a1: int, m1: int, a2: int, m2: int) -> int:
    g = gcd(m1, m2)
    if (a1 - a2) % g != 0:
        raise ValueError("Incompatible CRT residues")
    lcm = m1 // g * m2
    if m1 == 1:
        return a2 % m2
    if m2 == 1:
        return a1 % m1
    m1r = m1 // g
    m2r = m2 // g
    inv = pow(m1r, -1, m2r)
    t = ((a2 - a1) // g * inv) % m2r
    return (a1 + m1 * t) % lcm


def combine_mod_crt(
    point_sets: Sequence[list[tuple[int, ...]]],
    variables: Sequence[sp.Symbol],
    moduli: Sequence[int],
) -> list[tuple[int, ...]]:
    variables = tuple(variables)
    moduli = tuple(map(_coerce_modulus, moduli))
    if len(point_sets) != len(moduli):
        raise ValueError("point_sets and moduli must have the same length")
    if not point_sets:
        return [tuple()]
    combined = [tuple([0] * len(variables))]
    combined_modulus = 1
    for points, modulus in zip(point_sets, moduli, strict=True):
        new_combined: list[tuple[int, ...]] = []
        for prefix in combined:
            for pt in points:
                merged = []
                ok = True
                for a, b in zip(prefix, pt, strict=True):
                    try:
                        merged.append(_crt_pair(int(a), combined_modulus, int(b), modulus))
                    except ValueError:
                        ok = False
                        break
                if ok:
                    new_combined.append(tuple(merged))
        combined = sorted(set(new_combined))
        combined_modulus *= modulus
    return combined


def _rel_to_mod_expr(rel: sp.Expr, modulus: int) -> sp.Expr:
    modulus = _coerce_modulus(modulus)
    if isinstance(rel, Eq):
        return sp.Eq(sp.Mod(sp.expand(rel.lhs - rel.rhs), modulus), 0)
    if isinstance(rel, Ne):
        return sp.Ne(sp.Mod(sp.expand(rel.lhs - rel.rhs), modulus), 0)
    if rel in (sp.true, sp.false):
        return rel
    return sp.Eq(sp.Mod(sp.expand(rel), modulus), 0)


def norm_mod_form(expr: sp.Expr, modulus: int) -> sp.Expr:
    modulus = _coerce_modulus(modulus)
    if expr in (sp.true, sp.false):
        return expr
    if isinstance(expr, (Eq, Ne)):
        return _rel_to_mod_expr(expr, modulus)
    if isinstance(expr, sp.Not):
        return sp.Not(norm_mod_form(expr.args[0], modulus))
    if isinstance(expr, sp.And):
        return sp.And(*[norm_mod_form(arg, modulus) for arg in expr.args])
    if isinstance(expr, sp.Or):
        return sp.Or(*[norm_mod_form(arg, modulus) for arg in expr.args])
    return _rel_to_mod_expr(expr, modulus)


def _eval_modular_truth(expr: sp.Expr, assignment: dict[sp.Symbol, int], modulus: int) -> bool:
    val = norm_mod_form(expr, modulus).subs(assignment)
    try:
        if val in (True, False):
            return bool(val)
        return bool(sp.simplify(val))
    except Exception:
        return False


def _points_to_formula(points: list[tuple[int, ...]], variables: Sequence[sp.Symbol]) -> sp.Expr:
    if not points:
        return sp.false
    pieces = []
    for pt in points:
        pieces.append(sp.And(*[sp.Eq(v, int(a)) for v, a in zip(variables, pt, strict=True)]))
    return sp.Or(*pieces) if len(pieces) > 1 else pieces[0]


def _enumerate_points(
    variables: Sequence[sp.Symbol], modulus: int, max_points: int
) -> list[tuple[int, ...]]:
    total = modulus ** len(tuple(variables))
    if total > max_points:
        raise RuntimeError(f"Too many modular points to enumerate exactly: {total}")
    return [tuple(p) for p in product(range(modulus), repeat=len(tuple(variables)))]


def _split_eqs_neqs(expr: sp.Expr) -> tuple[list[sp.Expr], list[sp.Expr]]:
    if expr in (sp.true, sp.false):
        return [], []
    atoms = list(expr.args) if isinstance(expr, sp.And) else [expr]
    eqs = [a for a in atoms if isinstance(a, sp.Equality)]
    neqs = [a for a in atoms if isinstance(a, sp.Unequality)]
    others = [a for a in atoms if not isinstance(a, (sp.Equality, sp.Unequality))]
    return eqs + others, neqs


def simp_mod_ineqs(inequations: Sequence[sp.Expr], modulus: int) -> list[sp.Expr]:
    """Simplify modular inequations with easy gcd/unit/content rules."""
    modulus = _coerce_modulus(modulus)
    simplified: list[sp.Expr] = []
    merged_by_tail: dict[str, tuple[int, sp.Expr]] = {}
    for rel in inequations:
        if rel in (sp.true, sp.false):
            simplified.append(rel)
            continue
        if not isinstance(rel, Ne):
            simplified.append(rel)
            continue
        diff = sp.expand(rel.lhs - rel.rhs)
        if diff.is_Integer:
            simplified.append(sp.true if int(diff) % modulus != 0 else sp.false)
            continue
        poly = sp.Poly(diff)
        if poly is None:
            simplified.append(rel)
            continue
        cont = int(poly.content())
        primitive = sp.expand(diff // cont) if cont != 0 else diff
        key = sp.srepr(primitive)
        if cont != 0 and gcd(abs(cont), modulus) == 1:
            inv = pow(cont % modulus, -1, modulus)
            simplified_expr = sp.Ne(sp.Mod(sp.expand(diff * inv), modulus), 0)
            merged_by_tail[key] = (1, simplified_expr)
        else:
            prev = merged_by_tail.get(key)
            if prev is None:
                merged_by_tail[key] = (abs(cont), rel)
            else:
                merged = gcd(prev[0], abs(cont))
                merged_by_tail[key] = (merged, sp.Ne(sp.expand(merged * primitive), 0))
    for _, expr in merged_by_tail.values():
        simplified.append(expr)
    out = [s for s in simplified if s is not sp.true]
    if any(s is sp.false for s in out):
        return [sp.false]
    return out


def _supports_direct_solver(equalities: Sequence[sp.Expr], variables: Sequence[sp.Symbol]) -> bool:
    for rel in equalities:
        if not isinstance(rel, Eq):
            return False
        poly = sp.Poly(sp.expand(rel.lhs - rel.rhs), *variables)
        if poly.total_degree() > 1:
            return False
    return True


def _reduced_eqn_groebner(
    equalities: Sequence[sp.Expr], variables: Sequence[sp.Symbol], modulus: int
):
    polys = [sp.expand(eq.lhs - eq.rhs) for eq in equalities if isinstance(eq, Eq)]
    if not polys:
        return []
    try:
        gb = sp.groebner(polys, *reversed(tuple(variables)), modulus=modulus, order="lex")
        return [sp.expand(g.as_expr()) for g in gb.polys if sp.expand(g.as_expr()) != 0]
    except Exception:
        return polys


def _linear_system_matrix(
    equalities: Sequence[sp.Expr], variables: Sequence[sp.Symbol], modulus: int
):
    rows = []
    rhs = []
    for rel in equalities:
        if not isinstance(rel, Eq):
            return None
        poly = sp.Poly(sp.expand(rel.lhs - rel.rhs), *variables)
        if poly.total_degree() > 1:
            return None
        row = []
        for v in variables:
            row.append(int(sp.Mod(poly.coeff_monomial(v), modulus)))
        const = int(sp.Mod(poly.coeff_monomial(1), modulus))
        rows.append(row)
        rhs.append(int((-const) % modulus))
    return rows, rhs


def _row_reduce_prime_mod(matrix: list[list[int]], rhs: list[int], modulus: int):
    n_rows = len(matrix)
    n_cols = len(matrix[0]) if matrix else 0
    A = [row[:] + [rhs_i] for row, rhs_i in zip(matrix, rhs, strict=True)]
    pivots = []
    r = 0
    for c in range(n_cols):
        pivot = None
        for i in range(r, n_rows):
            if A[i][c] % modulus != 0:
                pivot = i
                break
        if pivot is None:
            continue
        A[r], A[pivot] = A[pivot], A[r]
        inv = pow(A[r][c] % modulus, -1, modulus)
        A[r] = [(inv * v) % modulus for v in A[r]]
        for i in range(n_rows):
            if i == r:
                continue
            factor = A[i][c] % modulus
            if factor:
                A[i] = [(a - factor * b) % modulus for a, b in zip(A[i], A[r], strict=True)]
        pivots.append((r, c))
        r += 1
        if r == n_rows:
            break
    for i in range(n_rows):
        if all(A[i][j] % modulus == 0 for j in range(n_cols)) and A[i][-1] % modulus != 0:
            return None, None, None
    return A, pivots, [c for c in range(n_cols) if c not in [pc for _, pc in pivots]]


def enum_lin_sols_prime(
    equalities: Sequence[sp.Expr], variables: Sequence[sp.Symbol], modulus: int, max_points: int
):
    data = _linear_system_matrix(equalities, variables, modulus)
    if data is None:
        return None
    matrix, rhs = data
    rref, pivots, free_cols = _row_reduce_prime_mod(matrix, rhs, modulus)
    if rref is None:
        return []
    if modulus ** len(free_cols) > max_points:
        raise RuntimeError("Too many modular free-variable assignments for exact linear solving")
    points = []
    for free_vals in product(range(modulus), repeat=len(free_cols)):
        sol = [0] * len(variables)
        for c, val in zip(free_cols, free_vals, strict=True):
            sol[c] = val
        for row_idx, col_idx in pivots:
            value = rref[row_idx][-1]
            for c in free_cols:
                value = (value - rref[row_idx][c] * sol[c]) % modulus
            sol[col_idx] = value % modulus
        points.append(tuple(sol))
    return points


def solve_lin_sys_mod_comp(
    equalities: Sequence[sp.Expr], variables: Sequence[sp.Symbol], modulus: int, max_points: int
):
    """Solve A x = b mod m using Smith normal form over the integers."""
    data = _linear_system_matrix(equalities, variables, modulus)
    if data is None:
        return None
    matrix, rhs = data
    A = sp.Matrix(matrix)
    b = sp.Matrix(rhs)
    D, U, V = smith_normal_decomp(A, domain=sp.ZZ)
    b2 = U * b
    n = A.shape[1]

    # Consistency and solution of diagonal congruences d_i z_i ≡ c_i mod m
    value_sets: list[list[int]] = []
    rank = min(D.rows, D.cols)
    for i in range(rank):
        d = int(D[i, i])
        c = int(b2[i, 0])
        if d == 0:
            if c % modulus != 0:
                return []
            value_sets.append(list(range(modulus)))
            continue
        g = gcd(d, modulus)
        if c % g != 0:
            return []
        m_red = modulus // g
        d_red = d // g
        c_red = c // g
        inv = pow(d_red % m_red, -1, m_red)
        base = (c_red * inv) % m_red
        vals = [(base + k * m_red) % modulus for k in range(g)]
        value_sets.append(sorted(set(vals)))
    for i in range(rank, A.rows):
        if int(b2[i, 0]) % modulus != 0:
            return []
    for _i in range(rank, n):
        value_sets.append(list(range(modulus)))
    total = 1
    for vals in value_sets:
        total *= len(vals)
    if total > max_points:
        raise RuntimeError("Too many composite-modulus linear solutions for exact enumeration")
    points = []
    for z in product(*value_sets):
        xvec = V * sp.Matrix(list(z))
        points.append(tuple(int(v % modulus) for v in xvec))
    return sorted(set(points))


def solve_mod_lin_sys(
    equalities: Sequence[sp.Expr],
    variables: Sequence[sp.Symbol],
    modulus: int,
    *,
    max_points: int = 10000,
) -> ModularSolveResult:
    modulus = _coerce_modulus(modulus)
    variables = tuple(variables)
    method = "linear_enumeration"
    points = None

    if modulus > 1 and not _is_prime_modulus(modulus):
        factors = _factor_prime_powers(modulus)
        if len(factors) > 1:
            component_results = [
                solve_mod_lin_sys(equalities, variables, q, max_points=max_points) for q in factors
            ]
            point_sets = [r.points for r in component_results]
            try:
                points = combine_mod_crt(point_sets, variables, factors)
                method = "crt_decomposed_linear_solver"
            except Exception:
                points = None
        if points is None:
            try:
                points = solve_lin_sys_mod_comp(equalities, variables, modulus, max_points)
                method = "smith_normal_form_mod_composite"
            except Exception:
                points = None
    elif _is_prime_modulus(modulus):
        try:
            points = enum_lin_sols_prime(equalities, variables, modulus, max_points)
            method = "linear_row_reduction_mod_prime"
        except Exception:
            points = None

    if points is None:
        points = []
        for pt in _enumerate_points(variables, modulus, max_points):
            subst = dict(zip(variables, pt, strict=True))
            if all(_eval_modular_truth(eq, subst, modulus) for eq in equalities):
                points.append(pt)
    return ModularSolveResult(
        modulus=modulus,
        variables=variables,
        points=points,
        formula=_points_to_formula(points, variables),
        method=method,
        complete=True,
        metadata={
            "equation_count": len(tuple(equalities)),
            "modulus_factors": _factor_prime_powers(modulus),
        },
    )


def _cand_values_poly(poly_expr: sp.Expr, variable: sp.Symbol, modulus: int) -> list[int]:
    return [
        a
        for a in range(modulus)
        if int(sp.Mod(sp.expand(poly_expr).subs(variable, a), modulus)) == 0
    ]


def score_branch_var(
    polys: Sequence[sp.Expr], variables: Sequence[sp.Symbol], modulus: int
) -> tuple[sp.Symbol, list[int], list[sp.Expr]]:
    best = None
    for var in variables:
        free_of_others = [p for p in polys if p.free_symbols.issubset({var})]
        candidate_vals = None
        if free_of_others:
            for p in free_of_others:
                vals = set(_cand_values_poly(p, var, modulus))
                candidate_vals = (
                    vals if candidate_vals is None else candidate_vals.intersection(vals)
                )
            candidate_vals = sorted(candidate_vals or [])
        else:
            candidate_vals = list(range(modulus))
        score = len(candidate_vals)
        tup = (score, var.sort_key(), var, candidate_vals, free_of_others)
        if best is None or tup < best:
            best = tup
    assert best is not None
    return best[2], best[3], best[4]


def rec_solve_basis(
    groebner_polynomials: Sequence[sp.Expr],
    variables: Sequence[sp.Symbol],
    modulus: int,
    *,
    max_points: int = 10000,
    inequations: Sequence[sp.Expr] | None = None,
    partial_assignment: dict[sp.Symbol, int] | None = None,
) -> list[tuple[int, ...]]:
    modulus = _coerce_modulus(modulus)
    variables = tuple(variables)
    inequations = tuple(inequations or ())
    assignment = dict(partial_assignment or {})
    polys = [
        sp.expand(p.subs(assignment))
        for p in groebner_polynomials
        if sp.expand(p.subs(assignment)) != 0
    ]
    neqs = [sp.expand(p.subs(assignment)) for p in inequations]

    # Early contradiction pruning
    for p in polys:
        if not p.free_symbols and int(sp.Mod(p, modulus)) != 0:
            return []
    for q in neqs:
        if not q.free_symbols and int(sp.Mod(q, modulus)) == 0:
            return []

    active_vars = tuple(v for v in variables if v not in assignment)
    if not active_vars:
        ordered = tuple(int(assignment[v] % modulus) for v in variables)
        return [ordered]

    if not polys:
        pts = []
        for pt in _enumerate_points(active_vars, modulus, max_points):
            subst = dict(zip(active_vars, pt, strict=True))
            if all(int(sp.Mod(q.subs(subst), modulus)) != 0 for q in neqs):
                full = dict(assignment)
                full.update(subst)
                pts.append(tuple(int(full[v] % modulus) for v in variables))
        return pts

    target, allowed, direct_univariates = score_branch_var(polys, active_vars, modulus)
    # Advanced pruning: use inequations depending only on target to filter residues
    target_ineqs = [q for q in neqs if q.free_symbols.issubset({target})]
    if target_ineqs:
        filtered = []
        for a in allowed:
            ok = True
            for q in target_ineqs:
                if int(sp.Mod(q.subs(target, a), modulus)) == 0:
                    ok = False
                    break
            if ok:
                filtered.append(a)
        allowed = filtered

    pts: list[tuple[int, ...]] = []
    budget_per_branch = max(1, max_points // max(1, len(allowed)))
    for val in allowed:
        new_assignment = dict(assignment)
        new_assignment[target] = val
        child = rec_solve_basis(
            polys,
            variables,
            modulus,
            max_points=budget_per_branch,
            inequations=neqs,
            partial_assignment=new_assignment,
        )
        pts.extend(child)
        if len(pts) > max_points:
            raise RuntimeError("Too many recursively constructed modular solutions")
    return sorted(set(pts))


def solve_mod_poly_sys(
    equalities: Sequence[sp.Expr],
    variables: Sequence[sp.Symbol],
    modulus: int,
    *,
    max_points: int = 10000,
) -> ModularSolveResult:
    modulus = _coerce_modulus(modulus)
    variables = tuple(variables)
    method = "polynomial_enumeration"
    points = None

    factors = _factor_prime_powers(modulus)
    if len(factors) > 1:
        try:
            parts = [
                solve_mod_poly_sys(equalities, variables, q, max_points=max_points) for q in factors
            ]
            points = combine_mod_crt([p.points for p in parts], variables, factors)
            method = "crt_decomposed_polynomial_solver"
        except Exception:
            points = None

    reduced = _reduced_eqn_groebner(equalities, variables, modulus)
    if points is None:
        try:
            points = rec_solve_basis(reduced, variables, modulus, max_points=max_points)
            method = "recursive_groebner_modular_solver"
        except Exception:
            points = None

    if points is None:
        eq_exprs = [sp.Eq(poly, 0) for poly in reduced]
        points = []
        for pt in _enumerate_points(variables, modulus, max_points):
            subst = dict(zip(variables, pt, strict=True))
            if all(_eval_modular_truth(eq, subst, modulus) for eq in eq_exprs):
                points.append(pt)
    return ModularSolveResult(
        modulus=modulus,
        variables=variables,
        points=points,
        formula=_points_to_formula(points, variables),
        method=method,
        complete=True,
        metadata={
            "groebner_polynomials": [sp.srepr(p) for p in reduced],
            "modulus_factors": factors,
        },
    )


def solve_quant_free_mod_sys(
    expr: sp.Expr, variables: Sequence[sp.Symbol], modulus: int, *, max_points: int = 10000
) -> ModularSolveResult:
    modulus = _coerce_modulus(modulus)
    variables = tuple(variables)
    normalized = norm_mod_form(expr, modulus)

    factors = _factor_prime_powers(modulus)
    if len(factors) > 1:
        try:
            parts = [
                solve_quant_free_mod_sys(normalized, variables, q, max_points=max_points)
                for q in factors
            ]
            points = combine_mod_crt([p.points for p in parts], variables, factors)
            return ModularSolveResult(
                modulus=modulus,
                variables=variables,
                points=points,
                formula=_points_to_formula(points, variables),
                method="crt_decomposed_quantifier_free_solver",
                complete=True,
                metadata={"component_moduli": factors},
            )
        except Exception:
            pass

    eqs, neqs = _split_eqs_neqs(normalized)
    neqs = simp_mod_ineqs(neqs, modulus)
    if any(n is sp.false for n in neqs):
        return ModularSolveResult(
            modulus, variables, [], sp.false, "infeasible_by_inequation_simplification", True, {}
        )

    if eqs and not neqs and _supports_direct_solver(eqs, variables):
        return solve_mod_lin_sys(eqs, variables, modulus, max_points=max_points)

    if eqs:
        reduced = _reduced_eqn_groebner([e for e in eqs if isinstance(e, Eq)], variables, modulus)
        try:
            points = rec_solve_basis(
                reduced,
                variables,
                modulus,
                max_points=max_points,
                inequations=[sp.expand(n.lhs - n.rhs) for n in neqs if isinstance(n, Ne)],
            )
            return ModularSolveResult(
                modulus=modulus,
                variables=variables,
                points=points,
                formula=_points_to_formula(points, variables),
                method="recursive_quantifier_free_modular_solver",
                complete=True,
                metadata={
                    "groebner_polynomials": [sp.srepr(p) for p in reduced],
                    "inequation_count": len(neqs),
                },
            )
        except Exception:
            pass

    points = []
    rebuilt = sp.And(*(list(eqs) + list(neqs))) if eqs or neqs else sp.true
    for pt in _enumerate_points(variables, modulus, max_points):
        subst = dict(zip(variables, pt, strict=True))
        if _eval_modular_truth(rebuilt, subst, modulus):
            points.append(pt)

    return ModularSolveResult(
        modulus=modulus,
        variables=variables,
        points=points,
        formula=_points_to_formula(points, variables),
        method="quantifier_free_enumeration",
        complete=True,
        metadata={"normalized_formula": sp.srepr(rebuilt)},
    )


def find_quant_free_mod_inst(
    expr: sp.Expr, variables: Sequence[sp.Symbol], modulus: int, *, max_points: int = 10000
):
    result = solve_quant_free_mod_sys(expr, variables, modulus, max_points=max_points)
    if not result.points:
        return None
    return dict(zip(result.variables, result.points[0], strict=True))


def eliminate_one_var(
    expr: sp.Expr,
    variables: Sequence[sp.Symbol],
    quantified_variable: sp.Symbol,
    modulus: int,
    *,
    max_points: int = 10000,
) -> sp.Expr:
    modulus = _coerce_modulus(modulus)
    variables = tuple(variables)
    remaining = tuple(v for v in variables if v != quantified_variable)
    if len(remaining) != len(variables) - 1:
        raise ValueError("Quantified variable must be one of the provided variables")

    normalized = norm_mod_form(expr, modulus)
    factors = _factor_prime_powers(modulus)
    if len(factors) > 1:
        part_formulas = [
            eliminate_one_var(normalized, variables, quantified_variable, q, max_points=max_points)
            for q in factors
        ]
        # Existential elimination over CRT factors corresponds to conjunction of factorwise feasibility.
        return sp.And(*part_formulas)

    eqs, neqs = _split_eqs_neqs(normalized)
    neqs = simp_mod_ineqs(neqs, modulus)

    if all((not isinstance(eq, Eq)) or quantified_variable not in eq.free_symbols for eq in eqs):
        surviving = []
        for pt in _enumerate_points(remaining, modulus, max_points):
            subst = dict(zip(remaining, pt, strict=True))
            if not all(_eval_modular_truth(eq, subst, modulus) for eq in eqs):
                continue
            forbidden = set()
            heuristic_ok = True
            for ne in neqs:
                if not isinstance(ne, Ne):
                    heuristic_ok = False
                    break
                diff = sp.expand((ne.lhs - ne.rhs).subs(subst))
                poly = sp.Poly(diff, quantified_variable)
                if poly.total_degree() <= 2:
                    forbidden.update(
                        a
                        for a in range(modulus)
                        if int(sp.Mod(poly.as_expr().subs(quantified_variable, a), modulus)) == 0
                    )
                else:
                    heuristic_ok = False
                    break
            if heuristic_ok:
                if len(forbidden) < modulus:
                    surviving.append(pt)
                continue
            # Fallback local search only when heuristic stalls.
            witnessed = False
            for q in range(modulus):
                subst_q = dict(subst)
                subst_q[quantified_variable] = q
                if _eval_modular_truth(normalized, subst_q, modulus):
                    witnessed = True
                    break
            if witnessed:
                surviving.append(pt)
        return _points_to_formula(surviving, remaining)

    # If qvar appears in equations, solve the quantified-free problem and project.
    solved = solve_quant_free_mod_sys(normalized, variables, modulus, max_points=max_points)
    projected = sorted(
        set(
            tuple(pt[i] for i, v in enumerate(variables) if v != quantified_variable)
            for pt in solved.points
        )
    )
    return _points_to_formula(projected, remaining)


def solve_quant_mod_sys(
    expr: sp.Expr,
    free_variables: Sequence[sp.Symbol],
    quantified_variables: Sequence[sp.Symbol],
    modulus: int,
    *,
    quantifier: str = "exists",
    max_points: int = 10000,
) -> ModularSolveResult:
    modulus = _coerce_modulus(modulus)
    free_variables = tuple(free_variables)
    quantified_variables = tuple(quantified_variables)
    if quantifier not in {"exists", "forall"}:
        raise ValueError("quantifier must be 'exists' or 'forall'")

    factors = _factor_prime_powers(modulus)
    if len(factors) > 1:
        try:
            parts = [
                solve_quant_mod_sys(
                    expr,
                    free_variables,
                    quantified_variables,
                    q,
                    quantifier=quantifier,
                    max_points=max_points,
                )
                for q in factors
            ]
            points = combine_mod_crt([p.points for p in parts], free_variables, factors)
            return ModularSolveResult(
                modulus=modulus,
                variables=free_variables,
                points=points,
                formula=_points_to_formula(points, free_variables),
                method=f"crt_recursive_{quantifier}_quantified_solver",
                complete=True,
                metadata={"component_moduli": factors},
            )
        except Exception:
            pass

    formula = norm_mod_form(expr, modulus)
    current_vars = free_variables + quantified_variables
    for qvar in reversed(quantified_variables):
        if quantifier == "exists":
            formula = eliminate_one_var(formula, current_vars, qvar, modulus, max_points=max_points)
        else:
            neg = sp.Not(formula)
            eliminated = eliminate_one_var(neg, current_vars, qvar, modulus, max_points=max_points)
            formula = sp.Not(eliminated)
        current_vars = tuple(v for v in current_vars if v != qvar)

    result = solve_quant_free_mod_sys(formula, free_variables, modulus, max_points=max_points)
    result.method = f"recursive_{quantifier}_modular_quantifier_elimination"
    result.metadata["quantified_variables"] = tuple(map(str, quantified_variables))
    result.metadata["eliminated_formula"] = sp.srepr(formula)
    return result


# Convenience wrappers for modular solving
def solve_modular_system(
    expr: sp.Expr, variables: Sequence[sp.Symbol], modulus: int, *, max_points: int = 10000
) -> ModularSolveResult:
    return solve_quant_free_mod_sys(expr, variables, modulus, max_points=max_points)


def find_modular_instance(
    expr: sp.Expr, variables: Sequence[sp.Symbol], modulus: int, *, max_points: int = 10000
):
    return find_quant_free_mod_inst(expr, variables, modulus, max_points=max_points)


__all__ = [
    "ModularSolveResult",
    "norm_mod_form",
    "combine_mod_crt",
    "simp_mod_ineqs",
    "solve_mod_lin_sys",
    "solve_mod_poly_sys",
    "rec_solve_basis",
    "solve_quant_free_mod_sys",
    "find_quant_free_mod_inst",
    "eliminate_one_var",
    "solve_quant_mod_sys",
    "solve_modular_system",
    "find_modular_instance",
]
