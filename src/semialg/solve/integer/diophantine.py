from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass, field
from itertools import product

import sympy as sp
from sympy import Eq

from .congruence import solve_quant_free_mod_sys
from .factorization import solve_int_recursion
from .groebner_recursion import rec_reduce_sys
from .linear_divisibility import detect_lin_reduction
from .linear_recursion import rec_reduce_int_lin_sys
from .special_families import solve_int_fams


@dataclass
class IntEqnSolveResult:
    variables: tuple[sp.Symbol, ...]
    solutions: list[tuple[sp.Expr, ...]] = field(default_factory=list)
    formula: sp.Expr = sp.false
    method: str = "unknown"
    complete: bool = False
    metadata: dict = field(default_factory=dict)


def _conjuncts(expr: sp.Expr) -> list[sp.Expr]:
    return list(expr.args) if isinstance(expr, sp.And) else [expr]


def sol_points_to_form(
    variables: Sequence[sp.Symbol], points: Sequence[Sequence[sp.Expr]]
) -> sp.Expr:
    pts = [tuple(p) for p in points]
    if not pts:
        return sp.false
    pieces = []
    for pt in pts:
        pieces.append(sp.And(*[sp.Eq(v, val) for v, val in zip(variables, pt, strict=True)]))
    return sp.Or(*pieces) if len(pieces) > 1 else pieces[0]


def _dedupe_points(points: Sequence[Sequence[sp.Expr]]) -> list[tuple[sp.Expr, ...]]:
    seen = []
    keys = set()
    for pt in points:
        t = tuple(sp.simplify(x) for x in pt)
        key = tuple(sp.srepr(x) for x in t)
        if key not in keys:
            keys.add(key)
            seen.append(t)
    return seen


def solve_int_fams2(expr: sp.Expr, variables: Sequence[sp.Symbol]) -> IntEqnSolveResult | None:
    result = solve_int_fams(expr, variables)
    if result is None:
        return None
    return IntEqnSolveResult(
        variables=result.variables,
        solutions=result.solutions,
        formula=result.formula,
        method=result.method,
        complete=result.complete,
        metadata=result.metadata,
    )


def solve_int_divis(
    expr: sp.Expr, variables: Sequence[sp.Symbol], *, max_depth: int = 8
) -> IntEqnSolveResult | None:
    result = rec_reduce_int_lin_sys(expr, variables, max_depth=max_depth)
    if result is None:
        return None
    return IntEqnSolveResult(
        variables=result.variables,
        solutions=result.solutions,
        formula=result.formula,
        method=result.method,
        complete=result.complete,
        metadata=result.metadata,
    )


def detect_int_lin_elim(expr: sp.Expr, variables: Sequence[sp.Symbol]):
    return detect_lin_reduction(expr, variables)


def reduce_int_divis(expr: sp.Expr, variables: Sequence[sp.Symbol]) -> IntEqnSolveResult | None:
    reduction = detect_int_lin_elim(expr, variables)
    if reduction is None:
        return None
    reduced_vars = tuple(v for v in variables if v != reduction.solved_variable)
    if len(reduced_vars) == 1:
        var = reduced_vars[0]
        try:
            solset = sp.solveset(reduction.reduced_formula, var, domain=sp.S.Integers)
            points = []
            if isinstance(solset, sp.FiniteSet):
                for val in solset:
                    subst = {var: val}
                    solved_val = sp.simplify(reduction.replacement.subs(subst))
                    if solved_val.is_integer is True:
                        pt_map = {var: val, reduction.solved_variable: solved_val}
                        points.append(tuple(pt_map[v] for v in variables))
                points = _dedupe_points(points)
                return IntEqnSolveResult(
                    variables=tuple(variables),
                    solutions=points,
                    formula=sol_points_to_form(variables, points),
                    method="linear_divisibility_reduction",
                    complete=True,
                    metadata={"reduction": reduction},
                )
        except Exception:
            pass
    return IntEqnSolveResult(
        variables=tuple(variables),
        solutions=[],
        formula=sp.And(
            sp.Eq(reduction.solved_variable, reduction.replacement), reduction.reduced_formula
        ),
        method="linear_divisibility_reduction_symbolic",
        complete=False,
        metadata={"reduction": reduction},
    )


def solve_int_branches(expr: sp.Expr, variables: Sequence[sp.Symbol]) -> IntEqnSolveResult | None:
    result = solve_int_recursion(expr, variables)
    if result is None:
        return None
    return IntEqnSolveResult(
        variables=result.variables,
        solutions=result.solutions,
        formula=result.formula,
        method=result.method,
        complete=result.complete,
        metadata=result.metadata,
    )


def solve_int_sys_via_factor(
    expr: sp.Expr, variables: Sequence[sp.Symbol]
) -> IntEqnSolveResult | None:
    variables = tuple(variables)
    atoms = _conjuncts(expr)
    eqs = [a for a in atoms if isinstance(a, Eq)]
    others = [a for a in atoms if not isinstance(a, Eq)]
    if len(eqs) != 1:
        return None
    diff = sp.expand(eqs[0].lhs - eqs[0].rhs)
    coeff, factors = sp.factor_list(diff)
    nonconstant = [f for f, e in factors for _ in range(e)]
    if len(nonconstant) <= 1:
        return None
    branch_points = []
    branch_formulas = []
    for factor in set([f for f, _e in factors]):
        subexpr = sp.And(*(others + [sp.Eq(factor, 0)]))
        sub = solve_int_methods(subexpr, variables)
        if sub is not None:
            if sub.complete and sub.solutions:
                branch_points.extend(sub.solutions)
            branch_formulas.append(sub.formula)
        else:
            branch_formulas.append(subexpr)
    branch_points = _dedupe_points(branch_points)
    formula = (
        sol_points_to_form(variables, branch_points) if branch_points else sp.Or(*branch_formulas)
    )
    return IntEqnSolveResult(
        variables=variables,
        solutions=branch_points,
        formula=formula,
        method="factorization_branching",
        complete=bool(branch_points)
        and all(isinstance(f, IntEqnSolveResult) and f.complete for f in []),
        metadata={"factorized_equation": diff},
    )


def solve_int_branch(
    expr: sp.Expr,
    variables: Sequence[sp.Symbol],
    *,
    max_depth: int = 8,
    max_branch_points: int = 64,
) -> IntEqnSolveResult | None:
    result = rec_reduce_sys(
        expr, variables, max_depth=max_depth, max_branch_points=max_branch_points
    )
    if result is None:
        return None
    return IntEqnSolveResult(
        variables=result.variables,
        solutions=result.solutions,
        formula=result.formula,
        method=result.method,
        complete=result.complete,
        metadata=result.metadata,
    )


def int_roots_of_univar_poly(poly_expr: sp.Expr, var: sp.Symbol) -> list[sp.Expr]:
    try:
        roots = sp.solveset(sp.Eq(poly_expr, 0), var, domain=sp.S.Integers)
        if isinstance(roots, sp.FiniteSet):
            return list(sorted(roots, key=sp.default_sort_key))
    except Exception:
        pass
    try:
        roots = sp.Poly(poly_expr, var).all_roots()
        out = []
        for r in roots:
            sr = sp.simplify(r)
            if sr.is_integer is True:
                out.append(sr)
        return list(sorted(set(out), key=sp.default_sort_key))
    except Exception:
        return []


def solve_int_recursion2(
    expr: sp.Expr, variables: Sequence[sp.Symbol], *, max_branch_points: int = 200
) -> IntEqnSolveResult | None:
    variables = tuple(variables)
    atoms = _conjuncts(expr)
    eqs = [a for a in atoms if isinstance(a, Eq)]
    if not eqs:
        return None
    polys = [sp.expand(eq.lhs - eq.rhs) for eq in eqs]
    try:
        gb = sp.groebner(polys, *reversed(variables), order="lex")
        basis = [sp.expand(p.as_expr()) for p in gb.polys if sp.expand(p.as_expr()) != 0]
    except Exception:
        basis = polys

    # Seek a univariate polynomial in the earliest remaining variable.
    chosen_var = None
    chosen_poly = None
    for var in variables:
        candidates = [p for p in basis if p.free_symbols.issubset({var})]
        if candidates:
            chosen_var = var
            chosen_poly = candidates[0]
            break
    if chosen_var is None:
        return None

    roots = int_roots_of_univar_poly(chosen_poly, chosen_var)
    if not roots:
        return IntEqnSolveResult(
            variables=variables,
            solutions=[],
            formula=sp.false,
            method="groebner_recursive_univariate_root_pruning",
            complete=True,
            metadata={"basis": basis},
        )
    if len(roots) > max_branch_points:
        return None

    points = []
    branch_formulas = []
    remaining = tuple(v for v in variables if v != chosen_var)
    for root in roots:
        substituted_atoms = [
            sp.simplify(a.subs(chosen_var, root)) for a in atoms if a != sp.Eq(chosen_poly, 0)
        ]
        substituted_expr = sp.And(*substituted_atoms) if substituted_atoms else sp.true
        if not remaining:
            truth = sp.simplify(substituted_expr)
            if truth is not sp.false:
                points.append(tuple(root for _ in [chosen_var]))
            continue
        sub = solve_int_methods(substituted_expr, remaining)
        if sub is not None and sub.complete and sub.solutions:
            for tail in sub.solutions:
                mapping = {
                    chosen_var: root,
                    **{v: val for v, val in zip(remaining, tail, strict=True)},
                }
                points.append(tuple(mapping[v] for v in variables))
        else:
            branch_formulas.append(sp.And(sp.Eq(chosen_var, root), substituted_expr))

    points = _dedupe_points(points)
    formula = (
        sol_points_to_form(variables, points)
        if points
        else (sp.Or(*branch_formulas) if branch_formulas else sp.false)
    )
    return IntEqnSolveResult(
        variables=variables,
        solutions=points,
        formula=formula,
        method="groebner_recursive_integer_solver",
        complete=bool(points),
        metadata={"basis": basis, "chosen_variable": chosen_var, "roots": roots},
    )


def mod_res_cands(
    expr: sp.Expr, variables: Sequence[sp.Symbol], moduli: Sequence[int]
) -> list[tuple[int, ...]] | None:
    residue_sets = []
    for m in moduli:
        try:
            res = solve_quant_free_mod_sys(expr, variables, m, max_points=5000)
            if not res.points:
                return []
            residue_sets.append((m, res.points))
        except Exception:
            return None
    if not residue_sets:
        return None
    # Incrementally CRT-combine residue classes.
    combined = [tuple([0] * len(tuple(variables)))]
    from .congruence import combine_mod_crt

    try:
        combined = combine_mod_crt(
            [pts for _m, pts in residue_sets], variables, [m for m, _pts in residue_sets]
        )
        return combined
    except Exception:
        return None


def solve_int_pruning(
    expr: sp.Expr,
    variables: Sequence[sp.Symbol],
    *,
    search_radius: int = 10,
    moduli: Sequence[int] = (2, 3, 5, 7),
) -> IntEqnSolveResult | None:
    variables = tuple(variables)
    residues = mod_res_cands(expr, variables, moduli)
    if residues == []:
        return IntEqnSolveResult(
            variables=variables,
            solutions=[],
            formula=sp.false,
            method="modular_pruning_detected_inconsistency",
            complete=True,
            metadata={"moduli": tuple(moduli)},
        )
    if residues is None:
        return None

    combined_modulus = 1
    for m in moduli:
        combined_modulus *= m
    # Bounded lattice search around residue classes.
    per_var_span = max(1, 2 * search_radius + 1)
    if (per_var_span ** len(variables)) * max(1, len(residues)) > 50000:
        return None

    points = []
    for res in residues:
        for offsets in product(range(-search_radius, search_radius + 1), repeat=len(variables)):
            pt = tuple(int(r + combined_modulus * k) for r, k in zip(res, offsets, strict=True))
            subst = dict(zip(variables, pt, strict=True))
            truth = sp.simplify(expr.subs(subst))
            if truth is sp.true or truth is True:
                points.append(pt)
    points = _dedupe_points(points)
    return IntEqnSolveResult(
        variables=variables,
        solutions=points,
        formula=sol_points_to_form(variables, points),
        method="modular_pruning_with_bounded_lift_search",
        complete=bool(points),
        metadata={"moduli": tuple(moduli), "combined_modulus": combined_modulus},
    )


def detect_sum_eqn(expr: sp.Expr, variables: Sequence[sp.Symbol]):
    variables = tuple(variables)
    if len(variables) != 2:
        return None
    atoms = _conjuncts(expr)
    eqs = [a for a in atoms if isinstance(a, Eq)]
    others = [a for a in atoms if not isinstance(a, Eq)]
    if len(eqs) != 1 or others:
        return None
    x, y = variables
    diff = sp.expand(eqs[0].lhs - eqs[0].rhs)
    poly = sp.Poly(diff, x, y)
    if poly is None:
        return None
    if (
        poly.coeff_monomial(x**2) == 1
        and poly.coeff_monomial(y**2) == 1
        and poly.coeff_monomial(x * y) == 0
    ):
        c = -poly.coeff_monomial(1)
        if c.is_integer:
            return int(c)
    return None


def solve_sum_of_two_squares(
    expr: sp.Expr, variables: Sequence[sp.Symbol], *, search_limit: int = 1000
) -> IntEqnSolveResult | None:
    n = detect_sum_eqn(expr, variables)
    if n is None or n < 0:
        return None
    x, y = tuple(variables)
    points = []
    bound = int(sp.floor(sp.sqrt(n)))
    if bound > search_limit:
        return None
    for a in range(-bound, bound + 1):
        b2 = n - a * a
        if b2 < 0:
            continue
        b = int(sp.isqrt(b2))
        if b * b == b2:
            points.extend([(a, b), (a, -b)])
    points = _dedupe_points(points)
    return IntEqnSolveResult(
        variables=(x, y),
        solutions=points,
        formula=sol_points_to_form((x, y), points),
        method="sum_of_two_squares_special_solver",
        complete=True,
        metadata={"target_norm": n},
    )


def detect_binary_homog_eqn(expr: sp.Expr, variables: Sequence[sp.Symbol]):
    variables = tuple(variables)
    if len(variables) != 2:
        return None
    atoms = _conjuncts(expr)
    eqs = [a for a in atoms if isinstance(a, Eq)]
    others = [a for a in atoms if not isinstance(a, Eq)]
    if len(eqs) != 1 or others:
        return None
    x, y = variables
    diff = sp.expand(eqs[0].lhs - eqs[0].rhs)
    poly = sp.Poly(diff, x, y)
    if poly is None:
        return None
    degs = {sum(mon) for mon, coeff in poly.terms() if coeff != 0}
    if len(degs) != 1:
        return None
    return poly, list(degs)[0]


def solve_binary_eqn(
    expr: sp.Expr, variables: Sequence[sp.Symbol], *, search_bound: int = 200
) -> IntEqnSolveResult | None:
    detected = detect_binary_homog_eqn(expr, variables)
    if detected is None:
        return None
    poly, degree = detected
    x, y = tuple(variables)
    points = []
    for a in range(-search_bound, search_bound + 1):
        for b in range(-search_bound, search_bound + 1):
            if poly.eval({x: a, y: b}) == 0:
                points.append((a, b))
    points = _dedupe_points(points)
    if not points:
        return None
    return IntEqnSolveResult(
        variables=(x, y),
        solutions=points,
        formula=sol_points_to_form((x, y), points),
        method="binary_homogeneous_thue_like_bruteforce",
        complete=False,
        metadata={"degree": degree, "search_bound": search_bound},
    )


def solve_int_methods(expr: sp.Expr, variables: Sequence[sp.Symbol]) -> IntEqnSolveResult | None:
    variables = tuple(variables)
    try:
        from .engine import run_int_solver_pipeline

        pipeline_result = run_int_solver_pipeline(expr, variables)
        if pipeline_result is not None:
            return IntEqnSolveResult(
                variables=pipeline_result.variables,
                solutions=pipeline_result.solutions,
                formula=pipeline_result.formula,
                method=pipeline_result.method,
                complete=pipeline_result.complete,
                metadata={"provenance": pipeline_result.provenance, **pipeline_result.metadata},
            )
    except Exception:
        pass

    for solver in (
        lambda e, vs: solve_sum_of_two_squares(e, vs),
        lambda e, vs: solve_int_fams2(e, vs),
        lambda e, vs: solve_binary_eqn(e, vs),
        lambda e, vs: solve_int_divis(e, vs),
        lambda e, vs: reduce_int_divis(e, vs),
        lambda e, vs: solve_int_branches(e, vs),
        lambda e, vs: solve_int_sys_via_factor(e, vs),
        lambda e, vs: solve_int_branch(e, vs),
        lambda e, vs: solve_int_recursion2(e, vs),
        lambda e, vs: solve_int_pruning(e, vs),
    ):
        try:
            result = solver(expr, variables)
        except Exception:
            result = None
        if result is not None:
            return result
    return None


__all__ = [
    "IntEqnSolveResult",
    "detect_int_lin_elim",
    "solve_int_fams2",
    "solve_int_divis",
    "reduce_int_divis",
    "solve_int_branches",
    "solve_int_sys_via_factor",
    "solve_int_branch",
    "solve_int_recursion2",
    "solve_int_pruning",
    "detect_sum_eqn",
    "solve_sum_of_two_squares",
    "detect_binary_homog_eqn",
    "solve_binary_eqn",
    "solve_int_methods",
]
