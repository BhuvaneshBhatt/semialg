from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass
from functools import lru_cache

import sympy as sp

from ._common import RECOVERABLE_ERRORS as _RECOVERABLE_ERRORS
from ._common import expr_complexity as _expr_complexity
from .formula_utils import (
    integer_roots as _integer_roots,
)
from .formula_utils import (
    integer_roots_with_completeness as _integer_roots_complete,
)
from .formula_utils import (
    split_equalities as _split_equalities,
)
from .output_normalization import canon_int_result, dedup_int_points


@dataclass(frozen=True)
class GroebnerRecursiveStep:
    chosen_variable: sp.Symbol
    chosen_polynomial: sp.Expr
    groebner_basis: tuple[sp.Expr, ...]
    branch_roots: tuple[sp.Expr, ...]
    leading_linear_relations: tuple[sp.Expr, ...] = ()
    eliminated_variables: tuple[sp.Symbol, ...] = ()
    consistency_constraints: tuple[sp.Expr, ...] = ()


@dataclass(frozen=True)
class GroebnerTriangAnalysis:
    basis: tuple[sp.Expr, ...]
    triangular_polynomials: tuple[sp.Expr, ...]
    variables: tuple[sp.Symbol, ...]
    linear_leading_relations: tuple[sp.Expr, ...] = ()
    eliminated_variables: tuple[sp.Symbol, ...] = ()
    consistency_constraints: tuple[sp.Expr, ...] = ()
    mod_obstr_primes: tuple[int, ...] = ()


def _compute_small_obstr(
    expr: sp.Expr, variables: Sequence[sp.Symbol], primes: Sequence[int] = (2, 3, 5, 7)
) -> dict:
    from .congruence import solve_quant_free_mod_sys

    obstructions = []
    feasible_counts = {}
    for p in primes:
        try:
            result = solve_quant_free_mod_sys(expr, variables, p, max_points=8000)
            feasible_counts[p] = len(result.points)
            if result.points == []:
                obstructions.append(p)
        except _RECOVERABLE_ERRORS:
            continue
    return {"obstructing_primes": tuple(obstructions), "prime_solution_counts": feasible_counts}


@lru_cache(maxsize=256)
def _groebner_basis_cached(polys_key, vars_key):
    polys = [sp.sympify(s) for s in polys_key]
    vars_ = [sp.Symbol(v) for v in vars_key]
    try:
        gb = sp.groebner(polys, *reversed(tuple(vars_)), order="lex")
        return tuple(sp.expand(p.as_expr()) for p in gb.polys if sp.expand(p.as_expr()) != 0)
    except _RECOVERABLE_ERRORS:
        return tuple(polys)


def compute_groebner_basis(
    expr: sp.Expr, variables: Sequence[sp.Symbol]
) -> tuple[sp.Expr, ...] | None:
    eqs, _others = _split_equalities(expr)
    if not eqs:
        return None
    polys = tuple(sp.expand(eq.lhs - eq.rhs) for eq in eqs)
    vars_key = tuple(v.name for v in variables)
    return _groebner_basis_cached(tuple(map(sp.srepr, polys)), vars_key)


def extract_lin_basis(
    basis: Sequence[sp.Expr], variables: Sequence[sp.Symbol]
) -> tuple[sp.Expr, ...]:
    variables = tuple(variables)
    relations = []
    for var in reversed(variables):
        best = None
        for poly_expr in basis:
            try:
                p1 = sp.Poly(sp.expand(poly_expr), var)
            except _RECOVERABLE_ERRORS:
                continue
            if p1.degree() != 1:
                continue
            a = sp.expand(p1.coeff_monomial(var))
            b = sp.expand(p1.coeff_monomial(1))
            if a == 0 or a.has(var) or b.has(var):
                continue
            # Eliminating an affine variable avoids introducing unnecessary algebraic
            # branches into the recursive polynomial system.
            score = (
                _expr_complexity(a) + _expr_complexity(b),
                len(poly_expr.free_symbols),
                sp.srepr(poly_expr),
            )
            rel = sp.Eq(var, sp.simplify(-b / a))
            if best is None or score < best[0]:
                best = (score, rel)
        if best is not None:
            relations.append(best[1])
    by_var = {}
    for rel in relations:
        if isinstance(rel, sp.Equality):
            by_var[rel.lhs] = rel
    return tuple(by_var[v] for v in variables if v in by_var)


def _extract_consist_basis(
    basis: Sequence[sp.Expr], variables: Sequence[sp.Symbol]
) -> tuple[sp.Expr, ...]:
    variables = tuple(variables)
    constraints = []
    for poly_expr in basis:
        free = tuple(v for v in variables if poly_expr.has(v))
        if len(free) == 1:
            var = free[0]
            roots = _integer_roots(poly_expr, var)
            if roots:
                constraints.append(sp.Or(*[sp.Eq(var, r) for r in roots]))
    uniq = []
    seen = set()
    for c in constraints:
        k = sp.srepr(c)
        if k not in seen:
            seen.add(k)
            uniq.append(c)
    return tuple(uniq)


def analyze_groebner_struct(
    expr: sp.Expr, variables: Sequence[sp.Symbol]
) -> GroebnerTriangAnalysis | None:
    variables = tuple(variables)
    basis = compute_groebner_basis(expr, variables)
    if basis is None:
        return None
    triangular = []
    for var in variables:
        candidates = [p for p in basis if p.free_symbols.issubset({var})]
        if candidates:
            triangular.append(min(candidates, key=_expr_complexity))
    linear_rels = extract_lin_basis(basis, variables)
    eliminated = tuple(rel.lhs for rel in linear_rels if isinstance(rel, sp.Equality))
    consistency = _extract_consist_basis(basis, variables)
    obstruction = _compute_small_obstr(expr, variables)
    return GroebnerTriangAnalysis(
        basis=basis,
        triangular_polynomials=tuple(triangular),
        variables=variables,
        linear_leading_relations=linear_rels,
        eliminated_variables=eliminated,
        consistency_constraints=consistency,
        mod_obstr_primes=tuple(obstruction["obstructing_primes"]),
    )


def _choose_recursive_step(
    expr: sp.Expr, variables: Sequence[sp.Symbol], max_branch_points: int
) -> GroebnerRecursiveStep | None:
    variables = tuple(variables)
    analysis = analyze_groebner_struct(expr, variables)
    if analysis is None:
        return None
    best = None
    for var in variables:
        candidates = [p for p in analysis.basis if p.free_symbols.issubset({var})]
        for poly in candidates:
            roots, roots_complete = _integer_roots_complete(poly, var)
            if not roots and roots_complete:
                return GroebnerRecursiveStep(
                    var,
                    poly,
                    analysis.basis,
                    tuple(),
                    analysis.linear_leading_relations,
                    analysis.eliminated_variables,
                    analysis.consistency_constraints,
                )
            if not roots_complete:
                continue
            if len(roots) > max_branch_points:
                continue
            key = (len(roots), _expr_complexity(poly), var.name)
            step = GroebnerRecursiveStep(
                var,
                poly,
                analysis.basis,
                tuple(roots),
                analysis.linear_leading_relations,
                analysis.eliminated_variables,
                analysis.consistency_constraints,
            )
            if best is None or key < best[0]:
                best = (key, step)
    return best[1] if best is not None else None


def recon_int_sol_fams(
    expr: sp.Expr,
    variables: Sequence[sp.Symbol],
    partial_points: Sequence[Sequence[sp.Expr]],
    linear_relations: Sequence[sp.Expr] = (),
) -> sp.Expr:
    pts = dedup_int_points(partial_points)
    if pts:
        return sp.Or(
            *[sp.And(*[sp.Eq(v, val) for v, val in zip(variables, pt, strict=True)]) for pt in pts]
        )
    parts = list(linear_relations)
    if expr not in (sp.true, sp.false):
        parts.append(sp.simplify(expr))
    return sp.And(*parts) if parts else sp.true


def apply_lin_rels_to_expr(expr: sp.Expr, relations: Sequence[sp.Expr]) -> sp.Expr:
    out = expr
    for rel in relations:
        if isinstance(rel, sp.Equality) and rel.lhs.is_Symbol:
            out = sp.simplify(out.subs(rel.lhs, rel.rhs))
    return out


def recon_points_from_rels(
    variables: Sequence[sp.Symbol],
    reduced_vars: Sequence[sp.Symbol],
    reduced_points: Sequence[Sequence[sp.Expr]],
    relations: Sequence[sp.Expr],
) -> list[tuple[sp.Expr, ...]]:
    out = []
    for tail in reduced_points:
        mapping = {v: val for v, val in zip(reduced_vars, tail, strict=True)}
        progress = True
        while progress:
            progress = False
            for rel in relations:
                if isinstance(rel, sp.Equality) and rel.lhs.is_Symbol and rel.lhs not in mapping:
                    rhs = sp.simplify(rel.rhs.subs(mapping))
                    if rhs.free_symbols <= set(mapping.keys()):
                        mapping[rel.lhs] = rhs
                        progress = True
        if all(v in mapping for v in variables):
            out.append(tuple(sp.simplify(mapping[v]) for v in variables))
    return dedup_int_points(out)


def _part_spec_rels(
    basis: Sequence[sp.Expr],
    variables: Sequence[sp.Symbol],
    relations: Sequence[sp.Expr],
):
    mapping = {}
    for rel in relations:
        if isinstance(rel, sp.Equality) and rel.lhs.is_Symbol:
            mapping[rel.lhs] = rel.rhs
    new_basis = tuple(sp.simplify(sp.expand(p.subs(mapping))) for p in basis)
    remaining = tuple(v for v in variables if v not in mapping)
    return new_basis, remaining, mapping


def _choose_coupled_rels(
    basis: Sequence[sp.Expr],
    remaining: Sequence[sp.Symbol],
):
    best = None
    remaining = tuple(remaining)
    for poly in basis:
        free = tuple(v for v in remaining if poly.has(v))
        if len(free) == 2:
            key = (_expr_complexity(poly), len(free), sp.srepr(poly))
            if best is None or key < best[0]:
                best = (key, (poly, free))
    return best[1] if best else None


def _scan_coupled_branches(
    poly: sp.Expr,
    pair_vars: Sequence[sp.Symbol],
    *,
    branch_bound: int = 12,
):
    x, y = tuple(pair_vars)
    pts = []
    for xv in range(-branch_bound, branch_bound + 1):
        univ = sp.Poly(sp.expand(poly.subs(x, xv)), y)
        try:
            roots = univ.all_roots()
        except _RECOVERABLE_ERRORS:
            continue
        for r in roots:
            sr = sp.simplify(r)
            if sr.is_integer is True:
                pts.append((xv, sr))
    return dedup_int_points(pts)


def rec_reduce_sys(
    expr: sp.Expr,
    variables: Sequence[sp.Symbol],
    *,
    max_depth: int = 8,
    max_branch_points: int = 64,
):
    """Recursively reduce an integer polynomial system with Groebner-derived constraints."""
    variables = tuple(variables)
    truth = sp.simplify(expr)
    if truth is sp.false:
        return canon_int_result(
            variables,
            formula=sp.false,
            solutions=[],
            method="groebner_recursive_false",
            complete=True,
            provenance=["groebner"],
        )
    if len(variables) == 0:
        return canon_int_result(
            variables,
            formula=truth,
            solutions=[tuple()] if truth is not sp.false else [],
            method="groebner_recursive_zero_var",
            complete=truth is not sp.false,
            provenance=["groebner"],
        )
    if len(variables) == 1:
        var = variables[0]
        try:
            solset = sp.solveset(expr, var, domain=sp.S.Integers)
            if isinstance(solset, sp.FiniteSet):
                pts = [(s,) for s in solset]
                return canon_int_result(
                    (var,),
                    solutions=pts,
                    method="groebner_recursive_univariate",
                    complete=True,
                    provenance=["groebner"],
                    metadata={"solset": solset},
                )
        except _RECOVERABLE_ERRORS:
            pass

    obstruction = _compute_small_obstr(expr, variables)
    if obstruction["obstructing_primes"]:
        return canon_int_result(
            variables,
            formula=sp.false,
            solutions=[],
            method="groebner_recursive_modular_obstruction",
            complete=True,
            provenance=["groebner"],
            metadata=obstruction,
        )
    if max_depth <= 0:
        return canon_int_result(
            variables,
            formula=sp.And(
                sp.Contains(sp.Tuple(*variables), sp.S.Integers ** len(variables)),
                sp.simplify(expr),
            ),
            method="groebner_recursive_depth_limit",
            complete=False,
            provenance=["groebner"],
        )

    step = _choose_recursive_step(expr, variables, max_branch_points=max_branch_points)
    if step is None:
        return None

    # Stronger linear-leading reconstruction.
    if step.leading_linear_relations:
        reduced_expr = apply_lin_rels_to_expr(expr, step.leading_linear_relations)
        reduced_vars = tuple(v for v in variables if v not in set(step.eliminated_variables))
        if reduced_vars != variables:
            from .engine import run_int_solver_pipeline

            sub = run_int_solver_pipeline(reduced_expr, reduced_vars, search_bound=60)
            if sub is not None:
                points = (
                    recon_points_from_rels(
                        variables, reduced_vars, sub.solutions, step.leading_linear_relations
                    )
                    if sub.solutions
                    else []
                )
                formula = recon_int_sol_fams(
                    variables,
                    points,
                    linear_relations=step.leading_linear_relations,
                    expr=sub.formula,
                )
                return canon_int_result(
                    variables,
                    formula=formula,
                    solutions=points,
                    method="groebner_recursive_linear_reconstruction",
                    complete=bool(sub.complete and points),
                    provenance=["groebner"] + list(getattr(sub, "provenance", [])),
                    metadata={"step": step, "subresult": sub, **obstruction},
                )

    specialized_basis, specialized_remaining, relation_map = _part_spec_rels(
        step.groebner_basis, variables, step.leading_linear_relations
    )
    coupled = _choose_coupled_rels(specialized_basis, specialized_remaining)
    if coupled is not None and len(specialized_remaining) >= 2:
        poly2, pair_vars = coupled
        pair_points = _scan_coupled_branches(
            poly2, pair_vars, branch_bound=min(12, max_branch_points)
        )
        if pair_points:
            rebuilt = []
            for pt in pair_points:
                mapping = {v: val for v, val in zip(pair_vars, pt, strict=True)}
                # fill in eliminated vars from affine relations
                progress = True
                while progress:
                    progress = False
                    for rel in step.leading_linear_relations:
                        if (
                            isinstance(rel, sp.Equality)
                            and rel.lhs.is_Symbol
                            and rel.lhs not in mapping
                        ):
                            rhs = sp.simplify(rel.rhs.subs(mapping))
                            if rhs.free_symbols <= set(mapping.keys()):
                                mapping[rel.lhs] = rhs
                                progress = True
                if all(v in mapping for v in variables if v in mapping or v in pair_vars):
                    full = []
                    ok = True
                    for v in variables:
                        if v in mapping:
                            full.append(sp.simplify(mapping[v]))
                        else:
                            ok = False
                            break
                    if ok:
                        rebuilt.append(tuple(full))
            rebuilt = dedup_int_points(rebuilt)
            if rebuilt:
                return canon_int_result(
                    variables,
                    solutions=rebuilt,
                    method="groebner_recursive_coupled_scan",
                    complete=False,
                    provenance=["groebner", "groebner_coupled_scan"],
                    metadata={"step": step, **obstruction},
                )

        if not step.branch_roots:
            return canon_int_result(
                variables,
                formula=sp.false,
                solutions=[],
                method="groebner_recursive_no_integer_roots",
                complete=True,
                provenance=["groebner"],
                metadata={"step": step, **obstruction},
            )

        from .engine import run_int_solver_pipeline

        points = []
        branch_formulas = []
        remaining = tuple(v for v in variables if v != step.chosen_variable)
        for root in step.branch_roots:
            substituted_expr = sp.simplify(expr.subs(step.chosen_variable, root))
            if step.consistency_constraints:
                substituted_expr = sp.And(substituted_expr, *step.consistency_constraints)
            if not remaining:
                truth = sp.simplify(substituted_expr)
                if truth is not sp.false:
                    points.append((root,))
                continue
            sub = run_int_solver_pipeline(substituted_expr, remaining, search_bound=60)
            if sub is not None:
                if sub.solutions:
                    for tail in sub.solutions:
                        mapping = {
                            step.chosen_variable: root,
                            **{v: val for v, val in zip(remaining, tail, strict=True)},
                        }
                        points.append(tuple(mapping[v] for v in variables))
                branch_formulas.append(sp.And(sp.Eq(step.chosen_variable, root), sub.formula))
            else:
                deeper = rec_reduce_sys(
                    substituted_expr,
                    remaining,
                    max_depth=max_depth - 1,
                    max_branch_points=max_branch_points,
                )
                if deeper is not None:
                    if deeper.solutions:
                        for tail in deeper.solutions:
                            mapping = {
                                step.chosen_variable: root,
                                **{v: val for v, val in zip(remaining, tail, strict=True)},
                            }
                            points.append(tuple(mapping[v] for v in variables))
                    branch_formulas.append(
                        sp.And(sp.Eq(step.chosen_variable, root), deeper.formula)
                    )
                else:
                    branch_formulas.append(
                        sp.And(sp.Eq(step.chosen_variable, root), substituted_expr)
                    )

        points = dedup_int_points(points)
        if points:
            return canon_int_result(
                variables,
                solutions=points,
                method="groebner_recursive_integer_solver",
                complete=False,
                provenance=["groebner"],
                metadata={"step": step, **obstruction},
            )
        return canon_int_result(
            variables,
            formula=sp.Or(*branch_formulas) if branch_formulas else sp.false,
            method="groebner_recursive_symbolic_branch_union",
            complete=False,
            provenance=["groebner"],
            metadata={"step": step, **obstruction},
        )


def find_int_recursion(
    expr: sp.Expr,
    variables: Sequence[sp.Symbol],
    *,
    max_depth: int = 8,
    max_branch_points: int = 64,
):
    result = rec_reduce_sys(
        expr, variables, max_depth=max_depth, max_branch_points=max_branch_points
    )
    if result is None or not result.solutions:
        return None
    pt = result.solutions[0]
    return {v: val for v, val in zip(result.variables, pt, strict=True)}


__all__ = [
    "GroebnerRecursiveStep",
    "GroebnerTriangAnalysis",
    "compute_groebner_basis",
    "extract_lin_basis",
    "analyze_groebner_struct",
    "recon_int_sol_fams",
    "rec_reduce_sys",
    "find_int_recursion",
]
