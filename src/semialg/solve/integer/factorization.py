from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass
from itertools import product
from math import isqrt

import sympy as sp
from sympy import Eq

from ._common import RECOVERABLE_ERRORS as _RECOVERABLE_ERRORS
from .formula_utils import split_equalities as _split_equalities
from .output_normalization import CanonIntSolveResult, canon_int_result, dedup_int_points


@dataclass(frozen=True)
class IntegerFactorBranch:
    factor: sp.Expr
    branch_formula: sp.Expr
    pruned_by: str | None = None


@dataclass(frozen=True)
class FactorizationAnalysis:
    primitive_equation: sp.Expr
    content: sp.Expr
    factors: tuple[sp.Expr, ...]
    multiplicities: tuple[int, ...]


def norm_factor_int_eqn(
    expr: sp.Expr, variables: Sequence[sp.Symbol]
) -> FactorizationAnalysis | None:
    eqs, _others = _split_equalities(expr)
    if len(eqs) != 1:
        return None
    diff = sp.expand(eqs[0].lhs - eqs[0].rhs)
    try:
        coeff, facs = sp.factor_list(diff)
    except _RECOVERABLE_ERRORS:
        return None
    symbolic_factors = tuple(sp.expand(f) for f, _m in facs if not sp.sympify(f).is_number)
    multiplicities = tuple(int(m) for f, m in facs if not sp.sympify(f).is_number)
    if len(symbolic_factors) <= 1:
        return None
    if coeff == 0:
        return None
    primitive = sp.expand(diff / coeff)
    return FactorizationAnalysis(
        primitive_equation=primitive,
        content=sp.sympify(coeff),
        factors=symbolic_factors,
        multiplicities=multiplicities,
    )


def _extract_simple_bounds(other_atoms: Sequence[sp.Expr], var: sp.Symbol):
    lower = None
    upper = None
    equals = set()
    parity_residue = None
    for atom in other_atoms:
        try:
            simp = sp.simplify(atom)
        except _RECOVERABLE_ERRORS:
            simp = atom
        if isinstance(simp, Eq):
            if simp.lhs == var and simp.rhs.is_integer:
                equals.add(int(simp.rhs))
            elif simp.rhs == var and simp.lhs.is_integer:
                equals.add(int(simp.lhs))
            else:
                expr = sp.expand(simp.lhs - simp.rhs)
                try:
                    poly = sp.Poly(expr, var)
                    if poly.degree() == 1:
                        a = poly.coeff_monomial(var)
                        b = poly.coeff_monomial(1)
                        if all(getattr(v, "is_integer", False) for v in (a, b)) and int(a) % 2 == 1:
                            parity_residue = (-int(b)) % 2
                except _RECOVERABLE_ERRORS:
                    pass
        elif isinstance(simp, sp.GreaterThan):
            if simp.lhs == var and simp.rhs.is_integer:
                lower = max(lower, int(simp.rhs)) if lower is not None else int(simp.rhs)
            elif simp.rhs == var and simp.lhs.is_integer:
                upper = min(upper, int(simp.lhs)) if upper is not None else int(simp.lhs)
        elif isinstance(simp, sp.StrictGreaterThan):
            if simp.lhs == var and simp.rhs.is_integer:
                bound = int(simp.rhs) + 1
                lower = max(lower, bound) if lower is not None else bound
            elif simp.rhs == var and simp.lhs.is_integer:
                bound = int(simp.lhs) - 1
                upper = min(upper, bound) if upper is not None else bound
        elif isinstance(simp, sp.LessThan):
            if simp.lhs == var and simp.rhs.is_integer:
                upper = min(upper, int(simp.rhs)) if upper is not None else int(simp.rhs)
            elif simp.rhs == var and simp.lhs.is_integer:
                lower = max(lower, int(simp.lhs)) if lower is not None else int(simp.lhs)
        elif isinstance(simp, sp.StrictLessThan):
            if simp.lhs == var and simp.rhs.is_integer:
                bound = int(simp.rhs) - 1
                upper = min(upper, bound) if upper is not None else bound
            elif simp.rhs == var and simp.lhs.is_integer:
                bound = int(simp.lhs) + 1
                lower = max(lower, bound) if lower is not None else bound
    if equals:
        valset = sorted(equals)
        return valset[0], valset[-1], equals, parity_residue
    return lower, upper, equals, parity_residue


def _prune_branch_cons(
    factor: sp.Expr, variables: Sequence[sp.Symbol], other_atoms: Sequence[sp.Expr]
) -> str | None:
    free = tuple(v for v in variables if factor.has(v))
    if len(free) != 1:
        return None
    var = free[0]
    try:
        poly = sp.Poly(sp.expand(factor), var)
    except _RECOVERABLE_ERRORS:
        return None
    if poly.degree() != 1:
        return None
    a = sp.simplify(poly.coeff_monomial(var))
    b = sp.simplify(poly.coeff_monomial(1))
    if not (a.is_integer and b.is_integer):
        return None
    lower, upper, equals, parity_residue = _extract_simple_bounds(other_atoms, var)
    try:
        root = sp.Rational(-b, a)
    except _RECOVERABLE_ERRORS:
        return None
    if equals:
        if any(sp.simplify(factor.subs(var, e)) == 0 for e in equals):
            return None
        return "incompatible_fixed_value"
    if lower is not None and upper is not None and lower > upper:
        return "inconsistent_bounds"
    if lower is not None and upper is not None and (root < lower or root > upper):
        return "root_outside_interval"
    if root.q != 1 and lower == upper and lower is not None:
        return "nonintegral_root_at_singleton_interval"
    if parity_residue is not None and root.q == 1 and int(root) % 2 != parity_residue:
        return "parity_conflict"
    return None


def enum_factor_branches(
    expr: sp.Expr, variables: Sequence[sp.Symbol]
) -> list[IntegerFactorBranch]:
    variables = tuple(variables)
    analysis = norm_factor_int_eqn(expr, variables)
    if analysis is None:
        return []
    _eqs, others = _split_equalities(expr)
    branches = []
    seen = set()
    for factor in analysis.factors:
        key = sp.srepr(sp.factor(factor))
        if key in seen:
            continue
        seen.add(key)
        pruned = _prune_branch_cons(factor, variables, others)
        branch_formula = sp.And(*(list(others) + [sp.Eq(factor, 0)]))
        branches.append(
            IntegerFactorBranch(factor=factor, branch_formula=branch_formula, pruned_by=pruned)
        )
    return branches


def _integer_divisors(n: int) -> list[int]:
    if n == 0:
        return []
    n = abs(int(n))
    out = set()
    for d in range(1, isqrt(n) + 1):
        if n % d == 0:
            out.add(d)
            out.add(n // d)
    vals = sorted(out)
    return sorted(set(vals + [-v for v in vals]))


def prune_divisor_assign(
    symbolic_factors: Sequence[sp.Expr],
    assigned_values: Sequence[int],
    other_atoms: Sequence[sp.Expr],
) -> str | None:
    # simple pruning for factors known to equal odd/even affine linear terms
    for factor, value in zip(symbolic_factors, assigned_values, strict=True):
        free = tuple(v for v in factor.free_symbols)
        if len(free) != 1:
            continue
        var = free[0]
        try:
            poly = sp.Poly(sp.expand(factor), var)
        except _RECOVERABLE_ERRORS:
            continue
        if poly.degree() != 1:
            continue
        a = poly.coeff_monomial(var)
        b = poly.coeff_monomial(1)
        if all(getattr(v, "is_integer", False) for v in (a, b)):
            if int(a) % 2 == 0 and (value - int(b)) % 2 != 0:
                return "parity_conflict_with_assignment"
        lower, upper, equals, _parity = _extract_simple_bounds(other_atoms, var)
        if equals:
            if not any(sp.simplify(factor.subs(var, e) - value) == 0 for e in equals):
                return "fixed_value_conflict"
        if lower is not None and upper is not None and lower == upper:
            if sp.simplify(factor.subs(var, lower) - value) != 0:
                return "singleton_interval_conflict"
    return None


def factor_thread_eqn(expr: sp.Expr, variables: Sequence[sp.Symbol]) -> CanonIntSolveResult | None:
    """
    Stronger factor-threading factor-threading for equations of the form f1*...*fk == c
    with integer constant c. Each factor is assigned a divisor of c and the
    resulting branch system is solved recursively by the integer pipeline.
    """
    variables = tuple(variables)
    eqs, others = _split_equalities(expr)
    if len(eqs) != 1:
        return None
    lhs = sp.expand(eqs[0].lhs)
    rhs = sp.expand(eqs[0].rhs)
    const_side = None
    product_side = None
    if rhs.is_integer:
        const_side = int(rhs)
        product_side = sp.factor(lhs)
    elif lhs.is_integer:
        const_side = int(lhs)
        product_side = sp.factor(rhs)
    else:
        return None
    try:
        coeff, facs = sp.factor_list(sp.expand(product_side))
    except _RECOVERABLE_ERRORS:
        return None
    symbolic_factors = [
        sp.expand(f) for f, m in facs for _ in range(int(m)) if not sp.sympify(f).is_number
    ]
    numeric_coeff = int(coeff) if sp.sympify(coeff).is_integer else None
    if numeric_coeff is None or len(symbolic_factors) < 2:
        return None
    target = const_side
    if target % numeric_coeff != 0:
        return canon_int_result(
            variables,
            formula=sp.false,
            solutions=[],
            method="factor_threading_constant_mismatch",
            complete=True,
            provenance=["factorization"],
            metadata={"content": numeric_coeff, "target": target},
        )
    reduced_target = target // numeric_coeff
    divisors = _integer_divisors(reduced_target)
    if not divisors or len(symbolic_factors) > 4 or len(divisors) > 96:
        return None

    from .engine import run_int_solver_pipeline

    points = []
    branch_formulas = []
    pruned_assignments = 0
    for assigned in product(divisors, repeat=len(symbolic_factors) - 1):
        prod_prefix = 1
        for a in assigned:
            prod_prefix *= a
        if prod_prefix == 0 or reduced_target % prod_prefix != 0:
            continue
        last = reduced_target // prod_prefix
        values = list(assigned) + [last]
        pruned = prune_divisor_assign(symbolic_factors, values, others)
        if pruned is not None:
            pruned_assignments += 1
            continue
        branch_eqs = [sp.Eq(f, v) for f, v in zip(symbolic_factors, values, strict=True)]
        branch_expr = sp.And(*(list(others) + branch_eqs))
        sub = run_int_solver_pipeline(branch_expr, variables, search_bound=50)
        if sub is not None:
            if sub.solutions:
                points.extend(sub.solutions)
            branch_formulas.append(sub.formula)
        else:
            branch_formulas.append(branch_expr)

    points = dedup_int_points(points)
    if points:
        return canon_int_result(
            variables,
            solutions=points,
            method="factor_threading_divisor_assignment",
            complete=False,
            provenance=["factorization"],
            metadata={
                "reduced_target": reduced_target,
                "factor_count": len(symbolic_factors),
                "pruned_assignments": pruned_assignments,
            },
        )
    if branch_formulas:
        return canon_int_result(
            variables,
            formula=sp.Or(*branch_formulas),
            solutions=[],
            method="factor_threading_symbolic_branches",
            complete=False,
            provenance=["factorization"],
            metadata={
                "reduced_target": reduced_target,
                "factor_count": len(symbolic_factors),
                "pruned_assignments": pruned_assignments,
            },
        )
    return None


def solve_int_recursion(
    expr: sp.Expr, variables: Sequence[sp.Symbol]
) -> CanonIntSolveResult | None:
    variables = tuple(variables)

    threaded = factor_thread_eqn(expr, variables)
    if threaded is not None:
        return threaded

    branches = enum_factor_branches(expr, variables)
    if not branches:
        return None

    from .engine import run_int_solver_pipeline

    points = []
    symbolic_branches = []
    pruned = []
    for branch in branches:
        if branch.pruned_by is not None:
            pruned.append((sp.srepr(branch.factor), branch.pruned_by))
            continue
        sub = run_int_solver_pipeline(branch.branch_formula, variables, search_bound=50)
        if sub is not None:
            if sub.solutions:
                points.extend(sub.solutions)
            symbolic_branches.append(sub.formula)
        else:
            symbolic_branches.append(branch.branch_formula)

    points = dedup_int_points(points)
    if points:
        return canon_int_result(
            variables,
            solutions=points,
            method="factorization_zero_branching",
            complete=False,
            provenance=["factorization"],
            metadata={"pruned_branches": pruned},
        )
    if symbolic_branches:
        return canon_int_result(
            variables,
            formula=sp.Or(*symbolic_branches),
            solutions=[],
            method="factorization_symbolic_zero_branches",
            complete=False,
            provenance=["factorization"],
            metadata={"pruned_branches": pruned},
        )
    return canon_int_result(
        variables,
        formula=sp.false,
        solutions=[],
        method="factorization_all_branches_pruned",
        complete=True,
        provenance=["factorization"],
        metadata={"pruned_branches": pruned},
    )


__all__ = [
    "IntegerFactorBranch",
    "FactorizationAnalysis",
    "norm_factor_int_eqn",
    "enum_factor_branches",
    "factor_thread_eqn",
    "solve_int_recursion",
]
