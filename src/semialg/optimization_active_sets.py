"""Active-set pruning and KKT construction for exact polynomial optimization."""

from __future__ import annotations

from collections.abc import Sequence
from itertools import combinations

import sympy as sp
from sympy.polys.polyerrors import PolynomialError

from .relations import split_relation


def _atoms(condition: sp.Expr) -> tuple[sp.Expr, ...]:
    if condition is sp.true or condition == sp.true:
        return ()
    if isinstance(condition, sp.And):
        result: list[sp.Expr] = []
        for arg in condition.args:
            result.extend(_atoms(arg))
        return tuple(result)
    return (condition,)


def jacobian_rank_equations(
    active: Sequence[sp.Expr], variables: Sequence[sp.Symbol]
) -> tuple[sp.Expr, ...]:
    """Return maximal minors forcing active-gradient rank deficiency."""

    row_count = len(active)
    col_count = len(variables)
    if row_count == 0:
        return ()
    jacobian = sp.Matrix([[sp.diff(g, x) for x in variables] for g in active])
    rank = min(row_count, col_count)
    if rank <= 0:
        return ()
    minors: list[sp.Expr] = []
    for rows in combinations(range(row_count), rank):
        for cols in combinations(range(col_count), rank):
            det = sp.expand(jacobian.extract(rows, cols).det())
            if det != 0:
                minors.append(det)
    return tuple(dict.fromkeys(minors))


def canonical_residual(expr: sp.Expr, variables: Sequence[sp.Symbol]) -> sp.Expr:
    """Canonicalize a polynomial zero set up to nonzero rational scale."""

    poly = sp.Poly(sp.expand(expr), *variables, domain=sp.QQ)
    if poly.is_zero:
        return sp.Integer(0)
    return sp.expand(poly.monic().as_expr())


def active_equations_consistent(active: Sequence[sp.Expr], variables: Sequence[sp.Symbol]) -> bool:
    """Return false only when a cheap exact Gröbner check proves inconsistency."""

    if not active:
        return True
    try:
        basis = sp.groebner(tuple(active), *variables, order="grevlex", domain=sp.QQ)
    except (PolynomialError, ValueError, TypeError):
        return True
    var_set = set(variables)
    return not any(
        poly.as_expr().free_symbols.isdisjoint(var_set) and poly.as_expr() != 0
        for poly in basis.polys
    )


def constant_gradient_rank(equalities: Sequence[sp.Expr], variables: Sequence[sp.Symbol]) -> int:
    """Return a globally certified lower bound on equality-gradient rank."""

    if not equalities or not variables:
        return 0
    jacobian = sp.Matrix([[sp.diff(eq, var) for var in variables] for eq in equalities])
    max_rank = min(jacobian.rows, jacobian.cols)
    for rank in range(max_rank, 0, -1):
        for rows in combinations(range(jacobian.rows), rank):
            for cols in combinations(range(jacobian.cols), rank):
                det = sp.expand(jacobian.extract(rows, cols).det())
                if det != 0 and not det.free_symbols:
                    return rank
    return 0


def strict_boundary_keys(condition: sp.Expr | None, variables: Sequence[sp.Symbol]) -> set[str]:
    """Return canonical keys for strict-inequality boundaries."""

    if condition is None:
        return set()
    keys: set[str] = set()
    var_set = set(variables)
    for atom in _atoms(condition):
        try:
            residual, op = split_relation(atom)
        except (TypeError, ValueError, NotImplementedError):
            continue
        if op not in {"<", ">"} or not residual.free_symbols <= var_set:
            continue
        try:
            residual = canonical_residual(residual, variables)
        except (PolynomialError, ValueError, TypeError):
            residual = sp.expand(residual)
        keys.add(sp.srepr(residual))
    return keys


def pruned_active_subsets(
    equalities: tuple[sp.Expr, ...],
    inequalities: tuple[sp.Expr, ...],
    variables: tuple[sp.Symbol, ...],
    condition: sp.Expr | None = None,
) -> tuple[tuple[sp.Expr, ...], ...]:
    """Return nonredundant, algebraically consistent inequality active sets."""

    canonical_eqs: list[sp.Expr] = []
    equality_keys: set[str] = set()
    for eq in equalities:
        try:
            canon = canonical_residual(eq, variables)
        except (PolynomialError, ValueError, TypeError):
            canon = sp.expand(eq)
        key = sp.srepr(canon)
        if canon != 0 and key not in equality_keys:
            canonical_eqs.append(canon)
            equality_keys.add(key)

    equality_basis = None
    if canonical_eqs:
        try:
            equality_basis = sp.groebner(canonical_eqs, *variables, order="grevlex", domain=sp.QQ)
        except (PolynomialError, ValueError, TypeError):
            equality_basis = None
    filtered: list[sp.Expr] = []
    seen: set[str] = set()
    for inequality in inequalities:
        try:
            canon = canonical_residual(inequality, variables)
        except (PolynomialError, ValueError, TypeError):
            canon = sp.expand(inequality)
        if canon == 0:
            continue
        if equality_basis is not None:
            try:
                _, remainder = equality_basis.reduce(canon)
                if sp.expand(remainder) == 0:
                    continue
            except (PolynomialError, ValueError, TypeError):
                pass
        key = sp.srepr(canon)
        if key not in seen:
            filtered.append(canon)
            seen.add(key)

    eq_rank = constant_gradient_rank(canonical_eqs, variables)
    max_active = min(max(0, len(variables) - eq_rank), len(filtered))
    strict_keys = strict_boundary_keys(condition, variables)
    subsets: list[tuple[sp.Expr, ...]] = []
    infeasible: list[frozenset[str]] = []
    for size in range(max_active + 1):
        for subset in combinations(filtered, size):
            keys = frozenset(sp.srepr(item) for item in subset)
            if keys & strict_keys or any(bad <= keys for bad in infeasible):
                continue
            active = tuple(canonical_eqs) + tuple(subset)
            if not active_equations_consistent(active, variables):
                infeasible.append(keys)
                continue
            subsets.append(tuple(subset))
    return tuple(subsets)


def kkt_system(
    objective: sp.Expr,
    variables: tuple[sp.Symbol, ...],
    active: tuple[sp.Expr, ...],
) -> tuple[tuple[sp.Expr, ...], tuple[sp.Symbol, ...]]:
    """Construct polynomial KKT equations for one active set."""

    multipliers = tuple(sp.Symbol(f"_semialg_lambda_{i}", real=True) for i in range(len(active)))
    stationarity: list[sp.Expr] = []
    for var in variables:
        rhs = sum(
            (lam * sp.diff(g, var) for lam, g in zip(multipliers, active, strict=True)),
            sp.Integer(0),
        )
        stationarity.append(sp.expand(sp.diff(objective, var) - rhs))
    return (*active, *stationarity), multipliers
