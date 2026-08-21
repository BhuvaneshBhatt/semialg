from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass, field

import sympy as sp

from ._common import RECOVERABLE_ERRORS as _RECOVERABLE_ERRORS
from ._common import expr_complexity as _expr_complexity
from .formula_utils import conjuncts as _conjuncts
from .output_normalization import CanonIntSolveResult, canon_int_result


@dataclass(frozen=True)
class IntLinElimCand:
    equation: sp.Expr
    solved_variable: sp.Symbol
    coefficient: sp.Expr
    numerator: sp.Expr
    denominator: sp.Expr
    replacement: sp.Expr
    divisibility_condition: sp.Expr
    substituted_atoms: tuple[sp.Expr, ...]
    score: tuple[int, int, str] = field(default_factory=tuple)


def _candidate_score(
    replacement: sp.Expr, substituted_atoms: Sequence[sp.Expr], solved_variable: sp.Symbol
) -> tuple[int, int, str]:
    return (
        _expr_complexity(replacement) + sum(_expr_complexity(a) for a in substituted_atoms),
        len(tuple(substituted_atoms)),
        solved_variable.name,
    )


def enum_int_lin_elim_cands(expr: sp.Expr, variables: Sequence[sp.Symbol]) -> list[IntLinElimCand]:
    variables = tuple(variables)
    conjuncts = _conjuncts(expr)
    candidates: list[IntLinElimCand] = []
    for atom in conjuncts:
        if not isinstance(atom, sp.Equality):
            continue
        diff = sp.expand(atom.lhs - atom.rhs)
        for var in reversed(variables):
            try:
                poly = sp.Poly(diff, var)
            except _RECOVERABLE_ERRORS:
                continue
            if poly.degree() != 1:
                continue
            coeff = sp.expand(poly.coeff_monomial(var))
            const = sp.expand(poly.coeff_monomial(1))
            if coeff == 0 or coeff.has(var) or const.has(var):
                continue
            numerator = sp.expand(-const)
            denominator = sp.expand(coeff)
            replacement = sp.simplify(numerator / denominator)
            try:
                divisibility = sp.Eq(sp.Mod(numerator, denominator), 0)
            except _RECOVERABLE_ERRORS:
                divisibility = sp.Eq(sp.Mod(numerator, sp.Abs(denominator)), 0)
            rest = [a for a in conjuncts if a is not atom]
            substituted = tuple(sp.simplify(a.subs(var, replacement)) for a in rest)
            candidates.append(
                IntLinElimCand(
                    equation=atom,
                    solved_variable=var,
                    coefficient=coeff,
                    numerator=numerator,
                    denominator=denominator,
                    replacement=replacement,
                    divisibility_condition=divisibility,
                    substituted_atoms=substituted,
                    score=_candidate_score(replacement, substituted, var),
                )
            )
    return sorted(candidates, key=lambda c: c.score)


def choose_best_int_lin_elim(
    expr: sp.Expr, variables: Sequence[sp.Symbol]
) -> IntLinElimCand | None:
    candidates = enum_int_lin_elim_cands(expr, variables)
    return candidates[0] if candidates else None


def apply_int_lin_elim(
    expr: sp.Expr, variables: Sequence[sp.Symbol], candidate: IntLinElimCand | None = None
):
    variables = tuple(variables)
    cand = candidate or choose_best_int_lin_elim(expr, variables)
    if cand is None:
        return expr, None, variables
    reduced_expr = (
        sp.And(cand.divisibility_condition, *cand.substituted_atoms)
        if cand.substituted_atoms
        else cand.divisibility_condition
    )
    reduced_vars = tuple(v for v in variables if v != cand.solved_variable)
    return sp.simplify(reduced_expr), cand, reduced_vars


def _terminal_integer_solve(
    expr: sp.Expr, variables: Sequence[sp.Symbol]
) -> CanonIntSolveResult | None:
    variables = tuple(variables)
    truth = sp.simplify(expr)
    if truth is sp.false:
        return canon_int_result(
            variables,
            formula=sp.false,
            solutions=[],
            method="linear_recursion_terminal_false",
            complete=True,
            provenance=["linear_recursion"],
        )
    if truth is sp.true and len(variables) == 0:
        return canon_int_result(
            variables,
            formula=sp.true,
            solutions=[tuple()],
            method="linear_recursion_terminal_true",
            complete=True,
            provenance=["linear_recursion"],
        )
    if len(variables) == 0:
        return canon_int_result(
            variables,
            formula=truth,
            solutions=[tuple()] if truth is not sp.false else [],
            method="linear_recursion_zero_var",
            complete=truth is not sp.false,
            provenance=["linear_recursion"],
        )
    if len(variables) == 1:
        var = variables[0]
        try:
            solset = sp.solveset(expr, var, domain=sp.S.Integers)
            if isinstance(solset, sp.FiniteSet):
                pts = [(v,) for v in solset]
                return canon_int_result(
                    (var,),
                    solutions=pts,
                    method="linear_recursion_univariate",
                    complete=True,
                    provenance=["linear_recursion"],
                    metadata={"solset": solset},
                )
        except _RECOVERABLE_ERRORS:
            pass
        return canon_int_result(
            (var,),
            formula=sp.And(sp.Contains(var, sp.S.Integers), sp.simplify(expr)),
            solutions=[],
            method="linear_recursion_symbolic_univariate",
            complete=False,
            provenance=["linear_recursion"],
        )
    return None


def rec_reduce_int_lin_sys(
    expr: sp.Expr,
    variables: Sequence[sp.Symbol],
    *,
    max_depth: int = 8,
) -> CanonIntSolveResult | None:
    variables = tuple(variables)
    terminal = _terminal_integer_solve(expr, variables)
    if terminal is not None and (terminal.complete or len(variables) <= 1):
        return terminal
    if max_depth <= 0:
        return canon_int_result(
            variables,
            formula=sp.And(
                sp.Contains(sp.Tuple(*variables), sp.S.Integers ** len(variables)),
                sp.simplify(expr),
            ),
            method="linear_recursion_depth_limit",
            complete=False,
            provenance=["linear_recursion"],
            metadata={"depth_limited": True},
        )
    candidate = choose_best_int_lin_elim(expr, variables)
    if candidate is None:
        return terminal

    reduced_expr, cand, reduced_vars = apply_int_lin_elim(expr, variables, candidate)
    sub = rec_reduce_int_lin_sys(reduced_expr, reduced_vars, max_depth=max_depth - 1)
    if sub is None:
        return None

    rebuilt_formula = sp.And(sp.Eq(cand.solved_variable, cand.replacement), sub.formula)
    if sub.solutions:
        rebuilt_points = []
        for tail in sub.solutions:
            mapping = {v: val for v, val in zip(reduced_vars, tail, strict=True)}
            solved_val = sp.simplify(cand.replacement.subs(mapping))
            mapping[cand.solved_variable] = solved_val
            if all((mapping[v]).is_integer is True for v in variables):
                rebuilt_points.append(tuple(mapping[v] for v in variables))
    else:
        rebuilt_points = []

    return canon_int_result(
        variables,
        formula=sp.simplify(rebuilt_formula),
        solutions=rebuilt_points,
        method="recursive_linear_divisibility_solver",
        complete=bool(sub.complete and (sub.solutions or len(reduced_vars) <= 1)),
        provenance=["linear_recursion"] + list(sub.provenance),
        metadata={"candidate": cand, "subresult": sub},
    )


def find_int_recursion2(expr: sp.Expr, variables: Sequence[sp.Symbol], *, max_depth: int = 8):
    result = rec_reduce_int_lin_sys(expr, variables, max_depth=max_depth)
    if result is None:
        return None
    if result.solutions:
        pt = result.solutions[0]
        return {v: val for v, val in zip(result.variables, pt, strict=True)}
    return None


__all__ = [
    "IntLinElimCand",
    "enum_int_lin_elim_cands",
    "choose_best_int_lin_elim",
    "apply_int_lin_elim",
    "rec_reduce_int_lin_sys",
    "find_int_recursion2",
]
