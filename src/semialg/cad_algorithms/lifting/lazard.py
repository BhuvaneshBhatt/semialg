from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass

import sympy as sp

from ...structural_keys import symbol_identity_key


@dataclass(frozen=True)
class NativeLazardEvaluation:
    """A native Taylor residue and its lexicographically minimal valuation."""

    polynomial: sp.Poly
    variables: tuple[sp.Symbol, ...]
    parent_samples: tuple[object, ...]
    valuation: tuple[int, ...]
    residue: sp.Poly


def _nonzero_over_prefix(poly, variables, prefix):
    from ...algebraic.signs import sign_at_sample

    assigned = variables[: len(prefix)]
    remaining = variables[len(prefix) :]
    coefficients = sp.Poly(poly.as_expr(), *remaining).coeffs()
    return any(
        sign_at_sample(sp.Poly(coefficient, *assigned, domain="EX"), prefix) != 0
        for coefficient in coefficients
    )


def native_lazard_evaluate(poly, variables, parent_samples):
    """Compute Lazard valuation using Hasse derivatives over a selected tower.

    The first nonzero Taylor coefficient in each assigned coordinate is found
    by native coefficient-sign decisions. Keeping the mixed derivative as a
    polynomial in the original variables gives a moving section selector on
    every base cell where the valuation vector is constant. No algebraic sample
    is substituted as a reconstructed primitive element.
    """
    variables = tuple(variables)
    prefix = tuple(parent_samples)
    if len(prefix) >= len(variables):
        raise ValueError("Lazard lifting needs an unassigned fiber variable")
    poly = sp.Poly(poly.as_expr() if isinstance(poly, sp.Poly) else poly, *variables)
    if poly.is_zero:
        raise ValueError("Lazard valuation is defined for a nonzero polynomial")
    current = poly
    orders = []
    for index, variable in enumerate(variables[: len(prefix)]):
        order = 0
        while not _nonzero_over_prefix(current, variables, prefix[: index + 1]):
            order += 1
            current = sp.Poly(
                current.diff(variable).as_expr() / order, *variables, domain=current.domain
            )
            if current.is_zero:
                raise ArithmeticError("nonzero Taylor coefficient could not be established")
        orders.append(order)
    return NativeLazardEvaluation(poly, variables, prefix, tuple(orders), current)


@dataclass(frozen=True)
class LazardStep:
    variable: sp.Symbol
    sample_value: sp.Expr
    cancelled_power: int
    resulting_expr: sp.Expr


@dataclass(frozen=True)
class LazardEvaluationResult:
    original_expr: sp.Expr
    final_expr: sp.Expr
    steps: tuple[LazardStep, ...]
    used_cancellation: bool

    @property
    def valuation(self) -> tuple[int, ...]:
        return tuple(step.cancelled_power for step in self.steps)


def _cancel_linear_factor(expr: sp.Expr, var: sp.Symbol, value: sp.Expr) -> tuple[sp.Expr, int]:
    """Repeatedly cancel ``var - value`` before substituting ``var = value``.

    This is the core local operation in Lazard-style evaluation. It records
    how many times the current expression vanishes on the current section and
    then substitutes into the first non-vanishing quotient. The implementation
    is exact over SymPy's expression domain and is conservative:
    if exact division cannot prove divisibility, cancellation stops.
    """
    factor = sp.expand(var - value)
    current = sp.expand(expr)
    cancelled = 0
    while current != 0:
        gens = tuple(sorted((current.free_symbols | factor.free_symbols), key=symbol_identity_key))
        if not gens:
            break
        try:
            quotient, remainder = sp.div(current, factor, *gens, extension=True)
        except (sp.PolynomialError, TypeError, ValueError, NotImplementedError, ArithmeticError):
            break
        if sp.expand(remainder) != 0:
            break
        current = sp.expand(quotient)
        cancelled += 1
    return current, cancelled


def lazard_evaluate(
    expr: sp.Expr,
    assigned_vars: Sequence[sp.Symbol],
    sample_prefix: Sequence[sp.Expr],
) -> LazardEvaluationResult:
    """Evaluate by Lazard valuation: cancel before each substitution.

    For each pair ``x_i = a_i``, factors ``x_i - a_i`` are cancelled as many
    times as exact division proves possible; the quotient is then specialized.
    The returned valuation tuple is the sequence of cancelled powers. This is
    the evaluation primitive used by Lazard lifting hooks and diagnostics.
    """
    current = sp.expand(expr)
    steps: list[LazardStep] = []
    used_cancellation = False

    for var, value in zip(assigned_vars, sample_prefix, strict=True):
        cancelled_expr, cancelled_power = _cancel_linear_factor(current, var, value)
        if cancelled_power:
            used_cancellation = True
        substituted = sp.expand(cancelled_expr.subs(var, value))
        steps.append(
            LazardStep(
                variable=var,
                sample_value=value,
                cancelled_power=cancelled_power,
                resulting_expr=substituted,
            )
        )
        current = substituted

    return LazardEvaluationResult(
        original_expr=sp.expand(expr),
        final_expr=sp.expand(current),
        steps=tuple(steps),
        used_cancellation=used_cancellation,
    )


def lazard_valuation(
    expr: sp.Expr, assigned_vars: Sequence[sp.Symbol], sample_prefix: Sequence[sp.Expr]
) -> tuple[int, ...]:
    """Return the Lazard valuation vector for a specialization path."""
    return lazard_evaluate(expr, assigned_vars, sample_prefix).valuation


__all__ = [
    "LazardEvaluationResult",
    "LazardStep",
    "NativeLazardEvaluation",
    "lazard_evaluate",
    "lazard_valuation",
    "native_lazard_evaluate",
]
