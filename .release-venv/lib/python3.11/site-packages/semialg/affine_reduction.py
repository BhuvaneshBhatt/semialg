"""Exact affine reductions with parameter-exception branching."""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import replace

import sympy as sp

from .conditional import ParameterStratifiedResult, conditional_result
from .normalization import conjuncts, normalize_formula, normalize_variables
from .presolve import presolve_semialgebraic


def parametric_affine_reduction(
    formula,
    variables: Sequence[sp.Symbol | str],
    parameters: Sequence[sp.Symbol | str],
) -> ParameterStratifiedResult:
    """Reduce affine equalities, branching exactly where a parameter pivot vanishes."""
    normalized_formula = normalize_formula(formula)
    normalized_variables = normalize_variables(
        variables, normalized_formula, append_context_symbols=False
    )
    normalized_parameters = normalize_variables(
        parameters, normalized_formula, append_context_symbols=False
    )
    branches: list[
        tuple[sp.Expr, sp.Expr, tuple[tuple[sp.Symbol, sp.Expr], ...], tuple[sp.Symbol, ...]]
    ] = [(sp.true, normalized_formula, (), normalized_variables)]
    changed = True
    while changed:
        changed = False
        next_branches = []
        for guard, branch_formula, substitutions, remaining_variables in branches:
            atoms = list(
                conjuncts(branch_formula)
                if isinstance(branch_formula, sp.And)
                else (branch_formula,)
            )
            split_done = False
            for atom in atoms:
                if not isinstance(atom, sp.Equality):
                    continue
                residual = sp.expand(atom.lhs - atom.rhs)
                for variable in remaining_variables:
                    try:
                        polynomial = sp.Poly(residual, variable, domain="EX")
                    except (sp.PolynomialError, TypeError, ValueError):
                        continue
                    if polynomial.degree() != 1:
                        continue
                    pivot = sp.factor(polynomial.coeff_monomial(variable))
                    constant_term = sp.factor(polynomial.coeff_monomial(1))
                    if pivot.free_symbols - set(normalized_parameters):
                        continue
                    solution = sp.cancel(-constant_term / pivot)
                    other_atoms = [other for other in atoms if other != atom]
                    substituted_formula = sp.And(
                        *(sp.simplify(other.subs(variable, solution)) for other in other_atoms)
                    )
                    reduced_variables = tuple(v for v in remaining_variables if v != variable)
                    next_branches.append(
                        (
                            sp.And(guard, sp.Ne(pivot, 0)),
                            substituted_formula,
                            substitutions + ((variable, solution),),
                            reduced_variables,
                        )
                    )
                    zero_pivot_formula = sp.And(*(other_atoms + [sp.Eq(constant_term, 0)]))
                    next_branches.append(
                        (
                            sp.And(guard, sp.Eq(pivot, 0)),
                            zero_pivot_formula,
                            substitutions,
                            remaining_variables,
                        )
                    )
                    split_done = True
                    changed = True
                    break
                if split_done:
                    break
            if not split_done:
                next_branches.append((guard, branch_formula, substitutions, remaining_variables))
        branches = next_branches
        if len(branches) > 64:
            break
    result_branches = []
    for guard, branch_formula, substitutions, remaining_variables in branches:
        presolved = presolve_semialgebraic(
            branch_formula, remaining_variables, eliminate=remaining_variables
        )
        combined_substitutions = substitutions + tuple(presolved.substitutions)
        value = replace(presolved, substitutions=combined_substitutions)
        result_branches.append((guard, value))
    return conditional_result(
        normalized_parameters,
        result_branches,
        method="parametric_affine_pivot_stratification",
        diagnostics={"branch_count": len(result_branches)},
    )


__all__ = ["parametric_affine_reduction"]
