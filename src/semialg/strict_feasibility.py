"""Exact strict-feasibility and affine-relative-interior primitives."""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field

import sympy as sp

from .decision import implies, is_satisfiable
from .normalization import conjuncts, normalize_formula, normalize_problem_variables


@dataclass(frozen=True)
class StrictFeasibilityResult:
    feasible: bool
    formula: sp.Expr
    strict_formula: sp.Expr
    variables: tuple[sp.Symbol, ...]
    witness: Mapping[sp.Symbol, sp.Expr] | None = None
    method: str = "affine_relative_strict_feasibility"
    diagnostics: Mapping[str, object] = field(default_factory=dict)

    def __bool__(self) -> bool:
        return self.feasible


def _residual(relation: sp.Rel) -> sp.Expr:
    return sp.expand(relation.lhs - relation.rhs)


def _strictify_relation(
    relation: sp.Expr, affine_equalities: sp.Expr, variables: tuple[sp.Symbol, ...]
) -> sp.Expr:
    if isinstance(relation, (sp.Equality, sp.Unequality)):
        return relation
    if not isinstance(
        relation, (sp.LessThan, sp.StrictLessThan, sp.GreaterThan, sp.StrictGreaterThan)
    ):
        return relation
    residual = _residual(relation)
    # An inequality that is identically tight on the explicit affine hull must
    # remain non-strict; otherwise strictification would incorrectly erase the
    # relative interior of a lower-dimensional polyhedron.
    if affine_equalities is not sp.true and implies(
        affine_equalities, sp.Eq(residual, 0), variables
    ):
        return sp.Eq(residual, 0)
    if isinstance(relation, (sp.LessThan, sp.StrictLessThan)):
        return sp.Lt(relation.lhs, relation.rhs)
    return sp.Gt(relation.lhs, relation.rhs)


def affine_relative_interior_formula(
    constraints,
    variables: Sequence[sp.Symbol | str] | None = None,
) -> sp.Expr:
    """Return the relative-interior formula for a conjunctive affine system.

    Explicit affine equalities define the affine hull. Inequalities that are
    identically tight on that hull remain equalities; all other affine
    inequalities are made strict.
    """

    formula = normalize_formula(constraints)
    normalized_variables = normalize_problem_variables(variables, formula)
    atoms = conjuncts(formula) if isinstance(formula, sp.And) else (formula,)
    for atom in atoms:
        if not getattr(atom, "is_Relational", False):
            raise NotImplementedError(
                "relative interior requires a conjunction of affine relations"
            )
        try:
            polynomial = sp.Poly(_residual(atom), *normalized_variables, domain="EX")
        except (sp.PolynomialError, TypeError, ValueError) as exc:
            raise NotImplementedError(
                "relative interior supports affine polynomial relations"
            ) from exc
        if polynomial.total_degree() > 1:
            raise NotImplementedError("relative interior supports affine polynomial relations")
    equality_atoms = tuple(atom for atom in atoms if isinstance(atom, sp.Equality))
    affine_equalities = sp.And(*equality_atoms) if equality_atoms else sp.true
    return sp.And(
        *(_strictify_relation(atom, affine_equalities, normalized_variables) for atom in atoms)
    )


def strict_feasible(
    constraints,
    variables: Sequence[sp.Symbol | str] | None = None,
    *,
    relative: bool = True,
    return_result: bool = False,
) -> bool | StrictFeasibilityResult:
    """Decide strict feasibility, using affine-relative semantics by default."""

    formula = normalize_formula(constraints)
    normalized_variables = normalize_problem_variables(variables, formula)
    if relative:
        strict_formula = affine_relative_interior_formula(formula, normalized_variables)
    else:
        atoms = conjuncts(formula) if isinstance(formula, sp.And) else (formula,)
        strict_formula = sp.And(
            *(
                sp.Lt(atom.lhs, atom.rhs)
                if isinstance(atom, sp.LessThan)
                else sp.Gt(atom.lhs, atom.rhs)
                if isinstance(atom, sp.GreaterThan)
                else atom
                for atom in atoms
            )
        )
    satisfiability = is_satisfiable(strict_formula, normalized_variables, return_result=True)
    result = StrictFeasibilityResult(
        bool(satisfiability),
        formula,
        strict_formula,
        normalized_variables,
        witness=satisfiability.witness,
        method=f"{satisfiability.method}+{'relative' if relative else 'ambient'}_strictification",
    )
    return result if return_result else result.feasible


__all__ = ["StrictFeasibilityResult", "affine_relative_interior_formula", "strict_feasible"]
