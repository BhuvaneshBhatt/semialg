"""Exact local geometry of polynomial semialgebraic sets."""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass

import sympy as sp

from ._zero_testing import certified_equal
from .decision import implies
from .normalization import normalize_formula, normalize_problem_variables
from .region_analysis import local_dimension
from .regions.operations import region_closure, region_dimension


@dataclass(frozen=True)
class PolynomialConstraint:
    polynomial: sp.Expr
    relation: str
    atom: sp.Expr


@dataclass(frozen=True)
class PolynomialConstraintClause:
    constraints: tuple[PolynomialConstraint, ...]


@dataclass(frozen=True)
class PolynomialConstraintSystem:
    variables: tuple[sp.Symbol, ...]
    clauses: tuple[PolynomialConstraintClause, ...]
    formula: sp.Expr


@dataclass(frozen=True)
class ActiveConstraintResult:
    point: tuple[sp.Expr, ...]
    constraints: tuple[PolynomialConstraint, ...]
    satisfied_clauses: tuple[int, ...]


def _constraint_from_atom(atom: sp.Expr, variables: tuple[sp.Symbol, ...]) -> PolynomialConstraint:
    if isinstance(atom, sp.Equality):
        relation = "eq"
    elif isinstance(atom, sp.Unequality):
        relation = "ne"
    elif isinstance(atom, sp.GreaterThan):
        relation = "ge"
    elif isinstance(atom, sp.StrictGreaterThan):
        relation = "gt"
    elif isinstance(atom, sp.LessThan):
        relation = "le"
    elif isinstance(atom, sp.StrictLessThan):
        relation = "lt"
    else:
        raise TypeError(f"not a polynomial relation: {atom!r}")
    polynomial = sp.expand(atom.lhs - atom.rhs)
    try:
        sp.Poly(polynomial, *variables)
    except sp.PolynomialError as exc:
        raise ValueError(f"constraint is not polynomial in the problem variables: {atom}") from exc
    return PolynomialConstraint(polynomial, relation, atom)


def polynomial_constraints(
    region, variables: Sequence[sp.Symbol | str] | None = None
) -> PolynomialConstraintSystem:
    """Return a DNF-preserving structured polynomial constraint description.

    Boolean disjunction is represented by separate clauses; conjunction is
    represented inside a clause. Negations are pushed to relational atoms.
    """
    formula = normalize_formula(region)
    vars_ = normalize_problem_variables(variables, formula)
    nnf = sp.to_nnf(formula, simplify=True)
    dnf = sp.to_dnf(nnf, simplify=False, force=True)
    pieces = dnf.args if isinstance(dnf, sp.Or) else (dnf,)
    clauses = []
    for piece in pieces:
        if piece is sp.false:
            continue
        atoms = piece.args if isinstance(piece, sp.And) else (piece,)
        constraints = []
        for atom in atoms:
            if atom is sp.true:
                continue
            constraints.append(_constraint_from_atom(atom, vars_))
        clauses.append(PolynomialConstraintClause(tuple(constraints)))
    return PolynomialConstraintSystem(vars_, tuple(clauses), formula)


def _point_substitution(point, variables):
    if isinstance(point, Mapping):
        return {v: sp.sympify(point[v]) for v in variables}
    values = tuple(map(sp.sympify, point))
    if len(values) != len(variables):
        raise ValueError("point dimension does not match variables")
    return dict(zip(variables, values, strict=True))


def _atom_truth(atom, substitution) -> bool:
    value = sp.simplify(atom.subs(substitution))
    if value in (True, sp.true):
        return True
    if value in (False, sp.false):
        return False
    raise ValueError(f"constraint truth is undecidable at the supplied exact point: {atom}")


def active_constraints(region, point, variables=None) -> ActiveConstraintResult:
    """Return polynomial constraints active at a feasible exact point."""
    system = polynomial_constraints(region, variables)
    subs = _point_substitution(point, system.variables)
    satisfied = []
    active = []
    seen = set()
    for index, clause in enumerate(system.clauses):
        if all(_atom_truth(c.atom, subs) for c in clause.constraints):
            satisfied.append(index)
            for constraint in clause.constraints:
                if constraint.relation == "ne":
                    continue
                value = sp.expand(constraint.polynomial.subs(subs))
                if certified_equal(value, 0) is True:
                    key = (constraint.polynomial, constraint.relation)
                    if key not in seen:
                        seen.add(key)
                        active.append(constraint)
    if not satisfied:
        raise ValueError("point is not in the semialgebraic region")
    return ActiveConstraintResult(
        tuple(subs[v] for v in system.variables), tuple(active), tuple(satisfied)
    )


def _certified_affine_hull(region, variables):
    """Recover an affine hull from certified globally valid linear equalities."""
    formula = normalize_formula(region)
    candidates = []
    for atom in formula.atoms(sp.Equality):
        polynomial = sp.expand(atom.lhs - atom.rhs)
        try:
            poly = sp.Poly(polynomial, *variables)
        except sp.PolynomialError:
            continue
        if (
            poly.total_degree() <= 1
            and implies(formula, sp.Eq(polynomial, 0), variables, strategy="cad") is True
        ):
            candidates.append(sp.Eq(polynomial, 0))
    hull = sp.And(*candidates) if candidates else sp.true
    if region_dimension(hull, variables) != region_dimension(formula, variables):
        raise NotImplementedError(
            "could not certify the affine hull from the region's polynomial equalities"
        )
    return hull


def relative_interior(region, variables=None) -> sp.Expr:
    """Return the interior relative to the certified affine hull of ``region``."""
    formula = normalize_formula(region)
    vars_ = normalize_problem_variables(variables, formula)
    hull = _certified_affine_hull(formula, vars_)
    relative_complement = sp.And(hull, sp.Not(formula))
    return sp.simplify(sp.And(formula, sp.Not(region_closure(relative_complement, vars_))))


def relative_boundary(region, variables=None) -> sp.Expr:
    """Return the boundary relative to the certified affine hull of ``region``."""
    formula = normalize_formula(region)
    vars_ = normalize_problem_variables(variables, formula)
    hull = _certified_affine_hull(formula, vars_)
    return sp.simplify(
        sp.And(
            hull,
            region_closure(sp.And(hull, formula), vars_),
            region_closure(sp.And(hull, sp.Not(formula)), vars_),
        )
    )


def semialgebraic_tangent_cone(region, point, variables=None) -> sp.Expr:
    """Return the exact Bouligand tangent cone at a point of a semialgebraic set.

    The cone is the closure of secant velocities ``(x-p)/t`` for ``x`` in the
    set and ``t > 0``. Quantifier elimination makes the construction exact.
    """
    from .geometry_queries import semialgebraic_projection

    formula = normalize_formula(region)
    vars_ = normalize_problem_variables(variables, formula)
    subs = _point_substitution(point, vars_)
    if not _atom_truth(formula, subs):
        raise ValueError("point is not in the semialgebraic region")
    velocities = sp.symbols(f"_v0:{len(vars_)}", real=True)
    t = sp.Dummy("t", real=True)
    replacement = {x: subs[x] + t * v for x, v in zip(vars_, velocities, strict=True)}
    secants = sp.And(t > 0, formula.xreplace(replacement))
    projected = semialgebraic_projection(secants, (t,), (*velocities, t))
    cone = region_closure(projected, velocities)
    return sp.simplify(cone.xreplace(dict(zip(velocities, vars_, strict=True))))


__all__ = [
    "PolynomialConstraint",
    "PolynomialConstraintClause",
    "PolynomialConstraintSystem",
    "ActiveConstraintResult",
    "polynomial_constraints",
    "active_constraints",
    "relative_interior",
    "relative_boundary",
    "semialgebraic_tangent_cone",
    "local_dimension",
]
