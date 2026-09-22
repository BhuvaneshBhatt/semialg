"""Exact critical-locus geometry for polynomial/rational parameterizations."""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass
from itertools import combinations

import sympy as sp

from .decision import is_satisfiable
from .geometry_queries import semialgebraic_image
from .normalization import normalize_formula, normalize_problem_variables
from .regions.operations import region_dimension


@dataclass(frozen=True)
class ParameterizationGeometry:
    mapping: tuple[sp.Expr, ...]
    parameters: tuple[sp.Symbol, ...]
    domain: sp.Expr
    generic_rank: int
    domain_dimension: int
    image_dimension: int
    generic_fiber_dimension: int
    critical_locus: sp.Expr
    critical_value_set: sp.Expr
    image_variables: tuple[sp.Symbol, ...]


def _rank_at_least_formula(jacobian: sp.Matrix, rank: int) -> sp.Expr:
    if rank == 0:
        return sp.true
    minors = []
    for rows in combinations(range(jacobian.rows), rank):
        for cols in combinations(range(jacobian.cols), rank):
            minors.append(sp.Ne(sp.expand(jacobian.extract(rows, cols).det()), 0))
    return sp.Or(*minors) if minors else sp.false


def parameterization_geometry(
    mapping,
    domain=sp.true,
    parameters: Sequence[sp.Symbol | str] | None = None,
    *,
    image_variables=None,
) -> ParameterizationGeometry:
    """Analyze generic rank, rank-drop locus, image, and generic fiber dimension."""
    maps = tuple(map(sp.sympify, mapping))
    condition = normalize_formula(domain)
    params = normalize_problem_variables(parameters, sp.Tuple(*maps), condition)
    denominators = []
    for expression in maps:
        numerator, denominator = sp.fraction(sp.cancel(expression))
        try:
            sp.Poly(numerator, *params)
            sp.Poly(denominator, *params)
        except sp.PolynomialError as exc:
            raise ValueError("mapping must be polynomial or rational in the parameters") from exc
        if denominator != 1:
            denominators.append(sp.Ne(denominator, 0))
    condition = sp.And(condition, *denominators)
    domain_dim = region_dimension(condition, params)
    if image_variables is None:
        targets = sp.symbols(f"y0:{len(maps)}", real=True)
    else:
        targets = tuple(image_variables)
    image = semialgebraic_image(maps, condition, params, image_variables=targets)
    image_dim = region_dimension(image, targets)
    generic_rank = image_dim

    map_jacobian = sp.Matrix(maps).jacobian(params)
    codimension = len(params) - domain_dim
    if codimension == 0:
        rank_drop = (
            sp.And(condition, sp.Not(_rank_at_least_formula(map_jacobian, generic_rank)))
            if generic_rank
            else sp.false
        )
    else:
        # Recover a certified equality presentation of the domain tangent space.
        # This covers algebraic parameter manifolds such as probability simplices.
        equality_polynomials = []
        for atom in condition.atoms(sp.Equality):
            polynomial = sp.expand(atom.lhs - atom.rhs)
            try:
                sp.Poly(polynomial, *params)
            except sp.PolynomialError:
                continue
            equality_polynomials.append(polynomial)
        if not equality_polynomials:
            raise NotImplementedError(
                "lower-dimensional parameter domains require polynomial equalities that "
                "certify their intrinsic tangent spaces"
            )
        constraint_jacobian = sp.Matrix(equality_polynomials).jacobian(params)
        regular_domain = _rank_at_least_formula(constraint_jacobian, codimension)
        if is_satisfiable(sp.And(condition, regular_domain), params, strategy="cad") is not True:
            raise NotImplementedError(
                "the supplied equality constraints do not certify the parameter-domain codimension"
            )
        augmented = constraint_jacobian.col_join(map_jacobian)
        required_rank = codimension + generic_rank
        rank_drop = sp.And(
            condition,
            sp.Or(
                sp.Not(regular_domain),
                sp.Not(_rank_at_least_formula(augmented, required_rank)),
            ),
        )
    critical_values = (
        semialgebraic_image(maps, rank_drop, params, image_variables=targets)
        if rank_drop is not sp.false
        else sp.false
    )
    return ParameterizationGeometry(
        maps,
        params,
        condition,
        generic_rank,
        domain_dim,
        image_dim,
        domain_dim - image_dim,
        sp.simplify(rank_drop),
        sp.simplify(critical_values),
        tuple(targets),
    )


def parameterization_critical_locus(mapping, domain=sp.true, parameters=None) -> sp.Expr:
    """Return the exact generic-rank-drop locus of a parameterization."""
    return parameterization_geometry(mapping, domain, parameters).critical_locus


def parameterization_critical_values(
    mapping, domain=sp.true, parameters=None, *, image_variables=None
) -> sp.Expr:
    """Return the exact image of the parameterization's critical locus."""
    return parameterization_geometry(
        mapping, domain, parameters, image_variables=image_variables
    ).critical_value_set


__all__ = [
    "ParameterizationGeometry",
    "parameterization_geometry",
    "parameterization_critical_locus",
    "parameterization_critical_values",
]
