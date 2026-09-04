"""Local algebraic-geometry utilities for real polynomial varieties."""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from itertools import combinations

import sympy as sp

from .algebraic_decomposition import (
    AlgebraicLocusPiece,
    DecompositionCertificate,
    DecompositionPieceCertificate,
    EquidimensionalDecomposition,
    equidimensional_decomposition,
    verify_decomposition_certificate,
)
from .internal_symbols import fresh_dummy
from .normalization import normalize_point, normalize_problem_variables
from .optimization_geometry import polynomial_locus_dimension


def _equations(equations) -> tuple[sp.Expr, ...]:
    if isinstance(equations, (sp.Expr, sp.Equality)):
        equations = (equations,)
    out = []
    for eq in equations:
        expr = sp.sympify(eq)
        if isinstance(expr, sp.Equality):
            expr = sp.expand(expr.lhs - expr.rhs)
        out.append(sp.expand(expr))
    return tuple(out)


def singular_locus(equations, variables=None, *, codimension: int | None = None) -> sp.Expr:
    """Return equations defining the singular locus of an algebraic variety.

    The Jacobian criterion is used with the expected codimension inferred from
    ``polynomial_locus_dimension`` unless supplied explicitly.
    """

    eqs = _equations(equations)
    vars_ = normalize_problem_variables(variables, sp.Tuple(*eqs))
    if codimension is None:
        dim = polynomial_locus_dimension(eqs, vars_)
        if dim is None:
            raise NotImplementedError(
                "could not infer variety dimension; pass codimension explicitly"
            )
        codimension = max(0, len(vars_) - dim)
    if codimension < 0 or codimension > len(vars_):
        raise ValueError("invalid codimension")
    base = [sp.Eq(eq, 0) for eq in eqs]
    if codimension == 0:
        return sp.false
    jac = sp.Matrix([[sp.diff(eq, v) for v in vars_] for eq in eqs])
    if codimension > min(jac.rows, jac.cols):
        return sp.And(*base)
    minors = []
    for rows in combinations(range(jac.rows), codimension):
        for cols in combinations(range(jac.cols), codimension):
            minors.append(sp.expand(jac.extract(rows, cols).det()))
    if not minors:
        return sp.And(*base)
    return sp.And(*base, *(sp.Eq(m, 0) for m in minors))


@dataclass(frozen=True)
class TangentSpaceResult:
    variables: tuple[sp.Symbol, ...]
    point: Mapping[sp.Symbol, sp.Expr]
    jacobian: sp.Matrix
    basis: tuple[sp.Matrix, ...]
    dimension: int
    equations: tuple[sp.Expr, ...]


def tangent_space(equations, point, variables=None) -> TangentSpaceResult:
    """Return the Zariski tangent space at a point as the Jacobian nullspace."""

    eqs = _equations(equations)
    vars_ = normalize_problem_variables(variables, sp.Tuple(*eqs))
    point_map = normalize_point(point, vars_, context=(sp.Tuple(*eqs),))
    for eq in eqs:
        if sp.simplify(eq.subs(point_map)) != 0:
            raise ValueError("point does not lie on the algebraic variety")
    jac = sp.Matrix([[sp.diff(eq, v).subs(point_map) for v in vars_] for eq in eqs])
    basis = tuple(jac.nullspace())
    directions = tuple(sp.Symbol(f"d{i}", real=True) for i in range(len(vars_)))
    linear = tuple(
        sp.expand(sum(jac[i, j] * directions[j] for j in range(len(vars_))))
        for i in range(jac.rows)
    )
    return TangentSpaceResult(vars_, point_map, jac, basis, len(basis), linear)


@dataclass(frozen=True)
class TangentConeResult:
    variables: tuple[sp.Symbol, ...]
    point: Mapping[sp.Symbol, sp.Expr]
    direction_variables: tuple[sp.Symbol, ...]
    initial_forms: tuple[sp.Expr, ...]
    formula: sp.Expr
    certified: bool
    method: str

    @property
    def ideal_generators(self) -> tuple[sp.Expr, ...]:
        """Reduced Groebner-basis generators of the exact tangent-cone ideal."""
        return self.initial_forms


def _lowest_total_degree(poly: sp.Expr, variables: Sequence[sp.Symbol]) -> int:
    p = sp.Poly(sp.expand(poly), *variables)
    nonzero = [(monom, coeff) for monom, coeff in p.terms() if coeff != 0]
    if not nonzero:
        raise ValueError("zero polynomial has no order of vanishing")
    return min(sum(monom) for monom, _ in nonzero)


def _tangent_deformation(
    poly: sp.Expr,
    variables: Sequence[sp.Symbol],
    t: sp.Symbol,
) -> sp.Expr:
    """Return ``t**(-ord(poly)) * poly(t*x)`` as an exact polynomial."""

    order = _lowest_total_degree(poly, variables)
    scaled = sp.expand(poly.subs({v: t * v for v in variables}))
    quotient = sp.cancel(scaled / t**order)
    result = sp.expand(quotient)
    # This division is exact by construction.  Check it instead of allowing a
    # rational expression to leak into the Groebner computation.
    sp.Poly(result, t, *variables)
    return result


def _groebner_exact(polys: Sequence[sp.Expr], variables: Sequence[sp.Symbol]):
    """Compute an exact Groebner basis, allowing algebraic coefficients."""

    try:
        return sp.groebner(tuple(polys), *variables, order="lex", extension=True)
    except (sp.polys.polyerrors.CoercionFailed, sp.polys.polyerrors.GeneratorsError):
        # ``EX`` remains exact for symbolic coefficient expressions when SymPy
        # cannot place all exact coordinates in one algebraic extension.
        return sp.groebner(tuple(polys), *variables, order="lex", domain=sp.EX)


def _saturate_by_parameter(
    generators: Sequence[sp.Expr],
    parameter: sp.Symbol,
    variables: Sequence[sp.Symbol],
) -> tuple[sp.Expr, ...]:
    """Return generators of ``<generators> : parameter**infinity``.

    The standard Rabinowitsch elimination identity is used::

        I : t^infinity = (I + <1 - u*t>) intersect k[t, x].

    With lexicographic order and ``u`` first, the Groebner-basis elements not
    containing ``u`` generate the elimination ideal exactly.
    """

    u = fresh_dummy("tangent_cone_saturation")
    all_variables = (u, parameter, *variables)
    basis = _groebner_exact((*generators, 1 - u * parameter), all_variables)
    eliminated = tuple(
        sp.expand(poly.as_expr()) for poly in basis.polys if u not in poly.as_expr().free_symbols
    )
    if not eliminated:
        return (sp.Integer(0),)
    return eliminated


def _exact_tangent_cone_ideal(
    translated_generators: Sequence[sp.Expr],
    directions: Sequence[sp.Symbol],
) -> tuple[sp.Expr, ...]:
    """Compute ``in_m(I)`` by saturated m-adic deformation."""

    expanded = tuple(sp.expand(f) for f in translated_generators)
    nonzero = tuple(f for f in expanded if f != 0)
    if not nonzero:
        return ()

    t = fresh_dummy("tangent_cone_parameter")
    deformation = tuple(_tangent_deformation(f, directions, t) for f in nonzero)
    saturated = _saturate_by_parameter(deformation, t, directions)

    # The special fibre t=0 is the associated-graded (tangent-cone) ideal.
    fibre = tuple(sp.expand(g.subs(t, 0)) for g in saturated)
    fibre = tuple(g for g in fibre if g != 0)
    if not fibre:
        return ()

    basis = _groebner_exact(fibre, tuple(directions))
    if basis.polys == [sp.Poly(1, *directions)]:
        return (sp.Integer(1),)
    return tuple(sp.expand(poly.as_expr()) for poly in basis.polys)


def tangent_cone(equations, point, variables=None) -> TangentConeResult:
    """Return the exact ideal-theoretic Zariski tangent cone at ``point``.

    If ``I`` is the translated defining ideal and ``m`` is the maximal ideal at
    the origin, the tangent-cone ideal is the associated-graded initial ideal
    ``in_m(I)``.  It is computed exactly by the flat ``m``-adic deformation

    ``t**(-ord(f)) * f(t*d)``

    followed by saturation with respect to ``t`` and specialization at
    ``t = 0``.  Saturation is essential: it captures initial forms arising from
    cancellations between the supplied generators, so the result depends on
    the ideal rather than on a particular generating set.

    ``initial_forms`` holds a
    reduced Groebner basis for the exact tangent-cone ideal.  ``certified`` is
    always ``True``.
    """

    eqs = _equations(equations)
    vars_ = normalize_problem_variables(variables, sp.Tuple(*eqs))
    point_map = normalize_point(point, vars_, context=(sp.Tuple(*eqs),))
    for eq in eqs:
        if sp.simplify(eq.subs(point_map)) != 0:
            raise ValueError("point does not lie on the algebraic variety")

    dirs = tuple(sp.Symbol(f"d{i}", real=True) for i, _ in enumerate(vars_))
    shift = {v: point_map[v] + d for v, d in zip(vars_, dirs, strict=True)}
    translated = tuple(sp.expand(eq.subs(shift)) for eq in eqs)
    cone_ideal = _exact_tangent_cone_ideal(translated, dirs)
    formula = sp.And(*(sp.Eq(f, 0) for f in cone_ideal)) if cone_ideal else sp.true
    return TangentConeResult(
        vars_,
        point_map,
        dirs,
        cone_ideal,
        formula,
        True,
        "saturated_m_adic_deformation",
    )


def is_singular(equations, point, variables=None, *, codimension: int | None = None) -> bool:
    """Return whether ``point`` is singular on the polynomial variety."""
    eqs = _equations(equations)
    vars_ = normalize_problem_variables(variables, sp.Tuple(*eqs))
    point_map = normalize_point(point, vars_, context=(sp.Tuple(*eqs),))
    for eq in eqs:
        if sp.simplify(eq.subs(point_map)) != 0:
            raise ValueError("point does not lie on the algebraic variety")
    locus = singular_locus(eqs, vars_, codimension=codimension)
    value = sp.simplify(locus.subs(point_map))
    if value in (sp.true, True):
        return True
    if value in (sp.false, False):
        return False
    raise ValueError("singularity test did not reduce to an exact Boolean value")


def is_smooth(equations, variables=None, *, codimension: int | None = None) -> bool:
    """Return whether the real polynomial variety has empty singular locus."""
    from .decision import is_satisfiable

    eqs = _equations(equations)
    vars_ = normalize_problem_variables(variables, sp.Tuple(*eqs))
    locus = singular_locus(eqs, vars_, codimension=codimension)
    return not is_satisfiable(locus, vars_)


def tangent_dimension(equations, point, variables=None) -> int:
    """Return the exact Zariski tangent-space dimension at ``point``."""
    return tangent_space(equations, point, variables).dimension


__all__ = [
    "singular_locus",
    "TangentSpaceResult",
    "tangent_space",
    "TangentConeResult",
    "tangent_cone",
    "is_singular",
    "is_smooth",
    "tangent_dimension",
    "AlgebraicLocusPiece",
    "DecompositionPieceCertificate",
    "DecompositionCertificate",
    "EquidimensionalDecomposition",
    "equidimensional_decomposition",
    "verify_decomposition_certificate",
]
