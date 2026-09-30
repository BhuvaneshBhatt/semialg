"""Certified local semialgebraic geometry for asymptotic analysis."""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
from math import gcd

import sympy as sp
from sympy.core.relational import Relational

from .decision import is_satisfiable
from .derived_geometry import connected_components, point_in_closure
from .geometry_queries import semialgebraic_image, semialgebraic_preimage
from .internal_symbols import collision_free_real_symbols
from .normalization import (
    normalize_formula,
    normalize_point,
    normalize_problem_variables,
    normalize_symbol_sequence,
)
from .optimization import semialgebraic_maximize, semialgebraic_minimize
from .region_analysis import local_dimension
from .regions.operations import region_closure
from .semialgebraic_local_geometry import semialgebraic_tangent_cone


@dataclass(frozen=True)
class LocalGerm:
    """A semialgebraic set viewed in arbitrarily small neighborhoods of a point."""

    formula: sp.Expr
    point: Mapping[sp.Symbol, sp.Expr]
    variables: tuple[sp.Symbol, ...]

    def representative(self, radius: sp.Expr | None = None) -> sp.Expr:
        if radius is None:
            return self.formula
        radius = sp.sympify(radius)
        if radius.is_positive is not True:
            raise ValueError("radius must be known positive")
        distance = sp.Add(*((v - self.point[v]) ** 2 for v in self.variables))
        return sp.And(self.formula, distance < radius**2)


@dataclass(frozen=True)
class LocalGeometry:
    """Dimension and tangent-cone data for a local semialgebraic set."""

    germ: LocalGerm
    dimension: int
    tangent_cone: sp.Expr


@dataclass(frozen=True)
class LocalSignStratum:
    """A local region carrying one exact sign vector."""

    formula: sp.Expr
    signs: tuple[int, ...]


@dataclass(frozen=True)
class CurveSelectionWitness:
    """A certified semialgebraic curve-selection result.

    ``mapping`` is present when an explicit polynomial/rational witness was
    constructed.  Otherwise ``existence_formula`` records the exact local
    closure statement and ``method='curve-selection-theorem'`` certifies
    existence by the semialgebraic curve-selection theorem without inventing
    an explicit parametrization.
    """

    parameter: sp.Symbol
    mapping: tuple[sp.Expr, ...] | None
    domain: sp.Expr
    target: Mapping[sp.Symbol, sp.Expr]
    variables: tuple[sp.Symbol, ...]
    region: sp.Expr
    certified: bool = True
    explicit: bool = True
    method: str = "polynomial-witness"
    existence_formula: sp.Expr | None = None


@dataclass(frozen=True)
class LocalRangeResult:
    """Exact range information on a chosen local neighborhood."""

    lower: sp.Expr
    upper: sp.Expr
    lower_attained: bool
    upper_attained: bool
    neighborhood: sp.Expr


@dataclass(frozen=True)
class LocalBoundCertificate:
    """A certified bound valid on some punctured local neighborhood."""

    bound: sp.Expr
    radius: sp.Expr
    neighborhood: sp.Expr
    condition: sp.Expr
    certified: bool = True
    method: str = "quantifier-elimination"


@dataclass(frozen=True)
class ParameterStratum:
    """A parameter condition paired with a requested invariant value."""

    condition: sp.Expr
    value: object


@dataclass(frozen=True)
class BlowupChart:
    """A polynomial weighted blow-up chart."""

    radial: sp.Symbol
    directions: tuple[sp.Symbol, ...]
    mapping: tuple[sp.Expr, ...]
    weights: tuple[int, ...]
    pivot: int


def local_germ(region, point, variables=None) -> LocalGerm:
    """Return the semialgebraic germ of ``region`` at ``point``.

    The point must lie in the closure of the region.  The returned :class:`LocalGerm`
    stores the normalized point and variable order together with the exact region
    formula; it does not replace the representative by a numerical neighborhood.
    """
    formula = normalize_formula(region)
    vars_ = normalize_problem_variables(variables, formula)
    point_map = normalize_point(point, vars_, context=(formula,))
    if not point_in_closure(formula, point_map, vars_):
        raise ValueError("point is not in the closure of the region")
    return LocalGerm(formula, point_map, tuple(vars_))


def local_geometry(region, point, variables=None) -> LocalGeometry:
    """Return local dimension and the Bouligand tangent cone at ``point``.

    Both results are computed from the exact semialgebraic germ.  Singular inputs
    can require CAD or algebraic decomposition, so this operation may be expensive.
    """
    germ = local_germ(region, point, variables)
    closed = region_closure(germ.formula, germ.variables)
    return LocalGeometry(
        germ,
        local_dimension(germ.formula, germ.point, germ.variables),
        semialgebraic_tangent_cone(closed, germ.point, germ.variables),
    )


def local_sign_strata(expressions, region, point, variables=None) -> tuple[LocalSignStratum, ...]:
    """Partition a local germ by realizable exact sign vectors.

    Sign choices are pruned incrementally by closure incidence, avoiding the
    unconditional ``3**n`` candidate materialization used by the earlier
    implementation.
    """
    germ = local_germ(region, point, variables)
    exprs = tuple(map(sp.sympify, expressions))
    active: list[tuple[sp.Expr, tuple[int, ...]]] = [(germ.formula, ())]
    for expr in exprs:
        next_active: list[tuple[sp.Expr, tuple[int, ...]]] = []
        for base, prefix in active:
            for sign, atom in ((-1, expr < 0), (0, sp.Eq(expr, 0)), (1, expr > 0)):
                candidate = sp.And(base, atom)
                if point_in_closure(candidate, germ.point, germ.variables):
                    next_active.append((candidate, (*prefix, sign)))
        active = next_active
        if not active:
            break
    return tuple(LocalSignStratum(formula, signs) for formula, signs in active)


def _factorized_plane_curve_components(germ: LocalGerm, half_width: sp.Expr):
    """Return local branches factor by factor for a reducible plane curve."""
    if len(germ.variables) != 2 or not isinstance(germ.formula, sp.Equality):
        return None
    residual = sp.expand(germ.formula.lhs - germ.formula.rhs)
    try:
        _, factors = sp.factor_list(residual, *germ.variables)
    except (sp.PolynomialError, TypeError, ValueError, NotImplementedError):
        return None
    factors = tuple(sp.expand(factor) for factor, _ in factors if factor.free_symbols)
    if len(factors) <= 1:
        return None
    offsets = tuple(v - germ.point[v] for v in germ.variables)
    box = sp.And(*(sp.And(offset > -half_width, offset < half_width) for offset in offsets))
    puncture = sp.Or(*(sp.Ne(offset, 0) for offset in offsets))
    pieces: list[sp.Expr] = []
    for factor in factors:
        branch = sp.And(sp.Eq(factor, 0), puncture, box)
        for component in connected_components(branch, germ.variables):
            if point_in_closure(component, germ.point, germ.variables, strategy="syntactic"):
                pieces.append(component)
    return tuple(pieces)


def _plane_curve_sector_components(germ: LocalGerm, half_width: sp.Expr):
    """Split a plane algebraic curve germ into open coordinate sectors.

    Strict sectors expose ordinary punctured branches of singular plane curves
    and substantially reduce CAD work.  Curves with a branch contained in a
    coordinate axis fall back to the general connectivity path.
    """
    if len(germ.variables) != 2 or not isinstance(germ.formula, sp.Equality):
        return None
    residual = sp.expand(germ.formula.lhs - germ.formula.rhs)
    try:
        poly = sp.Poly(residual, *germ.variables)
        if poly.is_zero:
            return None
        # An axis factor represents a branch omitted by strict sectors.
        if sp.rem(
            poly, sp.Poly(germ.variables[0] - germ.point[germ.variables[0]], *germ.variables)
        ).is_zero:
            return None
        if sp.rem(
            poly, sp.Poly(germ.variables[1] - germ.point[germ.variables[1]], *germ.variables)
        ).is_zero:
            return None
    except (sp.PolynomialError, TypeError, ValueError):
        return None

    offsets = tuple(v - germ.point[v] for v in germ.variables)
    box = sp.And(*(sp.And(offset > -half_width, offset < half_width) for offset in offsets))
    pieces: list[sp.Expr] = []
    for first_sign in (-1, 1):
        for second_sign in (-1, 1):
            atoms = (
                offsets[0] < 0 if first_sign < 0 else offsets[0] > 0,
                offsets[1] < 0 if second_sign < 0 else offsets[1] > 0,
            )
            sector = sp.And(germ.formula, box, *atoms)
            if not is_satisfiable(sector, germ.variables):
                continue
            for component in connected_components(sector, germ.variables):
                if point_in_closure(component, germ.point, germ.variables, strategy="syntactic"):
                    pieces.append(component)
    return tuple(pieces) if pieces else None


def local_components(region, point, variables=None) -> tuple[LocalGerm, ...]:
    """Return connected branches of the punctured germ incident at ``point``.

    Connectivity is computed after restricting to an exact sufficiently small
    neighborhood derived from critical-distance information.  This distinguishes
    branches that reconnect only away from the base point.  Singular multivariate
    cases can require expensive CAD and algebraic sign certification.
    """
    germ = local_germ(region, point, variables)
    puncture = sp.Or(*(sp.Ne(v, germ.point[v]) for v in germ.variables))
    punctured = sp.And(germ.formula, puncture)
    if len(germ.variables) == 1:
        variable = germ.variables[0]
        value = germ.point[variable]
        candidates = (sp.And(punctured, variable < value), sp.And(punctured, variable > value))
        pieces = tuple(
            candidate for candidate in candidates if is_satisfiable(candidate, germ.variables)
        )
    else:
        from .exact_arithmetic import compare_exact_reals
        from .geometry_queries import critical_values

        atoms = germ.formula.args if isinstance(germ.formula, sp.And) else (germ.formula,)
        linear = True
        for atom in atoms:
            if not isinstance(atom, Relational):
                linear = False
                break
            try:
                if sp.Poly(atom.lhs - atom.rhs, *germ.variables).total_degree() > 1:
                    linear = False
                    break
            except (sp.PolynomialError, TypeError, ValueError):
                linear = False
                break
        if linear:
            half_width = sp.Rational(1, 2)
            offsets = tuple(v - germ.point[v] for v in germ.variables)
            box = sp.And(*(sp.And(offset > -half_width, offset < half_width) for offset in offsets))
            pieces = connected_components(sp.And(punctured, box), germ.variables)
            return tuple(
                LocalGerm(piece, germ.point, germ.variables)
                for piece in pieces
                if point_in_closure(piece, germ.point, germ.variables, strategy="syntactic")
            )

        distance = sp.Add(*((v - germ.point[v]) ** 2 for v in germ.variables))
        dnf = sp.to_dnf(germ.formula, simplify=False, force=True)
        clauses = dnf.args if isinstance(dnf, sp.Or) else (dnf,)
        values = tuple(
            value
            for clause in clauses
            for value in critical_values(distance, clause, germ.variables)
        )
        positive = []
        for value in values:
            try:
                if compare_exact_reals(value, 0) > 0:
                    positive.append(value)
            except (TypeError, ValueError, NotImplementedError):
                continue
        if positive:
            threshold = positive[0]
            for value in positive[1:]:
                if compare_exact_reals(value, threshold) < 0:
                    threshold = value
            from .algebraic.roots import rational_between_algebraic_reals
            from .algebraic.samples import RationalSample

            radius_sq = rational_between_algebraic_reals(RationalSample(sp.Integer(0)), threshold)
        else:
            radius_sq = sp.Integer(1)
        # Use a rational axis-aligned neighborhood.  A box is locally
        # equivalent to a ball and avoids introducing square-root boundary
        # sections into CAD for otherwise rational input.
        half_width = sp.Rational(1, 2)
        if radius_sq.is_Rational and radius_sq > 0:
            half_width = min(half_width, radius_sq)
        box = sp.And(
            *(
                sp.And(v - germ.point[v] > -half_width, v - germ.point[v] < half_width)
                for v in germ.variables
            )
        )
        factor_pieces = _factorized_plane_curve_components(germ, half_width)
        if factor_pieces is not None:
            return tuple(LocalGerm(piece, germ.point, germ.variables) for piece in factor_pieces)
        sector_pieces = _plane_curve_sector_components(germ, half_width)
        if sector_pieces is not None:
            return tuple(LocalGerm(piece, germ.point, germ.variables) for piece in sector_pieces)
        local_region = sp.And(punctured, box)
        pieces = connected_components(local_region, germ.variables)
    return tuple(
        LocalGerm(piece, germ.point, germ.variables)
        for piece in pieces
        if point_in_closure(piece, germ.point, germ.variables)
    )


def _monomial_curve_candidates(germ: LocalGerm, max_order: int):
    t = sp.Dummy("t", positive=True)
    offsets = (-1, 1)
    for order in range(1, max_order + 1):
        for axis in range(len(germ.variables)):
            for sign in offsets:
                mapping = tuple(
                    germ.point[v] + (sign * t**order if i == axis else 0)
                    for i, v in enumerate(germ.variables)
                )
                yield t, mapping
    if len(germ.variables) == 2:
        for p in range(1, max_order + 1):
            for q in range(1, max_order + 1):
                for sx in offsets:
                    for sy in offsets:
                        yield (
                            t,
                            (
                                germ.point[germ.variables[0]] + sx * t**p,
                                germ.point[germ.variables[1]] + sy * t**q,
                            ),
                        )


def _certified_curve_domain(parameter: sp.Symbol, pulled: sp.Expr) -> sp.Expr | None:
    """Return a certified punctured interval on which a candidate stays in the set."""
    unit = sp.And(parameter > 0, parameter < 1)
    if not is_satisfiable(sp.And(unit, sp.Not(pulled)), (parameter,)):
        return unit

    from .qe import quantifier_eliminate
    from .solve import find_instance

    epsilon = sp.Dummy("epsilon", real=True)
    matrix = sp.Implies(sp.And(parameter > 0, parameter < epsilon), pulled)
    condition = quantifier_eliminate(
        matrix, (("forall", parameter),), variables=(epsilon, parameter), strategy="complete"
    )
    witness = find_instance(sp.And(epsilon > 0, condition), (epsilon,), return_result=True)
    if not witness.found or witness.first() is None:
        return None
    value = sp.sympify(witness.first()[epsilon])
    if value.is_positive is not True and sp.simplify(value > 0) is not sp.true:
        return None
    return sp.And(parameter > 0, parameter < value)


def curve_selection(region, point, variables=None, *, max_order: int = 4) -> CurveSelectionWitness:
    """Return a certified semialgebraic curve-selection result.

    Cheap polynomial witnesses are tried first.  If they do not suffice, the
    exact closure test certifies existence using the semialgebraic
    curve-selection theorem.  The fallback is explicit about not having an
    explicit parametrization, so downstream code cannot mistake existence for
    a usable pullback curve.
    """
    germ = local_germ(region, point, variables)
    for parameter, mapping in _monomial_curve_candidates(germ, max_order):
        pulled = sp.simplify(germ.formula.xreplace(dict(zip(germ.variables, mapping, strict=True))))
        domain = _certified_curve_domain(parameter, pulled)
        if domain is not None:
            return CurveSelectionWitness(
                parameter, mapping, domain, germ.point, germ.variables, germ.formula
            )

    parameter = sp.Dummy("t", positive=True)
    distance = sp.Add(*((v - germ.point[v]) ** 2 for v in germ.variables))
    from .quantifiers import Exists, ForAll

    nearby = sp.And(germ.formula, distance > 0, distance < parameter**2)
    existence = ForAll(
        parameter,
        sp.Implies(
            sp.And(parameter > 0, parameter < 1),
            Exists(germ.variables, nearby),
        ),
    )
    # local_germ already certified the closure condition.  For a semialgebraic
    # set this is exactly the hypothesis of the curve-selection theorem.
    return CurveSelectionWitness(
        parameter,
        None,
        sp.And(parameter > 0, parameter < 1),
        germ.point,
        germ.variables,
        germ.formula,
        explicit=False,
        method="curve-selection-theorem",
        existence_formula=existence,
    )


def local_range(expression, region, point, variables=None, *, radius=1) -> LocalRangeResult:
    """Return the exact range on a fixed local ball around the germ point.

    ``radius`` defines the neighborhood being certified.  The result reports exact
    extrema for the intersection of the germ representative with that ball; it is
    not an asymptotic estimate.
    """
    germ = local_germ(region, point, variables)
    radius = sp.sympify(radius)
    if radius.is_positive is not True:
        raise ValueError("radius must be known positive")
    distance = sp.Add(*((v - germ.point[v]) ** 2 for v in germ.variables))
    neighborhood = sp.And(germ.formula, distance <= radius**2)
    lo = semialgebraic_minimize(expression, neighborhood, germ.variables, return_result=True)
    hi = semialgebraic_maximize(expression, neighborhood, germ.variables, return_result=True)
    return LocalRangeResult(lo.value, hi.value, lo.attained, hi.attained, neighborhood)


def local_bound(
    expression,
    region,
    point,
    variables=None,
    *,
    radius=None,
    bound=None,
    return_certificate: bool = False,
):
    """Certify a local absolute bound, optionally finding a valid radius.

    With ``radius`` supplied this preserves the fixed-ball behavior.  With
    ``radius=None`` an exact quantified condition is solved for a positive
    radius.  ``bound`` may be supplied to test a requested bound; otherwise a
    finite bound and radius are sought together.
    """
    if radius is not None:
        result = local_range(
            sp.Abs(sp.sympify(expression)), region, point, variables, radius=radius
        )
        if return_certificate:
            return LocalBoundCertificate(
                result.upper, sp.sympify(radius), result.neighborhood, sp.true
            )
        return result.upper

    from .qe import quantifier_eliminate
    from .quantifiers import ForAll
    from .solve import find_instance

    germ = local_germ(region, point, variables)
    expr = sp.sympify(expression)
    r = sp.Dummy("radius", real=True)
    b = sp.Dummy("bound", real=True)
    chosen_bound = sp.sympify(bound) if bound is not None else b
    if bound is not None and chosen_bound.is_nonnegative is False:
        raise ValueError("bound must be nonnegative")
    distance = sp.Add(*((v - germ.point[v]) ** 2 for v in germ.variables))
    punctured = sp.And(germ.formula, distance > 0, distance < r**2)
    matrix = sp.Implies(punctured, sp.And(expr <= chosen_bound, expr >= -chosen_bound))
    condition = quantifier_eliminate(
        ForAll(germ.variables, matrix),
        variables=(r, b) if bound is None else (r,),
        strategy="complete",
    )
    search = sp.And(r > 0, b >= 0, condition) if bound is None else sp.And(r > 0, condition)
    search_vars = (r, b) if bound is None else (r,)
    witness = find_instance(search, search_vars, strict=True, return_result=True)
    if not witness.found or witness.first() is None:
        raise ValueError("no finite local bound could be certified")
    values = witness.first()
    radius_value = sp.sympify(values[r])
    bound_value = chosen_bound if bound is not None else sp.sympify(values[b])

    def violates(candidate_radius: sp.Expr) -> bool:
        candidate_region = sp.And(germ.formula, distance > 0, distance < candidate_radius**2)
        violation = sp.And(
            candidate_region,
            sp.Or(expr > bound_value, expr < -bound_value),
        )
        return is_satisfiable(violation, germ.variables, strategy="cad")

    # Never trust a parametric QE witness without checking the concrete local
    # claim independently.  If it lands on a boundary, shrink through exact
    # rational radii; local validity is monotone under shrinking.
    if violates(radius_value):
        found = None
        for exponent in range(1, 65):
            candidate = sp.Rational(1, 2**exponent)
            if not violates(candidate):
                found = candidate
                break
        if found is None:
            raise ValueError("no finite local bound could be independently certified")
        radius_value = found

    neighborhood = sp.And(germ.formula, distance > 0, distance < radius_value**2)
    certificate = LocalBoundCertificate(bound_value, radius_value, neighborhood, condition)
    return certificate if return_certificate else bound_value


def parameter_strata(
    formula,
    variables,
    parameters,
    *,
    conditions_by_value: Mapping[object, sp.Expr] | None = None,
) -> tuple[ParameterStratum, ...]:
    """Partition parameter space by feasibility or caller-defined invariants.

    Without ``conditions_by_value`` the two strata classify feasibility of
    ``formula``.  With it, each value is associated with an additional
    semialgebraic condition in the quantified variables and parameters.
    Conditions are made disjoint in insertion order; callers should order
    overlapping invariant cases from most specific to least specific.
    """
    from .parameters import solvability_conditions

    base = normalize_formula(formula)
    if conditions_by_value is None:
        condition = solvability_conditions(base, variables, parameters)
        return (
            ParameterStratum(condition, True),
            ParameterStratum(sp.Not(condition), False),
        )

    strata: list[ParameterStratum] = []
    covered = sp.false
    for value, invariant_condition in conditions_by_value.items():
        condition = solvability_conditions(
            sp.And(base, normalize_formula(invariant_condition)),
            variables,
            parameters,
        )
        exclusive = sp.And(condition, sp.Not(covered))
        if exclusive is not sp.false:
            strata.append(ParameterStratum(exclusive, value))
        covered = sp.Or(covered, condition)
    return tuple(strata)


def vanishing_order(expression, curve: CurveSelectionWitness) -> sp.Expr:
    """Return the exact vanishing order along an explicit curve witness.

    The witness must contain a usable parametrization.  A theorem-only
    curve-selection certificate proves existence but cannot define this pullback.
    """
    if not curve.explicit or curve.mapping is None:
        raise ValueError("vanishing_order requires an explicit curve-selection witness")
    expr = sp.sympify(expression).xreplace(dict(zip(curve.variables, curve.mapping, strict=True)))
    t = curve.parameter
    numerator, denominator = sp.fraction(sp.cancel(expr))
    try:
        num_poly = sp.Poly(numerator, t)
        den_poly = sp.Poly(denominator, t)
    except sp.PolynomialError as exc:
        raise NotImplementedError(
            "vanishing order currently requires a rational curve pullback"
        ) from exc

    def valuation(poly):
        if poly.is_zero:
            return sp.oo
        return min(monomial[0] for monomial, _ in poly.terms())

    return valuation(num_poly) - valuation(den_poly)


def contact_order(left, right, curve: CurveSelectionWitness) -> sp.Expr:
    """Return the exact contact order of two expressions along a curve.

    This is the vanishing order of ``left - right`` after pullback to the explicit
    curve witness.
    """
    return vanishing_order(sp.sympify(left) - sp.sympify(right), curve)


def blowup_charts(variables, point=None, *, weights=None) -> tuple[BlowupChart, ...]:
    """Return polynomial blow-up charts centered at ``point``.

    Unit weights give ordinary blow-up charts; positive integer ``weights`` expose
    anisotropic approach scales.  The charts describe exact polynomial coordinate
    maps and do not themselves perform an asymptotic expansion.
    """
    vars_ = tuple(variables)
    center = tuple(sp.sympify(v) for v in (point if point is not None else (0,) * len(vars_)))
    if len(center) != len(vars_):
        raise ValueError("point dimension must match variables")
    ws = tuple(int(w) for w in (weights if weights is not None else (1,) * len(vars_)))
    if len(ws) != len(vars_) or any(w <= 0 for w in ws):
        raise ValueError("weights must be positive integers matching variables")
    common = 0
    for w in ws:
        common = gcd(common, w)
    if common > 1:
        ws = tuple(w // common for w in ws)
    charts = []
    for pivot in range(len(vars_)):
        radial = sp.Dummy(f"r{pivot}", real=True)
        directions = tuple(sp.Dummy(f"u{pivot}_{i}", real=True) for i in range(len(vars_) - 1))
        it = iter(directions)
        mapping = []
        for i, (c, w) in enumerate(zip(center, ws, strict=True)):
            mapping.append(c + radial**w if i == pivot else c + radial**w * next(it))
        charts.append(BlowupChart(radial, directions, tuple(mapping), ws, pivot))
    return tuple(charts)


def local_image(mapping, germ: LocalGerm, *, image_variables=None) -> LocalGerm:
    """Return the exact image germ of ``germ`` under ``mapping``.

    The image is obtained by exact semialgebraic graph projection.  General maps can
    therefore inherit the cost and supported-fragment limits of quantifier
    elimination.
    """
    maps = tuple(map(sp.sympify, mapping))
    targets = (
        collision_free_real_symbols("y", len(maps), germ.formula, maps, germ.variables)
        if image_variables is None
        else normalize_symbol_sequence(image_variables)
    )
    formula = semialgebraic_image(maps, germ.formula, germ.variables, image_variables=targets)
    target_point = {
        y: sp.simplify(expr.subs(germ.point)) for y, expr in zip(targets, maps, strict=True)
    }
    return local_germ(formula, target_point, targets)


def local_preimage(target: LocalGerm, mapping, variables) -> LocalGerm:
    """Return the preimage germ when the target-point fiber has one source point.

    Existence and uniqueness of the real source point are certified exactly.  An
    empty, multiple, or positive-dimensional fiber is rejected instead of choosing
    an arbitrary representative.
    """
    source = tuple(variables)
    maps = tuple(map(sp.sympify, mapping))
    formula = semialgebraic_preimage(
        maps, target.formula, source, target_variables=target.variables
    )
    equations = tuple(
        sp.Eq(expr, target.point[y]) for expr, y in zip(maps, target.variables, strict=True)
    )
    fiber = sp.And(*equations)
    from .solve import find_instance

    witness = find_instance(fiber, source, return_result=True)
    point = witness.first()
    if not witness.found or point is None or any(v not in point for v in source):
        raise ValueError("target germ point has no exact source point")
    different = sp.Or(*(sp.Ne(v, point[v]) for v in source))
    if is_satisfiable(sp.And(fiber, different), source):
        raise ValueError("local_preimage requires a unique exact source point")
    return local_germ(formula, point, source)


def path_independent(expression, germ: LocalGerm, value) -> bool:
    """Certify local path independence of ``expression`` at the germ point.

    The equality to ``value`` need only hold on some sufficiently small punctured
    neighborhood of the germ; it is not required on the entire supplied
    representative.
    """
    expr = sp.sympify(expression) - sp.sympify(value)
    try:
        local_bound(
            expr,
            germ.formula,
            germ.point,
            germ.variables,
            bound=0,
            return_certificate=True,
        )
    except ValueError:
        return False
    return True


__all__ = [
    "BlowupChart",
    "CurveSelectionWitness",
    "LocalGeometry",
    "LocalGerm",
    "LocalRangeResult",
    "LocalBoundCertificate",
    "LocalSignStratum",
    "ParameterStratum",
    "blowup_charts",
    "contact_order",
    "curve_selection",
    "local_bound",
    "local_components",
    "local_geometry",
    "local_germ",
    "local_image",
    "local_preimage",
    "local_range",
    "local_sign_strata",
    "parameter_strata",
    "path_independent",
    "vanishing_order",
]
