"""Bounded, cacheable geometry certificates for coordinate and coarea charts.

This module proves geometry only. Null sets do not certify integrability.
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
from itertools import permutations
from math import comb, factorial

import sympy as sp
from sympy.core.relational import Relational
from sympy.logic.algorithms.lra_theory import UnhandledInput

__all__ = [
    "GeometryBudget",
    "GeometryProofError",
    "NullSetCertificate",
    "Chart",
    "ChartCertificate",
    "CoverCertificate",
    "certify_null_set",
    "certify_chart",
    "certify_cover",
    "polar_chart",
    "clear_chart_cache",
    "chart_cache_info",
    "LinearMetricCertificate",
    "certify_linear_metric",
    "chart_parameter_events",
    "parameter_geometry_events",
]


class GeometryProofError(NotImplementedError):
    """A geometric prerequisite could not be proved within the supplied budget."""


@dataclass(frozen=True)
class GeometryBudget:
    max_variables: int = 4
    max_degree: int = 8
    max_atoms: int = 24
    max_charts: int = 64
    max_operations: int = 1000

    def __post_init__(self):
        if any(type(v) is not int or v < 0 for v in vars(self).values()):
            raise ValueError("Geometry budget entries must be nonnegative integers")


DEFAULT_GEOMETRY_BUDGET = GeometryBudget()
_REAL_COORDINATE_FACTS = sp.Symbol("geometry_coordinate", real=True).assumptions0


@dataclass(frozen=True)
class NullSetCertificate:
    region: sp.Expr
    variables: tuple
    assumptions: tuple
    zero_polynomials: tuple
    basis: str = "subset of nonzero polynomial zero sets in each parameter fiber"


@dataclass(frozen=True)
class Chart:
    coordinates: tuple
    forward: tuple
    inverse: tuple
    ranges: tuple
    chart_exclusions: sp.Expr = sp.false
    source_exclusions: sp.Expr = sp.false
    source_domain: sp.Expr = sp.true
    kind: str = "rational"

    def __post_init__(self):
        for name in ("coordinates", "forward", "inverse"):
            object.__setattr__(self, name, tuple(map(sp.sympify, getattr(self, name))))
        object.__setattr__(self, "ranges", tuple(tuple(map(sp.sympify, r)) for r in self.ranges))
        for name in ("chart_exclusions", "source_exclusions", "source_domain"):
            object.__setattr__(self, name, sp.sympify(getattr(self, name)))


@dataclass(frozen=True)
class ChartCertificate:
    chart: Chart
    variables: tuple
    level_map: sp.Expr
    assumptions: tuple
    domain: sp.Expr
    image: sp.Expr
    interior_image: sp.Expr
    determinant: sp.Expr
    map_jacobian_squared: sp.Expr
    fiber_metric_squared: sp.Expr
    chart_null: NullSetCertificate
    source_null: NullSetCertificate
    basis: str = "inverse, regular rank and coarea metric identities"


@dataclass(frozen=True)
class CoverCertificate:
    region: sp.Expr
    variables: tuple
    assumptions: tuple
    charts: tuple
    source_null: NullSetCertificate
    basis: str = "coverage modulo certified null set and disjoint interiors"


def _assumptions(assumptions):
    if assumptions is None:
        return ()
    values = assumptions if isinstance(assumptions, (tuple, list)) else (assumptions,)
    return tuple(sp.sympify(a) for a in values)


def _local_ask(predicate, assumptions=sp.true, *, context=None):
    """Isolate every symbol from global facts, including recursive ask handlers.

    SymPy's compound-predicate handlers may omit the explicit context during
    recursion. Fresh symbols retain intrinsic facts and explicit cell assumptions
    while preventing a temporary global fact from entering a cached certificate.
    """
    predicate, assumptions = map(sp.sympify, (predicate, assumptions))
    replacements = {
        v: sp.Dummy("certificate_scope", **v.assumptions0)
        for v in predicate.free_symbols | assumptions.free_symbols
    }
    return sp.ask(
        predicate.xreplace(replacements),
        assumptions.xreplace(replacements),
        context=set() if context is None else context,
    )


def _predicates(condition):
    if isinstance(condition, Relational):
        e = condition.lhs - condition.rhs
        if condition.rel_op == "!=" and e.is_real is not True:
            return ~sp.Q.zero(e)
        return {
            ">": sp.Q.positive,
            ">=": sp.Q.nonnegative,
            "<": sp.Q.negative,
            "<=": sp.Q.nonpositive,
            "==": sp.Q.zero,
            "!=": sp.Q.nonzero,
        }[condition.rel_op](e)
    if condition.func in (sp.And, sp.Or, sp.Not):
        return condition.func(*map(_predicates, condition.args))
    return condition


def _bounded(expressions, variables, budget):
    exprs = tuple(map(sp.sympify, expressions))
    if sum(sp.count_ops(e) for e in exprs) > budget.max_operations:
        raise GeometryProofError("Geometry operation budget exceeded")
    for e in exprs:
        for power in e.atoms(sp.Pow):
            if power.exp.is_Integer and power.base.free_symbols:
                exponent = abs(int(power.exp))
                if exponent > budget.max_degree:
                    raise GeometryProofError(
                        "Geometry polynomial degree budget exceeded before expansion"
                    )
                if (
                    power.base.is_Add
                    and comb(exponent + len(power.base.args) - 1, len(power.base.args) - 1)
                    > budget.max_operations
                ):
                    raise GeometryProofError("Geometry expansion budget exceeded")
    if len(variables) > budget.max_variables:
        raise GeometryProofError("Geometry variable budget exceeded")
    if len(set().union(*(e.atoms(Relational) for e in exprs))) > budget.max_atoms:
        raise GeometryProofError("Geometry Boolean atom budget exceeded")
    for e in exprs:
        relations = e.atoms(Relational)
        predicates = e.atoms(sp.AppliedPredicate)
        terms = (
            ([a.lhs - a.rhs for a in relations] + [arg for a in predicates for arg in a.arguments])
            if relations or predicates or e.is_Boolean or e in (sp.true, sp.false)
            else [e]
        )
        for term in terms:
            for part in sp.cancel(term).as_numer_denom():
                if not variables:
                    continue
                try:
                    poly = sp.Poly(part, *variables)
                except sp.PolynomialError:
                    continue
                if poly.total_degree() > budget.max_degree:
                    raise GeometryProofError("Geometry polynomial degree budget exceeded")


def _determinant(matrix, budget):
    """Reject excessive symbolic product growth before determinant expansion."""
    if matrix.rows != matrix.cols or matrix.rows > budget.max_variables:
        raise GeometryProofError("Determinant dimension budget exceeded")
    if sum(sp.count_ops(e) for e in matrix) > budget.max_operations:
        raise GeometryProofError("Determinant operation budget exceeded")

    if matrix.is_upper or matrix.is_lower:
        return sp.prod(matrix[i, i] for i in range(matrix.rows))
    if factorial(matrix.rows) > budget.max_operations:
        raise GeometryProofError("Determinant expansion budget exceeded")

    def terms(e):
        if e == 0:
            return 0
        if e.is_Add:
            return sum(terms(a) for a in e.args)
        if e.is_Mul:
            count = 1
            for a in e.args:
                count *= terms(a)
                if count > budget.max_operations:
                    break
            return count
        if e.is_Pow and e.exp.is_Integer:
            count = terms(e.base)
            if count > 1:
                exponent = abs(int(e.exp))
                if exponent > budget.max_degree:
                    return budget.max_operations + 1
                return comb(count + exponent - 1, exponent)
        return 1

    estimate = 0
    for order in permutations(range(matrix.rows)):
        count = 1
        for row, col in enumerate(order):
            count *= terms(matrix[row, col])
        estimate += count
        if estimate > budget.max_operations:
            raise GeometryProofError("Determinant expansion budget exceeded")
    return matrix.det()


def _positive_constraints(premise):
    for atom in sp.And.make_args(premise):
        if isinstance(atom, Relational) and atom.rel_op in (">", ">=", "<", "<="):
            yield (
                (atom.lhs - atom.rhs if atom.rel_op in (">", ">=") else atom.rhs - atom.lhs),
                atom.rel_op in (">", "<"),
            )
        elif isinstance(atom, sp.AppliedPredicate):
            if atom.function in (sp.Q.positive, sp.Q.nonnegative):
                yield atom.arguments[0], atom.function == sp.Q.positive
            elif atom.function in (sp.Q.negative, sp.Q.nonpositive):
                yield -atom.arguments[0], atom.function == sp.Q.negative


def _positive(expr, premise, strict=False, depth=0):
    if depth > 8:
        return False
    expr = sp.factor(expr)
    if (
        strict
        and sp.Ne(expr, 0) in sp.And.make_args(premise)
        and _positive(expr, premise, False, depth + 1)
    ):
        return True
    if (expr.is_positive if strict else expr.is_nonnegative) is True:
        return True
    for d, atom_strict in _positive_constraints(premise):
        if d != 0:
            ratio = sp.cancel(expr / d)
            if ratio.is_positive is True and not ratio.free_symbols and (not strict or atom_strict):
                return True
            residual = sp.expand(expr - d)
            if (
                residual.is_positive if strict and not atom_strict else residual.is_nonnegative
            ) is True:
                return True
    if expr.is_Add:
        if all(_positive(a, premise, False, depth + 1) for a in expr.args):
            return not strict or any(_positive(a, premise, True, depth + 1) for a in expr.args)
    if expr.is_Mul:
        signs = []
        for a in expr.args:
            if _positive(a, premise, strict, depth + 1):
                signs.append(1)
            elif _positive(-a, premise, strict, depth + 1):
                signs.append(-1)
            else:
                break
        if len(signs) == len(expr.args) and signs.count(-1) % 2 == 0:
            return True
    if expr.is_Pow and expr.exp.is_even is True and _native(sp.Ne(expr.base, 0), premise):
        return True
    if expr.is_Pow and expr.exp.is_real is True and _positive(expr.base, premise, True, depth + 1):
        return True
    if not strict:
        try:
            return _local_ask(sp.Q.nonnegative(expr), _predicates(premise), context=set()) is True
        except (ValueError, TypeError, AssertionError, UnhandledInput):
            pass
    return False


def _native(conclusion, premise):
    if conclusion is sp.true or conclusion in sp.And.make_args(premise):
        return True
    if conclusion.func == sp.And:
        return all(_native(a, premise) for a in conclusion.args)
    if conclusion.func == sp.Or:
        return any(_native(a, premise) for a in conclusion.args)
    if isinstance(conclusion, Relational):
        e = sp.factor(conclusion.lhs - conclusion.rhs)
        op = conclusion.rel_op
        if op == "==":
            return e == 0
        if op == "!=":
            if sp.Q.nonzero(e) in sp.And.make_args(premise):
                return True
            if e.is_zero is False:
                return True
            if _positive(e, premise, True) or _positive(-e, premise, True):
                return True
            if e.is_Mul:
                return all(_native(sp.Ne(a, 0), premise) for a in e.args)
            if e.is_Pow:
                return _native(sp.Ne(e.base, 0), premise)
        if op in (">", ">=", "<", "<="):
            return _positive(e if op in (">", ">=") else -e, premise, op in (">", "<"))
    return False


def _prove(premise, conclusion, variables, budget):
    if sp.simplify(premise) is sp.false:
        return True
    conclusion = sp.simplify(conclusion)
    if premise is sp.false or _native(conclusion, premise):
        return True
    if conclusion.func == sp.And:
        return all(_prove(premise, a, variables, budget) for a in conclusion.args)

    # Clear rational relations only after their denominator's sign is proved.
    def clear(e):
        if isinstance(e, Relational):
            n, d = sp.cancel(e.lhs - e.rhs).as_numer_denom()
            if _positive(d, premise, True):
                return e.func(n, 0)
            if _positive(-d, premise, True):
                return e.func(-n, 0)
        if e.func in (sp.And, sp.Or, sp.Not):
            return e.func(*map(clear, e.args))
        return e

    premise, conclusion = clear(premise), clear(conclusion)
    universe = tuple(
        sorted(premise.free_symbols | conclusion.free_symbols, key=sp.default_sort_key)
    )
    _bounded((premise, conclusion), universe, budget)
    query = sp.And(premise, sp.Not(conclusion))
    try:
        if sp.satisfiable(query, use_lra_theory=True) is False:
            return True
    except (ValueError, AssertionError, TypeError, sp.PolynomialError, UnhandledInput):
        pass
    if any(not isinstance(a, Relational) for a in query.atoms(sp.AppliedPredicate)):
        return False
    if any(not (a.lhs - a.rhs).is_polynomial(*universe) for a in query.atoms(Relational)):
        return False
    from .decision import implies

    try:
        return implies(premise, conclusion, universe) is True
    except (NotImplementedError, ValueError, ArithmeticError):
        return False


def _domain(ranges, interior=False):
    return sp.And(
        *(sp.And(v > lo, v < hi) if interior else sp.And(v >= lo, v <= hi) for v, lo, hi in ranges)
    )


def certify_null_set(region, variables, *, assumptions=None, budget=DEFAULT_GEOMETRY_BUDGET):
    return _null_cached(sp.sympify(region), tuple(variables), _assumptions(assumptions), budget)


@lru_cache(maxsize=256)
def _null_cached(region, variables, assumptions, budget):
    if (
        not variables
        or len(set(variables)) != len(variables)
        or any(not isinstance(v, sp.Symbol) or v.is_real is not True for v in variables)
    ):
        raise ValueError("Null-set certification requires distinct real ambient coordinates")
    premise = sp.And(*assumptions)
    _bounded((region, premise), variables, budget)
    if region is sp.false:
        return NullSetCertificate(region, variables, assumptions, ())
    polynomials = []
    for a in region.atoms(sp.Equality):
        try:
            p = sp.Poly(a.lhs - a.rhs, *variables)
        except sp.PolynomialError:
            continue
        try:
            finite_real = all(
                _local_ask(sp.Q.real(c) & sp.Q.finite(c), _predicates(premise), context=set())
                is True
                for c in p.coeffs()
            )
        except (ValueError, AssertionError, UnhandledInput):
            finite_real = False
        if not finite_real:
            continue
        # An ambient polynomial must remain nonzero on every allowed parameter fiber.
        if not p.is_zero and any(_native(sp.Ne(c, 0), premise) for c in p.coeffs()):
            polynomials.append(p.as_expr())
    zeros = sp.Or(*(sp.Eq(p, 0) for p in polynomials))
    if not polynomials or not _prove(sp.And(premise, region), zeros, variables, budget):
        raise GeometryProofError(
            "Exceptional set is not proved ambient null on every parameter fiber"
        )
    return NullSetCertificate(
        region, variables, assumptions, tuple(sorted(set(polynomials), key=sp.default_sort_key))
    )


def _identity(expr, premise):
    value = sp.simplify(expr)
    replacements = {}
    for node in value.atoms(sp.Abs):
        arg = node.args[0]
        if _native(arg >= 0, premise):
            replacements[node] = arg
        elif _native(arg <= 0, premise):
            replacements[node] = -arg
    return sp.simplify(value.xreplace(replacements)) == 0


def _formula(expr, variables, domain, assumptions, kind, budget):
    expr = sp.sympify(expr)
    _bounded((expr,), variables, budget)
    if kind == "algebraic":
        roots = [a for a in expr.atoms(sp.Pow) if a.exp.is_Rational and a.exp.q > 1]
        for a in roots:
            if not _prove(sp.And(*assumptions, domain), a.base > 0, variables, budget):
                raise GeometryProofError("Algebraic root branch positivity is unproved")
        substituted = expr.xreplace({a: sp.Dummy(real=True) for a in roots})
    else:
        substituted = expr
    n, d = sp.cancel(substituted).as_numer_denom()
    if not n.is_polynomial(*variables) or not d.is_polynomial(*variables):
        raise GeometryProofError("Chart requires rational or selected positive algebraic formulas")
    if not _prove(
        sp.And(*assumptions, domain), sp.Ne(sp.denom(sp.cancel(expr)), 0), variables, budget
    ):
        raise GeometryProofError("Chart denominator exclusion is unproved")
    for p in (n, d):
        for c in sp.Poly(p, *variables).coeffs():
            if c.has(sp.Dummy):
                continue
            try:
                real = _local_ask(
                    sp.Q.real(c) & sp.Q.finite(c), _predicates(sp.And(*assumptions)), context=set()
                )
            except (ValueError, AssertionError):
                real = None
            if real is not True:
                raise GeometryProofError("Chart coefficients must be finite and real")


def certify_chart(chart, level_map, variables, *, assumptions=None, budget=DEFAULT_GEOMETRY_BUDGET):
    if not isinstance(chart, Chart):
        raise ValueError("Expected a Chart candidate")
    return _chart_cached(
        chart,
        sp.Tuple(*level_map)
        if isinstance(level_map, (tuple, list, sp.Tuple))
        else sp.sympify(level_map),
        tuple(variables),
        _assumptions(assumptions),
        budget,
    )


@lru_cache(maxsize=128)
def _chart_cached(chart, g, variables, assumptions, budget):
    coords = chart.coordinates
    maps = tuple(g) if isinstance(g, sp.Tuple) else (g,)
    if not 1 <= len(maps) < len(variables):
        raise ValueError("Coarea maps require 1 <= level dimension < source dimension")
    if (
        len(coords) != len(variables)
        or len(set(coords)) != len(coords)
        or set(coords) & set(variables)
        or any(not isinstance(v, sp.Symbol) or v.is_real is not True for v in (*coords, *variables))
    ):
        raise ValueError("Use distinct real chart and source coordinates of equal dimensions")
    if any(
        any(_REAL_COORDINATE_FACTS.get(k) != value for k, value in v.assumptions0.items())
        for v in (*coords, *variables)
    ):
        raise GeometryProofError(
            "Use plain real coordinate symbols; put sign and arithmetic restrictions in the domain"
        )
    if len(chart.forward) != len(coords) or len(chart.inverse) != len(coords):
        raise ValueError("Chart dimensions do not match")
    if tuple(r[0] for r in chart.ranges) != tuple(reversed(coords)):
        raise ValueError("Ranges must be inner-first, level coordinate last")
    external = (
        g.free_symbols
        | sp.And(*assumptions).free_symbols
        | chart.source_domain.free_symbols
        | chart.source_exclusions.free_symbols
    ) - set(variables)
    if set(coords) & external:
        raise ValueError("Chart coordinates cannot capture parameters")
    if sp.And(*assumptions).has(*variables, *coords):
        raise ValueError("Chart assumptions concern parameters only")
    _bounded((*chart.forward, _domain(chart.ranges)), coords, budget)
    _bounded((*chart.inverse, *maps, chart.source_domain), variables, budget)
    for i, (_, lo, hi) in enumerate(chart.ranges):
        if (lo.free_symbols | hi.free_symbols) & {r[0] for r in chart.ranges[: i + 1]}:
            raise ValueError("Bounds depend only on outer coordinates")
        outer = _domain(chart.ranges[i + 1 :])
        for b in (lo, hi):
            if b in (-sp.oo, sp.oo):
                continue
            try:
                real = _local_ask(
                    sp.Q.real(b) & sp.Q.finite(b),
                    _predicates(sp.And(*assumptions, outer)),
                    context=set(),
                )
            except (ValueError, AssertionError, UnhandledInput):
                real = None
            if real is not True:
                raise GeometryProofError(
                    "Chart interval bounds must be real and finite, except infinities"
                )
        inside = _domain(chart.ranges[i + 1 :], True)
        if not _prove(sp.And(*assumptions, outer), hi >= lo, coords, budget) or not _prove(
            sp.And(*assumptions, inside), hi > lo, coords, budget
        ):
            raise GeometryProofError("Chart intervals need proved positive width")
    full = _domain(chart.ranges)
    regular = sp.And(full, sp.Not(chart.chart_exclusions))
    premise = sp.And(*assumptions, regular)
    _bounded((*chart.forward, *chart.inverse, full, *maps), coords, budget)
    cn = certify_null_set(chart.chart_exclusions, coords, assumptions=assumptions, budget=budget)
    sn = certify_null_set(
        chart.source_exclusions, variables, assumptions=assumptions, budget=budget
    )
    if chart.kind.startswith("polar"):
        if len(maps) != 1:
            raise ValueError("Canonical polar charts require a scalar level map")
        return _polar_certificate(chart, g, variables, assumptions, budget, cn, sn)
    if chart.kind not in ("rational", "algebraic"):
        raise ValueError("Unknown chart kind")
    forward = dict(zip(variables, chart.forward, strict=True))
    inverse = dict(zip(coords, chart.inverse, strict=True))
    image = sp.And(
        full.subs(inverse, simultaneous=True), chart.source_domain, sp.Not(chart.source_exclusions)
    )
    source_premise = sp.And(*assumptions, image)
    for e in chart.forward:
        if e.has(*variables):
            raise ValueError("Forward formulas contain source coordinates")
        _formula(e, coords, regular, assumptions, chart.kind, budget)
    for e in chart.inverse:
        if e.has(*coords):
            raise ValueError("Inverse formulas contain chart coordinates")
        _formula(e, variables, image, assumptions, chart.kind, budget)
    if not _prove(
        premise, chart.source_domain.subs(forward, simultaneous=True), coords, budget
    ) or not _prove(
        premise, sp.Not(chart.source_exclusions.subs(forward, simultaneous=True)), coords, budget
    ):
        raise GeometryProofError("Forward chart leaves its regular source domain")
    if not _prove(
        source_premise,
        sp.Not(chart.chart_exclusions.subs(inverse, simultaneous=True)),
        variables,
        budget,
    ):
        raise GeometryProofError("Inverse chart leaves its regular parameter domain")
    if (
        any(
            not _identity(e.subs(forward, simultaneous=True) - v, premise)
            for e, v in zip(chart.inverse, coords, strict=True)
        )
        or any(
            not _identity(e.subs(inverse, simultaneous=True) - v, source_premise)
            for e, v in zip(chart.forward, variables, strict=True)
        )
        or any(
            not _identity(h.subs(forward, simultaneous=True) - t, premise)
            for h, t in zip(maps, coords, strict=False)
        )
    ):
        raise GeometryProofError("Chart inverse or level-map identity failed")
    matrix = sp.Matrix(chart.forward).jacobian(coords)
    det = sp.factor(_determinant(matrix, budget))
    if not _prove(premise, sp.Ne(det, 0), coords, budget):
        raise GeometryProofError("Chart full-rank prerequisite is unproved")
    fiber = matrix[:, len(maps) :]
    metric = sp.factor(_determinant(fiber.T * fiber, budget))
    gradient = sp.Matrix(maps).jacobian(variables)
    jac = sp.factor(_determinant(gradient * gradient.T, budget).subs(forward, simultaneous=True))
    if not _identity(metric - jac * det**2, premise):
        raise GeometryProofError("Coarea metric identity failed")
    interior = sp.And(
        _domain(chart.ranges, True).subs(inverse, simultaneous=True),
        chart.source_domain,
        sp.Not(chart.source_exclusions),
    )
    return ChartCertificate(
        chart, variables, g, assumptions, regular, image, interior, det, jac, metric, cn, sn
    )


def polar_chart(variables, coordinates, level_bounds, *, center=(0, 0), squared=False):
    """Canonical planar polar chart; the center and negative horizontal seam are excluded."""
    x, y = variables
    t, u = coordinates
    cx, cy = map(sp.sympify, center)
    lo, hi = map(sp.sympify, level_bounds)
    q = (x - cx) ** 2 + (y - cy) ** 2
    r = sp.sqrt(t) if squared else t
    return Chart(
        tuple(coordinates),
        (cx + r * sp.cos(u), cy + r * sp.sin(u)),
        (q if squared else sp.sqrt(q), sp.atan2(y - cy, x - cx)),
        ((u, -sp.pi, sp.pi), (t, lo, hi)),
        sp.Or(sp.Eq(t, 0), sp.Eq(u, -sp.pi), sp.Eq(u, sp.pi)),
        sp.And(sp.Eq(y, cy), x <= cx),
        sp.true,
        "polar_squared" if squared else "polar_radius",
    )


def _polar_certificate(chart, g, variables, assumptions, budget, cn, sn):
    if len(variables) != 2:
        raise GeometryProofError("Canonical polar charts are planar")
    t, u = chart.coordinates
    center = tuple(sp.simplify(e.subs({t: 0, u: 0})) for e in chart.forward)
    try:
        real_center = all(
            _local_ask(
                sp.Q.real(c) & sp.Q.finite(c), _predicates(sp.And(*assumptions)), context=set()
            )
            is True
            for c in center
        )
    except (ValueError, AssertionError, UnhandledInput):
        real_center = False
    if not real_center:
        raise GeometryProofError("Polar center must have finite real coordinates")
    canonical = polar_chart(
        variables,
        chart.coordinates,
        chart.ranges[-1][1:],
        center=center,
        squared=chart.kind == "polar_squared",
    )
    if chart != canonical:
        raise GeometryProofError("Polar chart differs from the checked canonical family")
    lo, hi = chart.ranges[-1][1:]
    premise = sp.And(*assumptions)
    if not _prove(premise, lo >= 0, variables, budget):
        raise GeometryProofError("Polar levels must be nonnegative")
    q = sum((v - c) ** 2 for v, c in zip(variables, center, strict=True))
    expected = q if chart.kind == "polar_squared" else sp.sqrt(q)
    if not _identity(g - expected, sp.And(*assumptions)):
        raise GeometryProofError("Polar level-map identity failed")
    # Injectivity for -pi < angle < pi follows from the canonical polar theorem.
    # Polynomial level inequalities avoid treating atan2 as semialgebraic.
    image = sp.And(
        q >= lo if chart.kind == "polar_squared" else q >= lo**2,
        q <= hi if chart.kind == "polar_squared" else q <= hi**2,
        sp.Not(chart.source_exclusions),
    )
    interior = sp.And(
        q > lo if chart.kind == "polar_squared" else q > lo**2,
        q < hi if chart.kind == "polar_squared" else q < hi**2,
        sp.Not(chart.source_exclusions),
    )
    det = sp.Rational(1, 2) if chart.kind == "polar_squared" else t
    jac = 4 * t if chart.kind == "polar_squared" else sp.S.One
    metric = t if chart.kind == "polar_squared" else t * t
    return ChartCertificate(
        chart,
        variables,
        g,
        assumptions,
        sp.And(_domain(chart.ranges), sp.Not(chart.chart_exclusions)),
        image,
        interior,
        det,
        jac,
        metric,
        cn,
        sn,
        "canonical polar inverse and seam theorem",
    )


def certify_cover(region, certificates, *, assumptions=None, budget=DEFAULT_GEOMETRY_BUDGET):
    return _cover_cached(sp.sympify(region), tuple(certificates), _assumptions(assumptions), budget)


@lru_cache(maxsize=128)
def _cover_cached(region, certificates, assumptions, budget):
    if not certificates or len(certificates) > budget.max_charts:
        raise GeometryProofError("Chart count is empty or exceeds budget")
    first = certificates[0]
    variables = first.variables
    _bounded((region, sp.And(*assumptions)), variables, budget)
    if any(
        c.variables != variables or c.level_map != first.level_map or c.assumptions != assumptions
        for c in certificates
    ):
        raise ValueError("Cover certificates must use identical coordinates, map and assumptions")
    # Do not trust manually constructed certificate records: replay cached certification.
    if any(
        certify_chart(c.chart, c.level_map, c.variables, assumptions=assumptions, budget=budget)
        != c
        for c in certificates
    ):
        raise GeometryProofError("Chart certificate replay failed")
    null = sp.Or(*(c.source_null.region for c in certificates))
    nc = certify_null_set(null, variables, assumptions=assumptions, budget=budget)
    union = sp.Or(*(c.image for c in certificates))
    premise = sp.And(*assumptions)
    if not _prove(sp.And(premise, region, sp.Not(null)), union, variables, budget) or not _prove(
        sp.And(premise, union), region, variables, budget
    ):
        raise GeometryProofError("Source coverage equivalence modulo null set is unproved")
    for i, c in enumerate(certificates):
        for other in certificates[i + 1 :]:
            if not _prove(
                sp.And(premise, c.interior_image, other.interior_image), sp.false, variables, budget
            ):
                raise GeometryProofError("Chart interiors overlap or disjointness is unproved")
    return CoverCertificate(region, variables, assumptions, certificates, nc)


def clear_chart_cache():
    for f in (_null_cached, _chart_cached, _cover_cached, _linear_cached):
        f.cache_clear()


def chart_cache_info():
    return tuple(
        f.cache_info() for f in (_null_cached, _chart_cached, _cover_cached, _linear_cached)
    )


@dataclass(frozen=True)
class LinearMetricCertificate:
    """Full-column-rank linear differential, usable for ambient or intrinsic measure."""

    matrix: sp.ImmutableMatrix
    gram_determinant: sp.Expr
    determinant: sp.Expr | None
    assumptions: tuple
    basis: str = "positive Gram determinant proves full column rank"


def certify_linear_metric(matrix, *, assumptions=None, budget=DEFAULT_GEOMETRY_BUDGET):
    return _linear_cached(sp.ImmutableMatrix(matrix), _assumptions(assumptions), budget)


@lru_cache(maxsize=256)
def _linear_cached(matrix, assumptions, budget):
    if matrix.cols < 1 or matrix.rows < matrix.cols:
        raise GeometryProofError("Linear map cannot have full column rank")
    if matrix.cols > budget.max_variables:
        raise GeometryProofError("Linear differential dimension budget exceeded")
    premise = sp.And(*assumptions)
    variables = tuple(
        sorted(set().union(*(v.free_symbols for v in matrix)), key=sp.default_sort_key)
    )
    _bounded(tuple(matrix), variables, budget)
    for c in matrix:
        try:
            real = _local_ask(sp.Q.real(c) & sp.Q.finite(c), _predicates(premise), context=set())
        except (ValueError, AssertionError, UnhandledInput):
            real = None
        if real is not True:
            raise GeometryProofError("Linear differential must have finite real entries")
    gram = sp.factor(_determinant(matrix.T * matrix, budget))
    if not _prove(premise, gram > 0, variables, budget):
        raise GeometryProofError("Linear differential full rank is unproved")
    determinant = sp.factor(_determinant(matrix, budget)) if matrix.rows == matrix.cols else None
    return LinearMetricCertificate(matrix, gram, determinant, assumptions)


def chart_parameter_events(
    chart, level_map, variables, parameters, *, budget=DEFAULT_GEOMETRY_BUDGET
):
    """Bounded parameter-only rank, pole and boundary-coincidence candidates.

    Returned expressions are event candidates, not a partition or a claim that
    their equality fibers are null. Caller cells must still be proved separately.
    """
    if not isinstance(chart, Chart):
        raise ValueError("Expected a Chart candidate")
    variables, parameters = tuple(variables), tuple(parameters)
    if set(parameters) & (set(variables) | set(chart.coordinates)) or any(
        not isinstance(p, sp.Symbol) or p.is_real is not True for p in parameters
    ):
        raise ValueError("Use distinct real parameter coordinates")
    _bounded(
        (*chart.forward, *chart.inverse, *[hi - lo for _, lo, hi in chart.ranges]),
        chart.coordinates,
        budget,
    )
    maps = (
        tuple(level_map)
        if isinstance(level_map, (tuple, list, sp.Tuple))
        else (sp.sympify(level_map),)
    )
    _bounded(maps, variables, budget)
    matrix = sp.Matrix(chart.forward).jacobian(chart.coordinates)
    expressions = [_determinant(matrix, budget), *[hi - lo for _, lo, hi in chart.ranges]]
    gradient = sp.Matrix(maps).jacobian(variables)
    expressions.append(_determinant(gradient * gradient.T, budget))
    return parameter_geometry_events(
        (*expressions, *chart.forward, *chart.inverse), parameters, budget=budget
    )


def parameter_geometry_events(expressions, parameters, *, budget=DEFAULT_GEOMETRY_BUDGET):
    """Factor bounded parameter-only endpoint, coefficient and pole candidates."""
    parameters = tuple(parameters)
    expressions = tuple(map(sp.sympify, expressions))
    if any(not isinstance(p, sp.Symbol) or p.is_real is not True for p in parameters) or len(
        set(parameters)
    ) != len(parameters):
        raise ValueError("Use distinct real parameters")
    _bounded(expressions, parameters, budget)
    events = set()
    for expression in expressions:
        for term in sp.fraction(sp.cancel(expression)):
            if not term.free_symbols or not term.free_symbols <= set(parameters):
                continue
            if term.is_Pow and term.exp.is_Rational:
                term = term.base
            try:
                factors = sp.factor_list(term, *parameters)[1]
            except sp.PolynomialError:
                continue
            events.update(sp.factor(factor) for factor, _ in factors if factor.has(*parameters))
    return tuple(sorted(events, key=sp.default_sort_key))
