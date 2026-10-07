from __future__ import annotations

from dataclasses import dataclass

import sympy as sp
from sympy import S
from sympy.polys.polyerrors import PolynomialError

from ..engine import simplify_condition
from ..exactness import certify_identity_equal, certify_identity_zero
from ..results import Homogeneity, Periodicity, Symmetry
from ..transcendental.periodic import common_period_values, positive_period
from .common import normalize_variables


@dataclass(frozen=True)
class PeriodicityResult:
    classification: Periodicity
    period: sp.Expr | None = None
    domain: object = S.Reals
    evidence: tuple = ()
    fundamental: bool | None = None
    conditions: sp.Expr = S.true


@dataclass(frozen=True)
class SymmetryResult:
    classification: Symmetry
    center: sp.Expr = S.Zero
    domain: object = S.Reals
    evidence: tuple = ()
    centers: object = ()
    center_period: sp.Expr | None = None
    central_value: sp.Expr | None = None
    center_families: tuple = ()


@dataclass(frozen=True)
class HomogeneityCase:
    degree: sp.Expr | None
    condition: sp.Expr
    weights: tuple


@dataclass(frozen=True)
class HomogeneityResult:
    classification: Homogeneity
    degree: sp.Expr | None = None
    variables: tuple = ()
    evidence: tuple = ()
    weights: tuple = ()
    conditions: sp.Expr = S.true
    cases: tuple[HomogeneityCase, ...] = ()
    assumptions: sp.Expr = S.true


def _zero(expr):
    return certify_identity_zero(expr)


def _period_from_frequency(arg, var):
    """Return exact positive period for sin/cos-like affine arguments."""
    try:
        a = sp.simplify(sp.diff(arg, var))
        if a.has(var) or a == 0:
            return None
        if certify_identity_zero(sp.diff(arg, var, 2)) is False:
            return None
        p = sp.simplify(2 * sp.pi / sp.Abs(a))
        return p if p.is_positive is True else None
    except (TypeError, ValueError, NotImplementedError):
        return None


def _period_of(expr, var):
    """Conservative exact period propagation; returned period need not be minimal."""
    expr = sp.sympify(expr)
    if not expr.has(var):
        return None
    p = positive_period(expr, var)
    if p is not None:
        return p
    if expr.func in (sp.sin, sp.cos) and len(expr.args) == 1:
        return _period_from_frequency(expr.args[0], var)
    if expr.func in (sp.tan, sp.cot) and len(expr.args) == 1:
        q = _period_from_frequency(expr.args[0], var)
        return sp.simplify(q / 2) if q is not None else None
    # Any pointwise outer function of a vector of periodic arguments inherits
    # every common period of those arguments. Constants do not constrain it.
    if expr.is_Function or (expr.func not in (sp.Add, sp.Mul, sp.Pow) and expr.args):
        varying = [a for a in expr.args if a.has(var)]
        if varying:
            ps = [_period_of(a, var) for a in varying]
            if all(q is not None for q in ps):
                return common_period_values(*ps)
    if expr.is_Add or expr.is_Mul:
        varying = [a for a in expr.args if a.has(var)]
        ps = [_period_of(a, var) for a in varying]
        if varying and all(q is not None for q in ps):
            return common_period_values(*ps)
    if expr.is_Pow and not expr.exp.has(var):
        return _period_of(expr.base, var)
    return None


def _visible_incommensurate_periods(expr, var):
    if not (expr.is_Add or expr.is_Mul):
        return False
    ps = []
    for a in expr.args:
        if not a.has(var):
            continue
        p = _period_of(a, var)
        if p is None:
            return False
        ps.append(p)
    return len(ps) > 1 and common_period_values(*ps) is None


def periodicity(expr, var, *, domain=True, assumptions=True, return_result=False):
    """Classify exact real periodicity and return a certified period in detail mode.

    A period propagated through composition/sums/products is certified but is
    not claimed fundamental unless the underlying exact periodicity engine
    establishes the same period directly.
    """
    expr = sp.sympify(expr)
    var = sp.sympify(var)
    if not expr.has(var):
        r = PeriodicityResult(
            Periodicity.CONSTANT, None, domain, ("independent-of-variable",), False
        )
        return r if return_result else r.classification
    direct = positive_period(expr, var)
    p = direct if direct is not None else _period_of(expr, var)
    # Preserve a symbolic affine-frequency theorem
    # Unknown merely because the frequency may vanish on a parameter stratum.
    if (
        p is None
        and expr.func in (sp.sin, sp.cos, sp.tan, sp.cot)
        and len(expr.args) == 1
    ):
        arg = expr.args[0]
        try:
            a = sp.simplify(sp.diff(arg, var))
            b = sp.simplify(arg - a * var)
            affine = (
                not a.has(var)
                and certify_identity_zero(sp.diff(arg, var, 2)) is True
                and not b.has(var)
            )
        except (TypeError, ValueError, NotImplementedError):
            affine = False
            a = None
        if affine and a is not None and a != 0:
            real_guard = simplify_condition(sp.Q.real(a), assumptions=assumptions)
            if real_guard is S.false:
                real_guard = S.false
            zero_guard = simplify_condition(
                sp.And(real_guard, sp.Eq(a, 0)), assumptions=assumptions
            )
            nonzero_guard = simplify_condition(
                sp.And(real_guard, sp.Ne(a, 0)), assumptions=assumptions
            )
            if zero_guard is S.true:
                r = PeriodicityResult(
                    Periodicity.CONSTANT,
                    None,
                    domain,
                    ("zero-affine-frequency",),
                    False,
                )
                return r if return_result else r.classification
            base = sp.pi if expr.func in (sp.tan, sp.cot) else 2 * sp.pi
            cp = sp.simplify(base / sp.Abs(a))
            if nonzero_guard is S.false:
                # No real nonzero-frequency stratum is compatible with the
                # assumptions, so this theorem does not apply.
                nonzero_guard = None
            cond = nonzero_guard
            if cond is not None:
                evidence = (
                    ("affine-frequency",)
                    if cond is S.true
                    else ("conditional-affine-frequency",)
                )
                r = PeriodicityResult(
                    Periodicity.PERIODIC, cp, domain, evidence, None, cond
                )
                return r if return_result else r.classification
    if p is not None:
        # Replay the shift identity before accepting structural propagation.
        z = _zero(expr.subs(var, var + p) - expr)
        if z is True:
            primitive_minimal = (
                expr.func in (sp.sin, sp.cos, sp.tan, sp.cot)
                and len(expr.args) == 1
                and _period_from_frequency(expr.args[0], var) is not None
            )
            fundamental = (
                True
                if primitive_minimal
                and direct is not None
                and certify_identity_equal(direct, p) is True
                else None
            )
            ev = (
                ("exact-periodicity",)
                if direct is not None
                else ("commensurate-period-propagation",)
            )
            r = PeriodicityResult(
                Periodicity.PERIODIC, sp.simplify(p), domain, ev, fundamental
            )
            return r if return_result else r.classification
    try:
        if expr.is_polynomial(var) and sp.Poly(expr, var).degree() > 0:
            result = PeriodicityResult(
                Periodicity.APERIODIC,
                None,
                domain,
                ("nonconstant-polynomial",),
                None,
            )
            return result if return_result else result.classification
    except PolynomialError:
        pass
    if _visible_incommensurate_periods(expr, var):
        r = PeriodicityResult(
            Periodicity.APERIODIC,
            None,
            domain,
            ("incommensurate-component-periods",),
            None,
        )
        return r if return_result else r.classification
    r = PeriodicityResult(Periodicity.UNKNOWN, None, domain, ())
    return r if return_result else r.classification


def _classify_about(expr, var, c, *, allow_central=False):
    reflected = sp.simplify(expr.subs(var, 2 * c - var))
    even = _zero(reflected - expr)
    odd = _zero(reflected + expr)
    k = sp.simplify(expr.subs(var, c))
    central = _zero(reflected + expr - 2 * k) if k is not None else None
    if even is True and odd is True:
        cls = Symmetry.BOTH
    elif even is True:
        cls = Symmetry.EVEN
    elif odd is True:
        cls = Symmetry.ODD
    elif allow_central and central is True:
        cls = Symmetry.CENTRAL
    elif even is False and odd is False and (central is False or not allow_central):
        cls = Symmetry.NEITHER
    else:
        cls = Symmetry.UNKNOWN
    return cls, reflected, k


def _polynomial_reflection_center(expr, var):
    try:
        p = sp.Poly(expr, var)
        n = p.degree()
        if n <= 0:
            return S.Zero
        if n == 1:
            return ("all-centers-central",)
        c = sp.Dummy("_center", real=True)
        t = sp.Dummy("_offset", real=True)
        shifted = sp.Poly(sp.expand(expr.subs(var, c + t)), t)
        # Even/reflection symmetry requires all odd powers to vanish. Central
        # symmetry requires all positive even powers to vanish. Solve either
        # exact coefficient system and replay on the original expression.
        systems = []
        odd = [coef for (m,), coef in shifted.terms() if m % 2]
        evenpos = [coef for (m,), coef in shifted.terms() if m > 0 and m % 2 == 0]
        if odd:
            systems.append(odd)
        if evenpos:
            systems.append(evenpos)
        vals = []
        for equations in systems:
            sol = sp.solve(equations, c, dict=False)
            if isinstance(sol, dict):
                sol = [sol.get(c)]
            if not isinstance(sol, (tuple, list, set)):
                sol = [sol]
            for q in sol:
                if isinstance(q, (tuple, list)):
                    q = q[0] if len(q) == 1 else None
                if q is not None and not sp.sympify(q).has(c):
                    vals.append(sp.simplify(q))
        vals = tuple(dict.fromkeys(vals))
        return vals[0] if len(vals) == 1 else vals
    except (PolynomialError, TypeError, ValueError, NotImplementedError):
        return None


def _trig_center_family(expr, var):
    # For a single affine sin/cos, return a representative reflection center
    # and period of equivalent centers. This remains exact and compact.
    if expr.func not in (sp.sin, sp.cos) or len(expr.args) != 1:
        return None
    arg = expr.args[0]
    try:
        a = sp.simplify(sp.diff(arg, var))
        b = sp.simplify(arg - a * var)
        if a.has(var) or a == 0:
            return None
        step = sp.simplify(sp.pi / sp.Abs(a))
        c0 = sp.simplify(-b / a)
        half = sp.simplify(sp.pi / (2 * a))
        if expr.func is sp.cos:
            cls = Symmetry.EVEN
            families = (
                (Symmetry.EVEN, c0, step),
                (Symmetry.ODD, sp.simplify(c0 + half), step),
            )
        else:
            cls = Symmetry.ODD
            families = (
                (Symmetry.ODD, c0, step),
                (Symmetry.EVEN, sp.simplify(c0 + half), step),
            )
        return c0, step, cls, families
    except (TypeError, ValueError, NotImplementedError):
        return None


def symmetry(expr, var, *, center=0, domain=True, return_result=False):
    """Classify reflection symmetry about ``center``.

    Passing ``center=None`` requests exact center discovery. Detail mode also
    reports periodic center families when a compact exact family is known.
    """
    expr = sp.sympify(expr)
    var = sp.sympify(var)
    if center is None:
        fam = _trig_center_family(expr, var)
        if fam is not None:
            c, step, cls, families = fam
            got, _, k = _classify_about(expr, var, c, allow_central=True)
            if got is cls:
                r = SymmetryResult(
                    cls,
                    c,
                    domain,
                    ("affine-trigonometric-center",),
                    (),
                    step,
                    k,
                    families,
                )
                return r if return_result else cls
        pc = _polynomial_reflection_center(expr, var)
        if isinstance(pc, tuple) and pc == ("all-centers-central",):
            c = S.Zero
            cls, ref, k = _classify_about(expr, var, c, allow_central=True)
            r = SymmetryResult(
                cls, c, domain, ("affine-all-centers-central",), S.Reals, None, k, ()
            )
            return r if return_result else cls
        candidates = pc if isinstance(pc, tuple) else (() if pc is None else (pc,))
        # A pointwise outer function preserves every even reflection center of
        # its varying argument. Verify on the full expression before accepting.
        if not candidates and len([a for a in expr.args if a.has(var)]) == 1:
            inner = next(a for a in expr.args if a.has(var))
            ir = symmetry(inner, var, center=None, domain=domain, return_result=True)
            if ir.classification in (Symmetry.EVEN, Symmetry.BOTH):
                candidates = (ir.center,)
        for c in candidates:
            cls, ref, k = _classify_about(expr, var, c, allow_central=True)
            if cls in (Symmetry.EVEN, Symmetry.ODD, Symmetry.BOTH, Symmetry.CENTRAL):
                r = SymmetryResult(
                    cls,
                    c,
                    domain,
                    (sp.Eq(ref, expr), sp.Eq(ref, -expr)),
                    candidates,
                    None,
                    k,
                    (),
                )
                return r if return_result else cls
        r = SymmetryResult(Symmetry.UNKNOWN, S.Zero, domain, (), candidates, None)
        return r if return_result else r.classification
    c = sp.sympify(center)
    cls, reflected, k = _classify_about(expr, var, c)
    r = SymmetryResult(
        cls,
        c,
        domain,
        (sp.Eq(reflected, expr), sp.Eq(reflected, -expr)),
        (c,),
        None,
        k,
        (),
    )
    return r if return_result else cls


def _normalize_weights(vars_, weights):
    if weights is None:
        return tuple(S.One for _ in vars_)
    if isinstance(weights, dict):
        return tuple(sp.sympify(weights.get(v, S.One)) for v in vars_)
    ww = tuple(sp.sympify(w) for w in weights)
    if len(ww) != len(vars_):
        raise ValueError("weights must match variables")
    return ww


def _weighted_poly_cases(expr, vars_, weights):
    """Return exact degree and parameter strata for a polynomial."""
    try:
        polynomial = sp.Poly(expr, *vars_)
    except PolynomialError:
        return None
    p = polynomial
    if p.is_zero:
        return (HomogeneityCase(None, S.true, weights),)
    grouped = {}
    for monom, coeff in p.terms():
        d = sp.simplify(
            sum(sp.Integer(m) * w for m, w in zip(monom, weights, strict=False))
        )
        grouped.setdefault(d, []).append(coeff)
    if len(grouped) == 1:
        d = next(iter(grouped))
        return (HomogeneityCase(d, S.true, weights),)
    cases = []
    all_coeffs = [sp.simplify(c) for coeffs in grouped.values() for c in coeffs]
    zero_cond = sp.And(*(sp.Eq(c, 0) for c in all_coeffs)) if all_coeffs else S.true
    cases.append(HomogeneityCase(None, zero_cond, weights))
    for d in grouped:
        vanish = []
        for od, coeffs in grouped.items():
            if od == d:
                continue
            vanish.extend(sp.Eq(sp.simplify(c), 0) for c in coeffs)
        # A weighted-degree block survives when *at least one* coefficient
        # in that block is nonzero.  Requiring every coefficient to be nonzero
        # incorrectly drops homogeneous parameter strata such as a*x + b*y
        # with a=0, b!=0.
        keep_terms = [sp.Ne(sp.simplify(c), 0) for c in grouped[d]]
        keep = sp.Or(*keep_terms) if keep_terms else S.false
        cond = sp.And(*([*vanish, keep])) if (vanish or keep_terms) else S.true
        cases.append(HomogeneityCase(d, cond, weights))
    return tuple(cases)


def _infer_homogeneous_degree(expr, vars_, weights):
    cases = _weighted_poly_cases(expr, vars_, weights)
    if cases is not None:
        if len(cases) == 1 and cases[0].condition is S.true:
            return cases[0].degree, cases
        return False, cases
    if expr == 0:
        return None, (HomogeneityCase(None, S.true, weights),)
    euler = sp.simplify(
        sum(w * v * sp.diff(expr, v) for v, w in zip(vars_, weights, strict=False))
    )
    ratio = sp.cancel(euler / expr)
    if not any(ratio.has(v) for v in vars_) and ratio.is_real is not False:
        if _zero(euler - ratio * expr) is True:
            return sp.simplify(ratio), (
                HomogeneityCase(sp.simplify(ratio), S.true, weights),
            )
    return None, ()


def homogeneity(
    expr, vars, *, degree=None, weights=None, assumptions=True, return_result=False
):
    """Classify positive weighted homogeneity in one or more variables.

    ``weights`` defaults to ordinary simultaneous scaling. Detail mode retains
    exact parameter strata when a polynomial becomes homogeneous only after
    some coefficients vanish.
    """
    expr = sp.sympify(expr)
    vv = normalize_variables(vars)
    ww = _normalize_weights(vv, weights)
    inferred, cases = _infer_homogeneous_degree(expr, vv, ww)
    assumptions = sp.sympify(assumptions)
    if expr == 0:
        if degree is not None:
            d = sp.sympify(degree)
            cls = Homogeneity.HOMOGENEOUS
            rd = d
        else:
            cls = Homogeneity.HOMOGENEOUS
            rd = None
        r = HomogeneityResult(
            cls, rd, vv, ("zero-function-all-degrees",), ww, S.true, cases, assumptions
        )
        return r if return_result else cls
    if degree is not None:
        d = sp.sympify(degree)
        if inferred not in (None, False):
            z = _zero(inferred - d)
            cls = (
                Homogeneity.HOMOGENEOUS
                if z is True
                else Homogeneity.NOT_HOMOGENEOUS
                if z is False
                else Homogeneity.UNKNOWN
            )
        else:
            t = sp.Dummy("_scale", positive=True)
            scaled = expr.xreplace({v: t**w * v for v, w in zip(vv, ww, strict=False)})
            z = _zero(scaled - t**d * expr)
            cls = (
                Homogeneity.HOMOGENEOUS
                if z is True
                else Homogeneity.NOT_HOMOGENEOUS
                if z is False
                else Homogeneity.UNKNOWN
            )
        rd = d if cls is Homogeneity.HOMOGENEOUS else None
    else:
        if inferred is False:
            applicable = []
            for case in cases:
                try:
                    implication = sp.simplify(sp.Implies(assumptions, case.condition))
                except (TypeError, ValueError, NotImplementedError):
                    implication = None
                if implication is S.true:
                    applicable.append(case)
            if len(applicable) == 1:
                cls = Homogeneity.HOMOGENEOUS
                rd = applicable[0].degree
            else:
                homogeneous_region = sp.Or(*(c.condition for c in cases))
                try:
                    excludes_all = sp.simplify(
                        sp.Implies(assumptions, sp.Not(homogeneous_region))
                    )
                except (TypeError, ValueError, NotImplementedError):
                    excludes_all = None
                cls = (
                    Homogeneity.NOT_HOMOGENEOUS
                    if excludes_all is S.true
                    else Homogeneity.UNKNOWN
                )
                rd = None
        elif inferred is None:
            cls = Homogeneity.UNKNOWN
            rd = None
        else:
            cls = Homogeneity.HOMOGENEOUS
            rd = sp.sympify(inferred)
    cond = S.true
    if cases and cls is Homogeneity.NOT_HOMOGENEOUS:
        cond = sp.Not(sp.Or(*(c.condition for c in cases)))
    r = HomogeneityResult(
        cls,
        rd,
        vv,
        ("weighted-Euler-homogeneity",) if cls is Homogeneity.HOMOGENEOUS else (),
        ww,
        cond,
        cases,
        assumptions,
    )
    return r if return_result else cls
