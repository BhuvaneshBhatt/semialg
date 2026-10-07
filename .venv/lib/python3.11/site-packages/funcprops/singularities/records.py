"""Internal singularity analysis implementation."""

from __future__ import annotations

from collections.abc import Iterable
from functools import lru_cache

import sympy as sp
from sympy import S
from sympy.core.relational import Relational
from sympy.logic.boolalg import And, Or

from ..engine import evaluate_property, simplify_condition
from ..exactness import certify_identity_equal, certify_zero
from ..exceptions import SYMBOLIC_QUERY_ERRORS
from ..results import (
    Singularity,
    SingularityKind,
)


def _variable_tuple(var):
    if isinstance(var, (tuple, list, sp.Tuple)):
        variables = tuple(sp.sympify(v) for v in var)
    else:
        variables = (sp.sympify(var),)
    if not variables:
        raise ValueError("at least one variable is required")
    return variables


def _is_empty(value) -> bool:
    return value in ({}, (), [], S.EmptySet, sp.Tuple())


def _replace_helper_heads(expr):
    """Normalize a few helper heads used only by the catalogue property dump."""
    if not isinstance(expr, sp.Basic):
        return expr
    name = getattr(expr.func, "__name__", "")
    if name == "RemoveQuadratics" and len(expr.args) == 1:
        a = expr.args[0]
        if getattr(a, "func", None) is sp.Abs:
            return _replace_helper_heads(a.args[0])
        if isinstance(a, sp.Pow):
            return _replace_helper_heads(a.args[0])
        return _replace_helper_heads(a)
    if name == "If" and len(expr.args) == 3:
        cond = _replace_helper_heads(expr.args[0])
        if cond in (S.true, True):
            return _replace_helper_heads(expr.args[1])
        if cond in (S.false, False):
            return _replace_helper_heads(expr.args[2])
    args = tuple(_replace_helper_heads(a) for a in expr.args)
    try:
        return expr.func(*args)
    except SYMBOLIC_QUERY_ERRORS:
        return expr


def _combine_guard(a, b, assumptions=S.true):
    return simplify_condition(And(a, b), assumptions=assumptions)


def _split_locus_guard(condition, var, *, assumptions=S.true):
    """Split a condition into analysis-locus and parameter-only parts.

    ``var`` may be one symbol or a collection of independent analysis variables.
    A Boolean atom belongs to the locus whenever it depends on at least one of
    those variables; only genuinely parameter-only atoms are moved to the guard.
    """
    c = sp.sympify(condition)
    variables = _variable_tuple(var) if var is not None else ()

    def depends_on_analysis_variable(node):
        return any(node.has(variable) for variable in variables)

    if c in (S.true, S.false):
        return c, S.true
    if c.func is And:
        locus = []
        guard = []
        for atom in c.args:
            (locus if depends_on_analysis_variable(atom) else guard).append(atom)
        lc = And(*locus) if locus else S.true
        gd = And(*guard) if guard else S.true
        return simplify_condition(lc, assumptions=assumptions), simplify_condition(
            gd, assumptions=assumptions
        )
    if depends_on_analysis_variable(c):
        return c, S.true
    return S.true, simplify_condition(c, assumptions=assumptions)


def _registry_loci(
    expr,
    prop: str,
    kind: SingularityKind,
    *,
    mode="complex",
    var=None,
    assumptions=S.true,
):
    result = evaluate_property(expr, prop, var=var, mode=mode, assumptions=assumptions)
    out = []
    for case in result.cases:
        value = _replace_helper_heads(case.value)
        if _is_empty(value):
            continue
        entries = (
            tuple(value) if isinstance(value, (sp.Tuple, tuple, list)) else (value,)
        )
        for entry in entries:
            entry = _replace_helper_heads(entry)
            if isinstance(entry, (sp.Tuple, tuple, list)) and len(entry) >= 2:
                locus, local_condition = entry[0], entry[1]
                local_locus, local_guard = _split_locus_guard(
                    local_condition, var, assumptions=assumptions
                )
                locus = _combine_guard(locus, local_locus, assumptions=assumptions)
                guard = _combine_guard(
                    case.condition, local_guard, assumptions=assumptions
                )
            else:
                locus, guard = entry, case.condition
            if guard is S.false or locus is S.false:
                continue
            out.append(
                Singularity(
                    kind, sp.sympify(locus), sp.sympify(guard), source="registry"
                )
            )
    return tuple(out)


def _unique(items: Iterable[Singularity]):
    out = []
    seen = set()
    for item in items:
        key = (
            item.kind,
            sp.srepr(item.condition),
            sp.srepr(item.guard),
            str(item.multiplicity),
        )
        if key in seen:
            continue
        seen.add(key)
        out.append(item)
    return tuple(out)


def order_at(expr, var, point):
    """Return the local integer/Puiseux order when SymPy can certify it.

    Positive means a zero, negative a pole, and zero a finite nonzero value.
    ``None`` means the order could not be established symbolically.
    """
    expr = sp.sympify(expr)
    var = sp.sympify(var)
    point = sp.sympify(point)
    t = sp.Dummy("funcprops_local")
    try:
        shifted = expr.subs(var, point + t)
        lead = shifted.as_leading_term(t)
        pd = lead.as_powers_dict()
        order = sp.simplify(pd.get(t, S.Zero))
        if order.is_number:
            # A formal leading power is useful only when its coefficient is
            # certified nonzero.  For example f(t)/t has apparent order -1,
            # but f(0) may vanish and cancel the denominator.
            coefficient = sp.simplify(lead / t**order)
            coefficient_at_point = sp.simplify(coefficient.subs(t, S.Zero))
            if certify_zero(coefficient_at_point) is False:
                return order
    except SYMBOLIC_QUERY_ERRORS:
        pass
    # A short series fallback catches cancellation such as sin(x)/x**3.
    try:
        ser = sp.series(expr, var, point, 8).removeO().expand()
        shifted = ser.subs(var, point + t).expand()
        terms = sp.Add.make_args(shifted)
        certified_orders = []
        for term in terms:
            power = term.as_powers_dict().get(t, S.Zero)
            if not power.is_number:
                continue
            coefficient = sp.simplify(term / t**power)
            coefficient_at_point = sp.simplify(coefficient.subs(t, S.Zero))
            if certify_zero(coefficient_at_point) is False:
                certified_orders.append(power)
        if certified_orders:
            return min(certified_orders, key=lambda q: float(q))
    except SYMBOLIC_QUERY_ERRORS:
        pass
    return None


@lru_cache(maxsize=256)
def _denominator_singularity_records(expr, var):
    """Classify explicit univariate denominator zeros in one structural pass.

    Returns ``(poles, removable)``.  Keeping the two classifications together
    avoids repeating ``together``/factorization/root solving/local-order work
    when a caller needs the complete singularity classification.
    """
    expr = sp.sympify(expr)
    var = sp.sympify(var)
    try:
        num, den = sp.fraction(sp.together(expr))
    except SYMBOLIC_QUERY_ERRORS:
        return (), ()
    if den == 1 or not den.has(var):
        return (), ()

    poles_out = []
    removable_out = []
    try:
        poly = sp.Poly(den, var)
        _, factors = sp.factor_list(poly.as_expr(), var)
    except SYMBOLIC_QUERY_ERRORS:
        return (
            Singularity(
                SingularityKind.POLE,
                sp.Eq(den, 0),
                S.true,
                None,
                source="denominator",
            ),
        ), ()

    for factor, mult in factors:
        try:
            roots = sp.solve(sp.Eq(factor, 0), var)
        except SYMBOLIC_QUERY_ERRORS:
            roots = []
        if not roots or any(root.has(var) for root in roots):
            poles_out.append(
                Singularity(
                    SingularityKind.POLE,
                    sp.Eq(factor, 0),
                    S.true,
                    mult,
                    source="denominator",
                )
            )
            continue

        for root in roots:
            order = order_at(expr, var, root)
            if order is not None and order.is_negative:
                multiplicity = -order
                poles_out.append(
                    Singularity(
                        SingularityKind.POLE,
                        sp.Eq(var, root),
                        S.true,
                        int(multiplicity) if multiplicity.is_Integer else multiplicity,
                        source="local-order",
                    )
                )
            elif order is not None and order.is_nonnegative:
                removable_out.append(
                    Singularity(
                        SingularityKind.REMOVABLE,
                        sp.Eq(var, root),
                        S.true,
                        source="local-order",
                    )
                )
            else:
                # A denominator zero alone does not certify a pole: an opaque or
                # parameter-dependent numerator may vanish there and cancel it.
                # Preserve the exact non-cancellation condition as a guard.
                try:
                    numerator_at_root = sp.simplify(num.subs(var, root))
                    nonzero = certify_zero(numerator_at_root)
                except SYMBOLIC_QUERY_ERRORS:
                    numerator_at_root = num.subs(var, root)
                    nonzero = None
                if nonzero is False:
                    guard = S.true
                elif nonzero is True:
                    # Cancellation is established, but its order is unknown.
                    # Do not invent either a pole or a removable classification.
                    continue
                else:
                    guard = sp.Ne(numerator_at_root, 0)
                poles_out.append(
                    Singularity(
                        SingularityKind.POLE,
                        sp.Eq(var, root),
                        guard,
                        mult,
                        source="denominator-noncancellation",
                    )
                )
    return _unique(poles_out), _unique(removable_out)


def _explicit_denominator_poles(expr, var):
    """Extract exact or conditional poles from an explicit denominator."""
    return _denominator_singularity_records(expr, var)[0]


def _recursive_argument_singularities(expr, var, kind, assumptions=S.true):
    out = []
    for arg in getattr(expr, "args", ()):
        if not getattr(arg, "has", lambda *_: False)(var):
            continue
        if kind is SingularityKind.POLE:
            out.extend(_function_pole_records(arg, var, assumptions=assumptions))
        elif kind is SingularityKind.BRANCH_POINT:
            out.extend(_function_branchpoint_records(arg, var, assumptions=assumptions))
        elif kind is SingularityKind.ESSENTIAL:
            out.extend(essential_singularities(arg, var, assumptions=assumptions))
    return out


def _validate_pole_candidates(expr, var, candidates):
    out = []
    for p in candidates:
        # Primitive registry/DLMF poles are already exact facts.  Local-order
        # probing is for algebraic denominator candidates and must not erase a
        # special-function pole merely because SymPy cannot expand locally.
        if p.source == "registry":
            out.append(p)
            continue
        cond = p.condition
        point = None
        if isinstance(cond, sp.Equality):
            if cond.lhs == var and not cond.rhs.has(var):
                point = cond.rhs
            elif cond.rhs == var and not cond.lhs.has(var):
                point = cond.lhs
        if point is not None:
            order = order_at(expr, var, point)
            if order is not None:
                if order.is_negative:
                    mult = -order
                    out.append(
                        Singularity(
                            SingularityKind.POLE,
                            cond,
                            p.guard,
                            int(mult) if mult.is_Integer else mult,
                            source="local-order",
                        )
                    )
                # nonnegative means cancellation/regularity, drop candidate
                continue
        out.append(p)
    return _unique(out)


def removable_singularities(expr, var):
    """Return denominator singularities certified as locally removable."""
    return _denominator_singularity_records(expr, var)[1]


def _function_pole_records(expr, var, *, assumptions=S.true, denominator_poles=None):
    expr = sp.sympify(expr)
    var = sp.sympify(var)
    out = []
    # Exact algebraic/rational structure is preferable to table data because it
    # sees cancellation in the concrete expression.
    if denominator_poles is None:
        denominator_poles = _explicit_denominator_poles(expr, var)
    out.extend(denominator_poles)
    if expr.is_Function:
        out.extend(
            _registry_loci(
                expr,
                "Poles",
                SingularityKind.POLE,
                mode="complex",
                var=var,
                assumptions=assumptions,
            )
        )
    # A pole of an argument is not normally a pole of a transcendental outer
    # function, so recurse only through algebraic Add/Mul/Pow containers.
    if expr.is_Add or expr.is_Mul or isinstance(expr, sp.Pow):
        out.extend(
            _recursive_argument_singularities(
                expr, var, SingularityKind.POLE, assumptions
            )
        )
    return _validate_pole_candidates(expr, var, _unique(out))


def _function_branchpoint_records(expr, var, *, assumptions=S.true):
    expr = sp.sympify(expr)
    var = sp.sympify(var)
    out = []
    if expr.is_Function:
        out.extend(
            _registry_loci(
                expr,
                "BranchPoints",
                SingularityKind.BRANCH_POINT,
                mode="complex",
                var=var,
                assumptions=assumptions,
            )
        )
    if isinstance(expr, sp.Pow):
        base, exponent = expr.args
        if exponent.is_integer is not True and base.has(var):
            out.append(
                Singularity(
                    SingularityKind.BRANCH_POINT, sp.Eq(base, 0), S.true, source="power"
                )
            )
    out.extend(
        _recursive_argument_singularities(
            expr, var, SingularityKind.BRANCH_POINT, assumptions
        )
    )
    return _unique(out)


def _pole_location(record, var):
    condition = sp.sympify(record.condition)
    if isinstance(condition, sp.Equality):
        if condition.lhs == var and not condition.rhs.has(var):
            return condition.rhs
        if condition.rhs == var and not condition.lhs.has(var):
            return condition.lhs
    return condition


def poles(expr, var, *, assumptions=S.true, return_result=False):
    """Return poles as ``(location, multiplicity)`` pairs.

    Locations that are not isolated remain exact Boolean conditions.  Use
    ``return_result=True`` for the underlying singularity records and provenance.
    """
    var = sp.sympify(var)
    records = _function_pole_records(expr, var, assumptions=assumptions)
    if return_result:
        return records
    return [(_pole_location(record, var), record.multiplicity) for record in records]


def branchpoints(expr, var, *, assumptions=S.true, return_result=False):
    """Return the Boolean branch-point locus of ``expr``.

    Use ``return_result=True`` for classified records, guards, and provenance.
    """
    records = _function_branchpoint_records(expr, var, assumptions=assumptions)
    if return_result:
        return records
    return _singularity_records_locus(records, assumptions=assumptions)


def definition_cuts(expr, var, *, mode="real", assumptions=S.true):
    expr = sp.sympify(expr)
    var = sp.sympify(var)
    out = []
    if expr.is_Function:
        out.extend(
            _registry_loci(
                expr,
                "DefCuts",
                SingularityKind.DEFINITION_CUT,
                mode=mode,
                var=var,
                assumptions=assumptions,
            )
        )
    # Elementary nonsmooth heads have incomplete DefCuts in parts of the source
    # database, so make their real definition boundaries explicit.
    if mode == "real":
        if expr.func in (sp.Abs, sp.sign) and expr.args and expr.args[0].has(var):
            out.append(
                Singularity(
                    SingularityKind.DEFINITION_CUT,
                    sp.Eq(expr.args[0], 0),
                    S.true,
                    source="structural",
                )
            )
        if expr.func in (sp.floor, sp.ceiling) and expr.args and expr.args[0].has(var):
            a = expr.args[0]
            out.append(
                Singularity(
                    SingularityKind.DEFINITION_CUT,
                    sp.Contains(a, S.Integers),
                    S.true,
                    source="structural",
                )
            )
        if isinstance(expr, sp.Pow):
            base, exponent = expr.args
            if exponent.is_integer is not True and base.has(var):
                out.append(
                    Singularity(
                        SingularityKind.DEFINITION_CUT,
                        sp.Eq(base, 0),
                        S.true,
                        source="power",
                    )
                )
    for arg in getattr(expr, "args", ()):
        if getattr(arg, "has", lambda *_: False)(var):
            out.extend(definition_cuts(arg, var, mode=mode, assumptions=assumptions))
    return _unique(out)


def essential_singularities(expr, var, *, assumptions=S.true):
    expr = sp.sympify(expr)
    var = sp.sympify(var)
    out = []
    if expr.is_Function:
        out.extend(
            _registry_loci(
                expr,
                "Essential",
                SingularityKind.ESSENTIAL,
                mode="complex",
                var=var,
                assumptions=assumptions,
            )
        )
        # A non-polynomial entire transcendental function composed with a pole
        # has an essential singularity at that pole (or, for periodic
        # meromorphic heads, a non-meromorphic accumulation).  These common
        # heads cover the elementary cases that matter for composition.
        transcendental_entire = {
            sp.exp,
            sp.sin,
            sp.cos,
            sp.sinh,
            sp.cosh,
            sp.erf,
            sp.erfc,
        }
        if expr.func in transcendental_entire:
            for arg in expr.args:
                if arg.has(var):
                    for pole in _function_pole_records(
                        arg, var, assumptions=assumptions
                    ):
                        out.append(
                            Singularity(
                                SingularityKind.ESSENTIAL,
                                pole.condition,
                                pole.guard,
                                source="composition",
                            )
                        )
    if expr.is_Add or expr.is_Mul or isinstance(expr, sp.Pow):
        out.extend(
            _recursive_argument_singularities(
                expr, var, SingularityKind.ESSENTIAL, assumptions
            )
        )
    return _unique(out)


def branchcuts(expr, var, *, assumptions=S.true, mode="conditions"):
    """Return principal branch-cut loci, unsolved by default."""
    from ..branchcuts import (
        IMPLEMENTED_RULE_NAMES,
        branch_cut_conditions,
        branch_cut_jumps,
    )

    expr = sp.sympify(expr)
    var = sp.sympify(var)
    if mode == "jumps":
        return branch_cut_jumps(expr, var, assumptions=assumptions)
    out = []
    for cond in branch_cut_conditions(expr, var, assumptions=assumptions):
        locus, guard = _split_locus_guard(cond, var, assumptions=assumptions)
        if locus is not S.false and guard is not S.false:
            out.append(
                Singularity(
                    SingularityKind.BRANCH_CUT, locus, guard, source="branch-rules"
                )
            )
    # Structural jump rules cover the primitive heads for which funcprops can
    # compute boundary values.  Catalogue loci supplement other function heads
    # without duplicating equivalent structural cuts.
    head = getattr(expr.func, "__name__", "").lower()
    if isinstance(expr, sp.Pow):
        head = "power"
    if expr.is_Function and head not in IMPLEMENTED_RULE_NAMES:
        out.extend(
            _registry_loci(
                expr,
                "BranchCuts",
                SingularityKind.BRANCH_CUT,
                mode="complex",
                var=var,
                assumptions=assumptions,
            )
        )

    # Formula-level DLMF structural cuts broader than the available exact jump
    # formulas.  Record the locus without inventing a jump formula at resonant
    # parameter values.
    if expr.func is sp.hyper and len(expr.args) == 3 and expr.args[2] == var:
        aa = tuple(expr.args[0])
        bb = tuple(expr.args[1])
        if len(aa) == len(bb) + 1:
            terminating = (
                Or(*[And(sp.Contains(a, S.Integers), a <= 0) for a in aa])
                if aa
                else S.false
            )
            bottom_bad = (
                Or(*[And(sp.Contains(b, S.Integers), b <= 0) for b in bb])
                if bb
                else S.false
            )
            guard = simplify_condition(
                And(sp.Not(terminating), sp.Not(bottom_bad)), assumptions=assumptions
            )
            if guard is not S.false:
                out.append(
                    Singularity(
                        SingularityKind.BRANCH_CUT,
                        And(sp.Contains(var, S.Reals), var >= 1),
                        guard,
                        source="DLMF 16.2(iii)",
                    )
                )
    if expr.func is sp.assoc_legendre and len(expr.args) == 3 and expr.args[2] == var:
        nu, mu = expr.args[:2]
        polynomial_exception = And(
            sp.Contains(nu, S.Naturals0),
            sp.Contains(mu, S.Integers),
            sp.Abs(mu) <= nu,
            sp.Eq(sp.Mod(mu, 2), 0),
        )
        zero_exception = And(
            sp.Contains(nu, S.Naturals0), sp.Contains(mu, S.Naturals), mu > nu
        )
        guard = simplify_condition(
            And(sp.Not(polynomial_exception), sp.Not(zero_exception)),
            assumptions=assumptions,
        )
        if guard is not S.false:
            out.append(
                Singularity(
                    SingularityKind.BRANCH_CUT,
                    And(sp.Contains(var, S.Reals), var <= 1),
                    guard,
                    source="DLMF 14.21",
                )
            )
    if (
        expr.func in (sp.elliptic_f, sp.elliptic_e)
        and len(expr.args) == 2
        and expr.args[1] == var
    ):
        phi = expr.args[0]
        q = 1 - var * sp.sin(phi) ** 2
        out.append(
            Singularity(
                SingularityKind.BRANCH_CUT,
                And(sp.Eq(sp.im(q), 0), sp.re(q) <= 0),
                S.true,
                source="DLMF 19.2.4-19.2.5",
            )
        )
    if expr.func is sp.elliptic_pi and len(expr.args) == 3:
        _n, phi, m = expr.args
        if m == var:
            q = 1 - var * sp.sin(phi) ** 2
            out.append(
                Singularity(
                    SingularityKind.BRANCH_CUT,
                    And(sp.Eq(sp.im(q), 0), sp.re(q) <= 0),
                    S.true,
                    source="DLMF 19.2.7",
                )
            )
    out = _unique(out)
    conds = tuple(
        s.condition if s.guard is S.true else And(s.condition, s.guard) for s in out
    )
    if mode == "conditions":
        return conds
    if mode == "records":
        return out
    if mode == "union":
        return Or(*conds) if conds else S.false
    raise ValueError("mode must be 'conditions', 'records', 'union', or 'jumps'")


def _piecewise_boundaries(expr):
    out = []
    if not isinstance(expr, sp.Piecewise):
        return out
    for _, cond in expr.args:
        if cond in (S.true, True):
            continue
        for rel in cond.atoms(Relational) if hasattr(cond, "atoms") else ():
            out.append(sp.Eq(rel.lhs, rel.rhs))
    return out


def _piecewise_branch_for_side(expr, var, point, side):
    t = sp.Dummy("funcprops_side", positive=True)
    probe = point + t if side == "+" else point - t
    for value, cond in expr.args:
        if cond in (S.true, True):
            return value
        try:
            test = sp.simplify(sp.refine(cond.subs(var, probe), sp.Q.positive(t)))
        except SYMBOLIC_QUERY_ERRORS:
            test = cond.subs(var, probe)
        if test is S.true or test is True:
            return value
        if test is S.false or test is False:
            continue
        # reduce simple one-variable inequalities in the positive probe.
        try:
            red = sp.reduce_inequalities(test, t)
            if red is S.true:
                return value
            if red is S.false:
                continue
        except SYMBOLIC_QUERY_ERRORS:
            pass
        return None
    return None


def _piecewise_discontinuous_at(expr, var, condition):
    """Return True/False/None for a simple explicitly solvable boundary."""
    try:
        sols = sp.solve(condition, var)
    except SYMBOLIC_QUERY_ERRORS:
        return None
    if not sols or any(s.has(var) for s in sols):
        return None
    any_unknown = False
    for p in sols:
        left_branch = _piecewise_branch_for_side(expr, var, p, "-")
        right_branch = _piecewise_branch_for_side(expr, var, p, "+")
        if left_branch is None or right_branch is None:
            any_unknown = True
            continue
        try:
            left = (
                sp.limit(left_branch, var, p, dir="-")
                if left_branch.has(var)
                else left_branch
            )
            right = (
                sp.limit(right_branch, var, p, dir="+")
                if right_branch.has(var)
                else right_branch
            )
            value = expr.subs(var, p)
            side_match = certify_identity_equal(left, right)
            value_match = certify_identity_equal(value, left)
            if side_match is False or value_match is False:
                return True
            if side_match is None or value_match is None:
                any_unknown = True
        except SYMBOLIC_QUERY_ERRORS:
            any_unknown = True
    return None if any_unknown else False


def _registry_discontinuity_records(
    expr, var, *, assumptions=S.true, registry_var=None
):
    """Return registry discontinuities for scalar or joint analysis variables.

    ``registry_var`` is the distinguished scalar variable used by
    univariate catalogue rules.  Joint multivariate transport leaves it unset while still
    using ``var`` to separate analysis geometry from parameter guards.
    """
    query_var = (
        var if registry_var is None and len(_variable_tuple(var)) == 1 else registry_var
    )
    result = evaluate_property(
        expr, "Discontinuities", var=query_var, assumptions=assumptions
    )
    out = []
    for case in result.cases:
        value = case.value
        if _is_empty(value):
            continue
        entries = (
            tuple(value) if isinstance(value, (sp.Tuple, tuple, list)) else (value,)
        )
        for entry in entries:
            condition = (
                entry[0]
                if isinstance(entry, (sp.Tuple, tuple, list)) and entry
                else entry
            )
            locus, local_guard = _split_locus_guard(
                condition, var, assumptions=assumptions
            )
            guard = _combine_guard(case.condition, local_guard, assumptions=assumptions)
            if locus is S.false or guard is S.false:
                continue
            out.append(
                Singularity(
                    SingularityKind.DISCONTINUITY,
                    sp.sympify(locus),
                    sp.sympify(guard),
                    source="registry",
                )
            )
    return _unique(out)


def _function_discontinuity_records(expr, var, *, assumptions=S.true):
    """Return unsolved real discontinuity conditions."""
    expr = sp.sympify(expr)
    var = sp.sympify(var)
    out = []
    if expr.is_Function:
        out.extend(_registry_discontinuity_records(expr, var, assumptions=assumptions))
    if isinstance(expr, sp.Piecewise):
        from ..real_analysis.piecewise import (
            certify_piecewise_boundaries,
            piecewise_decomposition,
        )

        try:
            decomposition = piecewise_decomposition(expr, var, assumptions)
            certificates = certify_piecewise_boundaries(decomposition, max_order=0)
        except SYMBOLIC_QUERY_ERRORS:
            certificates = ()
        if certificates:
            for certificate in certificates:
                if certificate.continuous is not True:
                    out.append(
                        Singularity(
                            SingularityKind.DISCONTINUITY,
                            sp.Eq(var, certificate.boundary.point),
                            S.true,
                            source="piecewise-boundary-certificate",
                        )
                    )
        else:
            for boundary in _piecewise_boundaries(expr):
                verdict = _piecewise_discontinuous_at(expr, var, boundary)
                if verdict is not False:
                    out.append(
                        Singularity(
                            SingularityKind.DISCONTINUITY,
                            boundary,
                            S.true,
                            source="piecewise",
                        )
                    )
    # Discontinuities of inner arguments propagate through ordinary composition
    # unless an outer rule proves a cancellation; conservatively retain them.
    recursive_args = (
        [v for v, _ in expr.args]
        if isinstance(expr, sp.Piecewise)
        else getattr(expr, "args", ())
    )
    for arg in recursive_args:
        if getattr(arg, "has", lambda *_: False)(var):
            out.extend(
                _function_discontinuity_records(arg, var, assumptions=assumptions)
            )
    return _unique(out)


_SINGULARITY_LOCUS_SIMPLIFY_OP_LIMIT = 12


def _bounded_locus_simplify(condition, *, assumptions=S.true):
    """Simplify small singularity loci without exploding warmed Boolean states."""
    condition = sp.sympify(condition)
    assumptions = sp.sympify(assumptions)
    if condition in (S.true, S.false):
        return condition
    if assumptions is not S.true:
        condition = And(condition, assumptions)
    # Rebuild top-level Boolean nodes structurally so an already-expanded warmed
    # expression is not sent through the global simplifier again.
    if condition.func is And:
        condition = And(*condition.args, evaluate=False)
    elif condition.func is Or:
        condition = Or(*condition.args, evaluate=False)
    if sp.count_ops(condition) > _SINGULARITY_LOCUS_SIMPLIFY_OP_LIMIT:
        return condition
    return simplify_condition(condition)


def discontinuities(expr, var, *, assumptions=S.true, return_result=False):
    """Return the Boolean real-discontinuity locus of ``expr``.

    Use ``return_result=True`` for classified records and proof provenance.
    """
    records = _function_discontinuity_records(expr, var, assumptions=assumptions)
    if return_result:
        return records
    return _singularity_records_locus(records, assumptions=assumptions)


def _singularity_record_locus(record, *, assumptions=S.true):
    return _bounded_locus_simplify(
        And(record.condition, record.guard), assumptions=assumptions
    )


def _singularity_records_locus(records, *, assumptions=S.true):
    conditions = tuple(
        _singularity_record_locus(record, assumptions=assumptions) for record in records
    )
    conditions = tuple(
        condition for condition in conditions if condition is not S.false
    )
    if not conditions:
        return S.false
    return _bounded_locus_simplify(Or(*conditions), assumptions=assumptions)


def _restrict_boolean_locus(locus, restriction, var, *, assumptions=S.true):
    """Intersect a Boolean locus with a function restriction conservatively."""
    locus = sp.sympify(locus)
    restriction = simplify_condition(restriction, assumptions=assumptions)
    if restriction is S.false or locus is S.false:
        return S.false
    if restriction is S.true:
        if len(_variable_tuple(var)) == 1:
            return simplify_condition(locus, assumptions=assumptions)
        return locus
    locus = simplify_condition(locus, assumptions=assumptions)
    combined = simplify_condition(And(locus, restriction), assumptions=assumptions)
    if combined in (S.true, S.false):
        return combined

    # ``simplify_logic`` does not reason about inequalities across DNF branches.
    # Resolve one-variable real clauses separately so restrictions such as x > 0
    # can certify that the principal logarithm's nonpositive-real cut is absent.
    vv = _variable_tuple(var)
    if len(vv) != 1:
        try:
            from .._semialg import SemialgUnavailable, is_satisfiable

            if not is_satisfiable(combined, vv):
                return S.false
        except (SemialgUnavailable, *SYMBOLIC_QUERY_ERRORS):
            pass
        return combined
    variable = vv[0]
    try:
        dnf = sp.to_dnf(combined, simplify=True, force=True)
    except SYMBOLIC_QUERY_ERRORS:
        return combined
    clauses = dnf.args if dnf.func is Or else (dnf,)
    reduced = []
    for clause in clauses:
        atoms = clause.args if clause.func is And else (clause,)
        contains_reals = any(
            isinstance(atom, sp.Contains) and atom.args == (variable, S.Reals)
            for atom in atoms
        )
        inequalities = tuple(
            atom
            for atom in atoms
            if atom.has(variable) and isinstance(atom, Relational)
        )
        residual = tuple(
            atom
            for atom in atoms
            if atom not in inequalities
            and not (isinstance(atom, sp.Contains) and atom.args == (variable, S.Reals))
        )
        if inequalities:
            try:
                real_part = sp.reduce_inequalities(inequalities, variable)
            except SYMBOLIC_QUERY_ERRORS:
                real_part = And(*inequalities)
            if real_part is S.false:
                continue
            clause_result = And(real_part, *residual)
            if contains_reals and not inequalities:
                clause_result = And(sp.Contains(variable, S.Reals), clause_result)
        else:
            clause_result = clause
        clause_result = simplify_condition(clause_result, assumptions=assumptions)
        if clause_result is not S.false:
            reduced.append(clause_result)
    if not reduced:
        return S.false
    return simplify_condition(Or(*reduced), assumptions=assumptions)


def _merge_zero_boundary_relations(condition):
    """Merge ``g == 0`` with an adjacent strict half-line exactly."""
    condition = sp.sympify(condition)
    if condition.func not in (And, Or):
        return condition
    args = [_merge_zero_boundary_relations(arg) for arg in condition.args]
    if condition.func is not Or:
        return And(*args)

    for index, equality in enumerate(args):
        if not isinstance(equality, sp.Equality):
            continue
        delta = sp.expand(equality.lhs - equality.rhs)
        for other_index, relation in enumerate(args):
            if index == other_index:
                continue
            relation_delta = (
                sp.expand(relation.lhs - relation.rhs)
                if isinstance(relation, Relational)
                else None
            )
            replacement = None
            if relation_delta == delta:
                if isinstance(relation, sp.StrictLessThan):
                    replacement = sp.LessThan(delta, 0)
                elif isinstance(relation, sp.StrictGreaterThan):
                    replacement = sp.GreaterThan(delta, 0)
            elif relation_delta == -delta:
                if isinstance(relation, sp.StrictLessThan):
                    replacement = sp.GreaterThan(delta, 0)
                elif isinstance(relation, sp.StrictGreaterThan):
                    replacement = sp.LessThan(delta, 0)
            if replacement is not None:
                remaining = [
                    arg
                    for pos, arg in enumerate(args)
                    if pos not in (index, other_index)
                ]
                return Or(replacement, *remaining)
    return Or(*args)


def _explicit_complex_locus(locus, var):
    """Write simple complex loci explicitly in ``re(var)``, ``im(var)``.

    Ordering relations on a complex analysis variable describe subsets of the
    real axis.  Preserve that geometry in public complex-domain results rather
    than allowing simplification to collapse it to an apparently real-only
    relation such as ``x <= 0``.
    """
    locus = sp.sympify(locus)
    if locus in (S.true, S.false):
        return locus
    vv = _variable_tuple(var)
    if len(vv) != 1:
        return locus
    variable = vv[0]

    def rewrite(node):
        if isinstance(node, sp.Contains) and node.args == (variable, S.Reals):
            return sp.Eq(sp.im(variable), 0)
        if node.func in (And, Or):
            return node.func(*(rewrite(arg) for arg in node.args))
        if isinstance(node, Relational):
            lhs, rhs = node.lhs, node.rhs
            if lhs == variable and not rhs.has(variable) and rhs.is_real is True:
                real_relation = node.func(sp.re(variable), rhs)
                return And(sp.Eq(sp.im(variable), 0), real_relation)
            if rhs == variable and not lhs.has(variable) and lhs.is_real is True:
                real_relation = node.func(lhs, sp.re(variable))
                return And(sp.Eq(sp.im(variable), 0), real_relation)
        return node

    rewritten = _merge_zero_boundary_relations(rewrite(locus))
    simplified = sp.simplify_logic(rewritten, force=True)
    return _merge_zero_boundary_relations(simplified)
