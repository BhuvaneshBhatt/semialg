"""Internal singularity analysis implementation."""

from __future__ import annotations

from functools import lru_cache

import sympy as sp
from sympy import S
from sympy.core.relational import Relational
from sympy.logic.boolalg import And, Or

from ..engine import simplify_condition
from ..exceptions import SYMBOLIC_QUERY_ERRORS
from ..results import (
    Singularity,
    SingularityKind,
)
from .records import (
    _denominator_singularity_records,
    _function_branchpoint_records,
    _function_discontinuity_records,
    _function_pole_records,
    _merge_zero_boundary_relations,
    _registry_discontinuity_records,
    _registry_loci,
    _unique,
    branchcuts,
    definition_cuts,
    essential_singularities,
    order_at,
    removable_singularities,
)


def _depends_on_any(expr, variables):
    return any(getattr(expr, "has", lambda *_: False)(v) for v in variables)


def _pullback_singularity_record(record, dummy, argument, *, source):
    """Pull one univariate singularity record back through ``dummy = argument``."""
    condition = sp.sympify(record.condition).subs(dummy, argument)
    guard = sp.sympify(record.guard).subs(dummy, argument)
    multiplicity = record.multiplicity
    if isinstance(multiplicity, sp.Basic):
        multiplicity = multiplicity.subs(dummy, argument)
    return Singularity(
        record.kind,
        condition,
        guard,
        multiplicity,
        source=source if record.source is None else f"{source}:{record.source}",
    )


def _multivariate_denominator_records(expr, variables):
    """Classify explicit denominator hypersurfaces without assuming isolation.

    When a denominator factor can itself be used as a local coordinate, replacing
    that factor by a dummy variable reduces the transverse behavior to the univariate
    local-order test.  This certifies common removable hypersurfaces such as
    ``sin(x + y)/(x + y)`` while retaining an unsolved pole condition when no
    transverse order can be proved.
    """
    vv = tuple(variables)
    expr = sp.sympify(expr)
    try:
        _num, den = sp.fraction(sp.together(expr))
    except SYMBOLIC_QUERY_ERRORS:
        return ()
    if den == 1 or not _depends_on_any(den, vv):
        return ()

    factors = []
    try:
        _constant, factor_data = sp.factor_list(den)
        factors = [(sp.sympify(factor), mult) for factor, mult in factor_data]
    except SYMBOLIC_QUERY_ERRORS:
        factors = [(den, None)]
    if not factors:
        factors = [(den, None)]

    out = []
    for factor, denominator_mult in factors:
        if not _depends_on_any(factor, vv):
            continue
        dummy = sp.Symbol("_funcprops_transverse")
        try:
            pulled = expr.xreplace({factor: dummy})
            order = order_at(pulled, dummy, S.Zero)
        except SYMBOLIC_QUERY_ERRORS:
            order = None
        condition = sp.Eq(factor, 0)
        if order is not None and order.is_nonnegative:
            out.append(
                Singularity(
                    SingularityKind.REMOVABLE,
                    condition,
                    S.true,
                    source="multivariate-transverse-order",
                )
            )
            continue
        if order is not None and order.is_negative:
            multiplicity = -order
            if multiplicity.is_Integer:
                multiplicity = int(multiplicity)
            out.append(
                Singularity(
                    SingularityKind.POLE,
                    condition,
                    S.true,
                    multiplicity,
                    source="multivariate-transverse-order",
                )
            )
            continue
        # An unresolved denominator hypersurface is not unconditionally a pole:
        # the numerator may vanish on it.  Certify only the non-cancelling stratum.
        guard = S.true
        if denominator_mult is not None:
            try:
                scaled = sp.cancel(pulled * dummy**denominator_mult)
                coefficient = sp.simplify(scaled.subs(dummy, S.Zero))
                if coefficient.has(dummy):
                    coefficient = sp.simplify(sp.limit(scaled, dummy, S.Zero))
                if coefficient == 0:
                    continue
                if coefficient.is_zero is not False:
                    guard = sp.Ne(coefficient, 0)
            except SYMBOLIC_QUERY_ERRORS:
                # Without a transverse coefficient there is no sound pole proof.
                continue
        else:
            continue
        out.append(
            Singularity(
                SingularityKind.POLE,
                condition,
                guard,
                denominator_mult,
                source="multivariate-denominator-noncancellation",
            )
        )
    return _unique(out)


@lru_cache(maxsize=128)
def _canonical_multivar_records(expr, canonical_variables, complex_mode, assumptions):
    """Transport joint registry singularity facts into multivariate loci.

    One-coordinate function heads are handled by the univariate pullback path.
    This routine handles the complementary case: native registry rules for a function
    whose singularity geometry depends on two or more varying arguments.  The
    catalogue rule is instantiated directly on the concrete node, so joint
    conditions such as ``m*sin(phi)**2 == 1`` retain their full geometry and
    parameter guards without introducing an artificial scalar coordinate.
    """
    vv = tuple(canonical_variables)
    analysis_var = sp.Tuple(*vv)
    out = []
    properties = [
        ("Poles", SingularityKind.POLE, "complex"),
        ("Essential", SingularityKind.ESSENTIAL, "complex"),
        ("BranchPoints", SingularityKind.BRANCH_POINT, "complex"),
        (
            "DefCuts",
            SingularityKind.DEFINITION_CUT,
            "complex" if complex_mode else "real",
        ),
    ]
    if complex_mode:
        properties.append(("BranchCuts", SingularityKind.BRANCH_CUT, "complex"))

    for node in sp.preorder_traversal(sp.sympify(expr)):
        if not getattr(node, "is_Function", False):
            continue
        varying = [argument for argument in node.args if _depends_on_any(argument, vv)]
        if len(varying) < 2:
            continue
        for prop, kind, mode in properties:
            out.extend(
                _registry_loci(
                    node,
                    prop,
                    kind,
                    mode=mode,
                    var=analysis_var,
                    assumptions=assumptions,
                )
            )
        if not complex_mode:
            out.extend(
                _registry_discontinuity_records(
                    node,
                    analysis_var,
                    assumptions=assumptions,
                    registry_var=None,
                )
            )
    return _unique(out)


def _multivar_pullbacks(expr, variables, *, complex_mode, assumptions=S.true):
    """Reuse the univariate singularity classifiers on one-coordinate pullbacks.

    A multivariate function node with exactly one varying argument is a pullback
    of a univariate function along that argument.  Running the established univariate
    classifiers on a dummy coordinate keeps the two paths classification-
    symmetric, including registry-backed special-function facts.
    """
    vv = tuple(variables)
    out = []
    elementary_structural_heads = {
        sp.log,
        sp.exp,
        sp.sin,
        sp.cos,
        sp.sinh,
        sp.cosh,
        sp.erf,
        sp.erfc,
        sp.Abs,
        sp.sign,
        sp.floor,
        sp.ceiling,
    }
    for node in sp.preorder_traversal(expr):
        if not getattr(node, "is_Function", False):
            continue
        if node.func in elementary_structural_heads:
            # These heads are handled by explicit multivariate rules below.
            # Avoid running the full univariate registry pipeline at every elementary
            # node; that would make simple multivariate expressions needlessly
            # expensive while adding no classification information.
            continue
        varying = [arg for arg in node.args if _depends_on_any(arg, vv)]
        if len(varying) != 1:
            continue
        argument = varying[0]
        dummy = sp.Symbol("_funcprops_pullback")
        try:
            template = node.xreplace({argument: dummy})
        except SYMBOLIC_QUERY_ERRORS:
            continue
        if template == node or template.has(*vv):
            continue

        groups = [
            _function_pole_records(template, dummy, assumptions=assumptions),
            removable_singularities(template, dummy),
            essential_singularities(template, dummy, assumptions=assumptions),
            _function_branchpoint_records(template, dummy, assumptions=assumptions),
            definition_cuts(
                template,
                dummy,
                mode="complex" if complex_mode else "real",
                assumptions=assumptions,
            ),
        ]
        if complex_mode:
            groups.append(
                branchcuts(template, dummy, assumptions=assumptions, mode="records")
            )
        else:
            groups.append(
                _function_discontinuity_records(
                    template, dummy, assumptions=assumptions
                )
            )

        for records in groups:
            for record in records:
                out.append(
                    _pullback_singularity_record(
                        record,
                        dummy,
                        argument,
                        source="multivariate-pullback",
                    )
                )
    return _unique(out)


def _multivariate_registry_records(
    expr, variables, *, complex_mode, assumptions=S.true
):
    """Alpha-normalized cached transport of joint multivariate registry facts."""
    vv = tuple(variables)
    canonical = tuple(sp.Symbol(f"_funcprops_alpha_{i}") for i in range(len(vv)))
    to_canonical = dict(zip(vv, canonical, strict=True))
    from_canonical = dict(zip(canonical, vv, strict=True))
    canonical_expr = sp.sympify(expr).xreplace(to_canonical)
    canonical_assumptions = sp.sympify(assumptions).xreplace(to_canonical)
    records = _canonical_multivar_records(
        canonical_expr, canonical, bool(complex_mode), canonical_assumptions
    )
    out = []
    for record in records:
        multiplicity = record.multiplicity
        if isinstance(multiplicity, sp.Basic):
            multiplicity = multiplicity.xreplace(from_canonical)
        out.append(
            Singularity(
                record.kind,
                sp.sympify(record.condition).xreplace(from_canonical),
                sp.sympify(record.guard).xreplace(from_canonical),
                multiplicity,
                record.source,
            )
        )
    return tuple(out)


def _multivar_structural_singularities(
    expr, variables, *, complex_mode, assumptions=S.true
):
    """Return classified exact singularity loci for multivariate expressions.

    The multivariate path deliberately mirrors all supported singularity classes.  It uses
    exact denominator hypersurfaces, univariate pullbacks for function heads,
    and structural composition rules where singular behavior comes from an
    argument.
    """
    vv = tuple(variables)
    expr = sp.sympify(expr)
    out = list(_multivariate_denominator_records(expr, vv))
    out.extend(
        _multivariate_registry_records(
            expr, vv, complex_mode=complex_mode, assumptions=assumptions
        )
    )
    out.extend(
        _multivar_pullbacks(
            expr, vv, complex_mode=complex_mode, assumptions=assumptions
        )
    )

    transcendental_entire = {
        sp.exp,
        sp.sin,
        sp.cos,
        sp.sinh,
        sp.cosh,
        sp.erf,
        sp.erfc,
    }
    for node in sp.preorder_traversal(expr):
        args = getattr(node, "args", ())
        if not complex_mode and args and _depends_on_any(args[0], vv):
            argument = args[0]
            if node.func in (sp.Abs, sp.sign):
                out.append(
                    Singularity(
                        SingularityKind.DEFINITION_CUT,
                        sp.Eq(argument, 0),
                        S.true,
                        source="multivariate-structural",
                    )
                )
                if node.func is sp.sign:
                    out.append(
                        Singularity(
                            SingularityKind.DISCONTINUITY,
                            sp.Eq(argument, 0),
                            S.true,
                            source="multivariate-structural",
                        )
                    )
            elif node.func in (sp.floor, sp.ceiling):
                integer_locus = sp.Contains(argument, S.Integers)
                out.append(
                    Singularity(
                        SingularityKind.DEFINITION_CUT,
                        integer_locus,
                        S.true,
                        source="multivariate-structural",
                    )
                )
                out.append(
                    Singularity(
                        SingularityKind.DISCONTINUITY,
                        integer_locus,
                        S.true,
                        source="multivariate-structural",
                    )
                )

        if node.func in transcendental_entire:
            for argument in args:
                if not _depends_on_any(argument, vv):
                    continue
                for pole in _multivariate_denominator_records(argument, vv):
                    if pole.kind is SingularityKind.POLE:
                        out.append(
                            Singularity(
                                SingularityKind.ESSENTIAL,
                                pole.condition,
                                pole.guard,
                                source="multivariate-composition",
                            )
                        )

        if node.func is sp.log and args and _depends_on_any(args[0], vv):
            argument = args[0]
            out.append(
                Singularity(
                    SingularityKind.BRANCH_POINT,
                    sp.Eq(argument, 0),
                    S.true,
                    source="multivariate-log",
                )
            )
            if complex_mode:
                out.append(
                    Singularity(
                        SingularityKind.BRANCH_CUT,
                        And(sp.Eq(sp.im(argument), 0), sp.re(argument) < 0),
                        S.true,
                        source="multivariate-log",
                    )
                )
            continue
        if isinstance(node, sp.Pow):
            base, exponent = node.args
            if exponent.is_integer is not True and _depends_on_any(base, vv):
                out.append(
                    Singularity(
                        SingularityKind.BRANCH_POINT,
                        sp.Eq(base, 0),
                        S.true,
                        source="multivariate-power",
                    )
                )
                if complex_mode:
                    out.append(
                        Singularity(
                            SingularityKind.BRANCH_CUT,
                            And(sp.Eq(sp.im(base), 0), sp.re(base) < 0),
                            S.true,
                            source="multivariate-power",
                        )
                    )
                else:
                    out.append(
                        Singularity(
                            SingularityKind.DEFINITION_CUT,
                            sp.Eq(base, 0),
                            S.true,
                            source="multivariate-power",
                        )
                    )

    return _unique(out)


def _relation_closure(atom):
    if isinstance(atom, sp.StrictGreaterThan):
        return sp.GreaterThan(atom.lhs, atom.rhs)
    if isinstance(atom, sp.StrictLessThan):
        return sp.LessThan(atom.lhs, atom.rhs)
    return atom


def _relation_interior(atom):
    if isinstance(atom, sp.GreaterThan):
        return sp.StrictGreaterThan(atom.lhs, atom.rhs)
    if isinstance(atom, sp.LessThan):
        return sp.StrictLessThan(atom.lhs, atom.rhs)
    if isinstance(atom, sp.Equality):
        return S.false
    return atom


def _domain_geometry(condition, variables, *, assumptions=S.true):
    """Return certified ambient interior, boundary and exterior conditions.

    This handles the semialgebraic support forms used throughout funcprops:
    conjunctions of equalities and weak/strict inequalities, plus simple
    disjunctions.  Unsupported Boolean atoms are retained conservatively.
    """
    condition = simplify_condition(condition, assumptions=assumptions)
    if condition is S.true:
        return S.true, S.false, S.false
    if condition is S.false:
        return S.false, S.false, S.true

    def rec(node):
        if node.func is And:
            parts = [rec(a) for a in node.args]
            closure = And(*(p[1] for p in parts))
            interior = And(*(p[0] for p in parts))
            return interior, closure
        if node.func is Or:
            parts = [rec(a) for a in node.args]
            interior = Or(*(p[0] for p in parts))
            closure = Or(*(p[1] for p in parts))
            return interior, closure
        if isinstance(node, Relational):
            return _relation_interior(node), _relation_closure(node)
        return node, node

    interior, closure = rec(condition)
    interior = simplify_condition(interior, assumptions=assumptions)
    closure = simplify_condition(closure, assumptions=assumptions)
    boundary = simplify_condition(
        And(closure, sp.Not(interior)), assumptions=assumptions
    )
    exterior = simplify_condition(sp.Not(closure), assumptions=assumptions)
    return interior, boundary, exterior


def _complexify_equality_locus_impl(locus, variables):
    """Expand complex equalities into explicit real/imaginary equations."""
    vv = tuple(variables)

    def rewrite(node):
        if node.func in (And, Or):
            rewritten = node.func(*(rewrite(a) for a in node.args))
            if rewritten.func is Or and len(rewritten.args) == 2:
                a, b = rewritten.args
                pairs = ((a, b), (b, a))
                for equality, relation in pairs:
                    if not isinstance(equality, sp.Equality):
                        continue
                    delta = sp.expand(equality.lhs - equality.rhs)
                    if delta == 0:
                        continue
                    if isinstance(relation, sp.StrictLessThan):
                        if sp.expand(relation.lhs - relation.rhs) == delta:
                            return sp.LessThan(delta, 0)
                    if isinstance(relation, sp.StrictGreaterThan):
                        if sp.expand(relation.lhs - relation.rhs) == delta:
                            return sp.GreaterThan(delta, 0)
            return rewritten
        if isinstance(node, sp.Equality):
            delta = sp.expand(node.lhs - node.rhs)
            if any(delta.has(v) for v in vv):
                # Registry branch-cut rules are frequently already expressed in
                # explicit real/imaginary coordinates.  Re-complexifying an
                # equality such as ``im(g(x, y)) == 0`` is redundant and can
                # trigger very large expansions for special functions.
                if delta.is_real is True:
                    return sp.Eq(delta, 0)
                if sp.count_ops(delta) > 8 or delta.has(
                    sp.sin,
                    sp.cos,
                    sp.tan,
                    sp.sinh,
                    sp.cosh,
                    sp.tanh,
                ):
                    real_part = sp.re(delta, evaluate=False)
                    imag_part = sp.im(delta, evaluate=False)
                else:
                    real_part = sp.re(delta)
                    imag_part = sp.im(delta)
                return And(sp.Eq(real_part, 0), sp.Eq(imag_part, 0))
        return node

    locus = sp.sympify(locus)
    rewritten = rewrite(locus)
    # Small elementary loci benefit from Boolean normalization (for example,
    # the log branch point plus negative-real cut becomes one nonpositive-real
    # condition).  Large special-function preimages can expand exponentially
    # under simplify_logic; their structurally exact Boolean form is already a
    # canonical public answer, so avoid distributive normalization there.
    if sp.count_ops(rewritten) > 30:
        return rewritten
    rewritten = _merge_zero_boundary_relations(rewritten)
    result = simplify_condition(rewritten)
    return _merge_zero_boundary_relations(simplify_condition(rewrite(result)))


@lru_cache(maxsize=256)
def _cached_complexify_equality_locus(locus, canonical_variables):
    return _complexify_equality_locus_impl(locus, canonical_variables)


def _complexify_equality_locus(locus, variables):
    """Alpha-normalized cached complexification of multivariate equality loci."""
    vv = tuple(variables)
    canonical = tuple(sp.Symbol(f"_funcprops_alpha_{i}") for i in range(len(vv)))
    to_canonical = dict(zip(vv, canonical, strict=True))
    from_canonical = dict(zip(canonical, vv, strict=True))
    canonical_locus = sp.sympify(locus).xreplace(to_canonical)
    result = _cached_complexify_equality_locus(canonical_locus, canonical)
    return sp.sympify(result).xreplace(from_canonical)


def _collect_singularity_records(expr, variables, *, complex_mode, assumptions):
    if len(variables) == 1:
        variable = variables[0]
        denominator_poles, denominator_removable = _denominator_singularity_records(
            expr, variable
        )
        return {
            "poles": _function_pole_records(
                expr,
                variable,
                assumptions=assumptions,
                denominator_poles=denominator_poles,
            ),
            "removable": denominator_removable,
            "essential": essential_singularities(
                expr, variable, assumptions=assumptions
            ),
            "branch_points": _function_branchpoint_records(
                expr, variable, assumptions=assumptions
            ),
            "branch_cuts": branchcuts(
                expr, variable, assumptions=assumptions, mode="records"
            )
            if complex_mode
            else (),
            "definition_cuts": definition_cuts(
                expr,
                variable,
                mode="complex" if complex_mode else "real",
                assumptions=assumptions,
            ),
            "discontinuities": ()
            if complex_mode
            else _function_discontinuity_records(
                expr, variable, assumptions=assumptions
            ),
        }

    structural = _multivar_structural_singularities(
        expr, variables, complex_mode=complex_mode, assumptions=assumptions
    )
    return {
        "poles": tuple(r for r in structural if r.kind is SingularityKind.POLE),
        "removable": tuple(
            r for r in structural if r.kind is SingularityKind.REMOVABLE
        ),
        "essential": tuple(
            r for r in structural if r.kind is SingularityKind.ESSENTIAL
        ),
        "branch_points": tuple(
            r for r in structural if r.kind is SingularityKind.BRANCH_POINT
        ),
        "branch_cuts": tuple(
            r for r in structural if r.kind is SingularityKind.BRANCH_CUT
        )
        if complex_mode
        else (),
        "definition_cuts": tuple(
            r for r in structural if r.kind is SingularityKind.DEFINITION_CUT
        ),
        "discontinuities": tuple(
            r for r in structural if r.kind is SingularityKind.DISCONTINUITY
        )
        if not complex_mode
        else (),
    }
