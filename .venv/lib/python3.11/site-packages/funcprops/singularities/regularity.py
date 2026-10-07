"""Internal singularity analysis implementation."""

from __future__ import annotations

import sympy as sp
from sympy import S
from sympy.calculus.util import continuous_domain
from sympy.core.relational import Relational
from sympy.logic.boolalg import And, Boolean, Or

from ..engine import evaluate_property, simplify_condition
from ..exceptions import SYMBOLIC_QUERY_ERRORS
from ..results import (
    Analyticity,
    ConditionalExpression,
    Continuity,
    Meromorphicity,
    PropertyResult,
    Singularity,
    SingularityKind,
    SingularityResult,
)
from ..stratification import evaluate_parameter_strata
from .geometry import (
    _collect_singularity_records,
    _complexify_equality_locus,
    _domain_geometry,
)
from .records import (
    _bounded_locus_simplify,
    _explicit_complex_locus,
    _function_branchpoint_records,
    _function_discontinuity_records,
    _merge_zero_boundary_relations,
    _piecewise_boundaries,
    _restrict_boolean_locus,
    _singularity_record_locus,
    _singularity_records_locus,
    _variable_tuple,
    branchcuts,
    definition_cuts,
    essential_singularities,
)


def singularities(
    expr,
    var,
    domain=S.Complexes,
    *,
    assumptions=S.true,
    return_result=False,
):
    """Return the Boolean locus where ``expr`` is singular.

    For several analysis variables, singularities remain as exact Boolean loci
    such as ``Eq(x + y, 0)``.
    ``ConditionalExpression`` restrictions define the effective function domain.
    Detailed results additionally classify the unrestricted singular locus as
    interior, boundary, or exterior to that domain.
    """
    restriction = S.true
    while isinstance(expr, ConditionalExpression):
        restriction = simplify_condition(
            And(restriction, expr.condition), assumptions=assumptions
        )
        expr = expr.value
    expr = sp.sympify(expr)
    variables = _variable_tuple(var)
    analysis_var = variables[0] if len(variables) == 1 else sp.Tuple(*variables)
    complex_mode = domain in (S.Complexes, "complex", "complexes", "C", "CC")

    classified_records = _collect_singularity_records(
        expr, variables, complex_mode=complex_mode, assumptions=assumptions
    )

    unrestricted_loci = {
        name: _singularity_records_locus(records, assumptions=assumptions)
        for name, records in classified_records.items()
    }
    unrestricted = _bounded_locus_simplify(
        Or(*unrestricted_loci.values()), assumptions=assumptions
    )

    classified_loci = {
        name: _restrict_boolean_locus(
            locus, restriction, analysis_var, assumptions=assumptions
        )
        for name, locus in unrestricted_loci.items()
    }
    # Intersection distributes over the classified union.  Reuse the already
    # restricted component loci
    # second time.
    singularities = _bounded_locus_simplify(
        Or(*classified_loci.values()), assumptions=assumptions
    )

    domain_interior, domain_boundary, domain_exterior = _domain_geometry(
        restriction, variables, assumptions=assumptions
    )
    interior_singularities = _restrict_boolean_locus(
        unrestricted, domain_interior, analysis_var, assumptions=assumptions
    )
    boundary_singularities = _restrict_boolean_locus(
        unrestricted, domain_boundary, analysis_var, assumptions=assumptions
    )
    exterior_singularities = _restrict_boolean_locus(
        unrestricted, domain_exterior, analysis_var, assumptions=assumptions
    )

    if complex_mode:
        if len(variables) == 1:
            singularities = _explicit_complex_locus(singularities, variables[0])
            unrestricted = _explicit_complex_locus(unrestricted, variables[0])
            classified_loci = {
                name: _explicit_complex_locus(locus, variables[0])
                for name, locus in classified_loci.items()
            }
            interior_singularities = _explicit_complex_locus(
                interior_singularities, variables[0]
            )
            boundary_singularities = _explicit_complex_locus(
                boundary_singularities, variables[0]
            )
            exterior_singularities = _explicit_complex_locus(
                exterior_singularities, variables[0]
            )
        else:
            complex_locus_cache = {}

            def complexify(locus):
                key = sp.srepr(locus)
                if key not in complex_locus_cache:
                    complex_locus_cache[key] = _complexify_equality_locus(
                        locus, variables
                    )
                return complex_locus_cache[key]

            classified_loci = {
                name: complexify(locus) for name, locus in classified_loci.items()
            }
            if restriction is S.true:
                aggregate = Or(*classified_loci.values())
                unrestricted = (
                    _merge_zero_boundary_relations(simplify_condition(aggregate))
                    if sp.count_ops(aggregate) <= 80
                    else aggregate
                )
                singularities = unrestricted
                interior_singularities = unrestricted
                boundary_singularities = S.false
                exterior_singularities = S.false
            else:
                singularities = complexify(singularities)
                unrestricted = complexify(unrestricted)
                interior_singularities = complexify(interior_singularities)
                boundary_singularities = complexify(boundary_singularities)
                exterior_singularities = complexify(exterior_singularities)

    if not return_result:
        return singularities

    restricted_records = []
    if restriction is S.true:
        for records in classified_records.values():
            restricted_records.extend(records)
    else:
        for records in classified_records.values():
            for record in records:
                locus = _restrict_boolean_locus(
                    _singularity_record_locus(record, assumptions=assumptions),
                    restriction,
                    analysis_var,
                    assumptions=assumptions,
                )
                if locus is S.false:
                    continue
                restricted_records.append(
                    Singularity(
                        record.kind,
                        locus,
                        S.true,
                        record.multiplicity,
                        record.source,
                    )
                )
    return SingularityResult(
        singularities=singularities,
        poles=classified_loci["poles"],
        removable=classified_loci["removable"],
        essential=classified_loci["essential"],
        branch_points=classified_loci["branch_points"],
        branch_cuts=classified_loci["branch_cuts"],
        definition_cuts=classified_loci["definition_cuts"],
        discontinuities=classified_loci["discontinuities"],
        records=tuple(restricted_records),
        domain=domain,
        restriction=restriction,
        unrestricted_singularities=unrestricted,
        interior_singularities=interior_singularities,
        boundary_singularities=boundary_singularities,
        exterior_singularities=exterior_singularities,
        domain_interior=domain_interior,
        domain_boundary=domain_boundary,
        domain_exterior=domain_exterior,
    )


def _real_domain_condition(expr, var):
    from ..api import domain

    try:
        d = domain(expr, var, domain=S.Reals)
        if isinstance(d, PropertyResult):
            return None
        return sp.sympify(d)
    except SYMBOLIC_QUERY_ERRORS:
        return None


def _condition_intersects_domain(cond, domain_cond, assumptions=S.true):
    if domain_cond is None:
        return None
    c = simplify_condition(And(cond, domain_cond), assumptions=assumptions)
    if c is S.false:
        return False
    if c is S.true:
        return True
    # For one-dimensional explicit relational conditions SymPy often resolves
    # conjunction consistency via reduce_inequalities.
    syms = tuple(c.free_symbols)
    if len(syms) == 1:
        try:
            red = sp.reduce_inequalities(c, syms[0])
            return red is not S.false
        except SYMBOLIC_QUERY_ERRORS:
            pass
    return None


def _conditional_boolean_status(result):
    """Return True/False when conditional Boolean cases cover all possibilities."""
    if result.primary is S.true:
        return True
    if result.primary is S.false:
        return False
    if not result.cases:
        return None
    if all(c.value is S.true for c in result.cases):
        try:
            if sp.simplify_logic(Or(*[c.condition for c in result.cases])) is S.true:
                return True
        except SYMBOLIC_QUERY_ERRORS:
            pass
    if all(c.value is S.false for c in result.cases):
        try:
            if sp.simplify_logic(Or(*[c.condition for c in result.cases])) is S.true:
                return False
        except SYMBOLIC_QUERY_ERRORS:
            pass
    return None


def _multivariate_structural_analytic(expr, variables, *, complex_mode=False):
    """Conservative joint analyticity classifier for several variables."""
    vv = tuple(variables)
    expr = sp.sympify(expr)
    if not any(expr.has(v) for v in vv):
        return Analyticity.ANALYTIC
    if expr.func in (sp.Abs, sp.sign, sp.floor, sp.ceiling) or isinstance(
        expr, sp.Piecewise
    ):
        return Analyticity.NOT_ANALYTIC
    try:
        if expr.is_polynomial(*vv) or expr.is_rational_function(*vv):
            return Analyticity.ANALYTIC
    except SYMBOLIC_QUERY_ERRORS:
        pass
    if expr.is_Add or expr.is_Mul or (isinstance(expr, sp.Pow) and expr.exp.is_integer):
        vals = [
            _multivariate_structural_analytic(a, vv, complex_mode=complex_mode)
            for a in expr.args
        ]
        if any(v is Analyticity.NOT_ANALYTIC for v in vals):
            return Analyticity.NOT_ANALYTIC
        if all(v is Analyticity.ANALYTIC for v in vals):
            return Analyticity.ANALYTIC
        return Analyticity.UNKNOWN
    analytic_heads = (
        sp.exp,
        sp.sin,
        sp.cos,
        sp.sinh,
        sp.cosh,
        sp.tanh,
        sp.erf,
        sp.erfc,
    )
    if expr.func in analytic_heads:
        vals = [
            _multivariate_structural_analytic(a, vv, complex_mode=complex_mode)
            for a in expr.args
            if any(a.has(v) for v in vv)
        ]
        if all(v is Analyticity.ANALYTIC for v in vals):
            return Analyticity.ANALYTIC
        if any(v is Analyticity.NOT_ANALYTIC for v in vals):
            return Analyticity.NOT_ANALYTIC
    return Analyticity.UNKNOWN


def _multivariate_structural_continuity(expr, variables):
    """Conservative joint continuity classifier on the natural real domain."""
    vv = tuple(variables)
    expr = sp.sympify(expr)
    if not any(expr.has(v) for v in vv):
        return Continuity.CONTINUOUS
    if expr.func in (sp.floor, sp.ceiling, sp.sign):
        return Continuity.DISCONTINUOUS
    try:
        if expr.is_polynomial(*vv) or expr.is_rational_function(*vv):
            return Continuity.CONTINUOUS
    except SYMBOLIC_QUERY_ERRORS:
        pass
    if expr.func is sp.Abs:
        vals = [_multivariate_structural_continuity(a, vv) for a in expr.args]
        return (
            Continuity.CONTINUOUS
            if all(v is Continuity.CONTINUOUS for v in vals)
            else Continuity.UNKNOWN
        )
    if expr.is_Add or expr.is_Mul or isinstance(expr, sp.Pow):
        vals = [
            _multivariate_structural_continuity(a, vv)
            for a in expr.args
            if any(a.has(v) for v in vv)
        ]
        if any(v is Continuity.DISCONTINUOUS for v in vals):
            return Continuity.DISCONTINUOUS
        if all(v is Continuity.CONTINUOUS for v in vals):
            return Continuity.CONTINUOUS
        return Continuity.UNKNOWN
    continuous_heads = (
        sp.exp,
        sp.log,
        sp.sin,
        sp.cos,
        sp.tan,
        sp.cot,
        sp.sinh,
        sp.cosh,
        sp.tanh,
        sp.asin,
        sp.acos,
        sp.atan,
        sp.asinh,
        sp.acosh,
        sp.atanh,
        sp.erf,
        sp.erfc,
        sp.gamma,
    )
    if expr.func in continuous_heads:
        vals = [
            _multivariate_structural_continuity(a, vv)
            for a in expr.args
            if any(a.has(v) for v in vv)
        ]
        if all(v is Continuity.CONTINUOUS for v in vals):
            return Continuity.CONTINUOUS
    return Continuity.UNKNOWN


def _complex_analyticity_stage(expr, var, *, assumptions, return_result):
    # "Analytic" means holomorphic on the function's natural complex
    # domain, matching the source LocallyAnalytic property.
    if expr.func in (sp.Abs, sp.sign, sp.floor, sp.ceiling) or isinstance(
        expr, sp.Piecewise
    ):
        result = PropertyResult(Analyticity.NOT_ANALYTIC, source="nonholomorphic-head")
        return result if return_result else result.primary
    # Generalized hypergeometric functions: DLMF 16.2 states that p<=q
    # is entire in z, while p=q+1 has a principal analytic branch off its
    # standard cut.  A terminating top parameter makes the series a
    # polynomial.  We require bottom parameters to be provably outside
    # the excluded nonpositive integers.
    if expr.func is sp.hyper and len(expr.args) == 3 and expr.args[2] == var:
        aa = tuple(expr.args[0])
        bb = tuple(expr.args[1])

        def _bottom_ok(b):
            return b.is_integer is False or b.is_positive is True

        bottoms_ok = all(_bottom_ok(b) for b in bb)
        terminating = any(a.is_integer is True and a.is_nonpositive is True for a in aa)
        if bottoms_ok and (len(aa) <= len(bb) + 1 or terminating):
            src = (
                "DLMF 16.2(ii)"
                if len(aa) <= len(bb)
                else (
                    "DLMF 16.2 terminating polynomial"
                    if terminating
                    else "DLMF 16.2(iii) principal branch"
                )
            )
            result = PropertyResult(Analyticity.ANALYTIC, source=src)
            return result if return_result else result.primary
    local = evaluate_property(
        expr, "LocallyAnalytic", mode="complex", var=var, assumptions=assumptions
    )
    local_status = _conditional_boolean_status(local)
    if local_status is False:
        return (
            PropertyResult(Analyticity.NOT_ANALYTIC, source="registry")
            if return_result
            else Analyticity.NOT_ANALYTIC
        )
    if local_status is True:
        return (
            PropertyResult(Analyticity.ANALYTIC, source="registry")
            if return_result
            else Analyticity.ANALYTIC
        )
    entire = evaluate_property(
        expr, "Entire", mode="complex", var=var, assumptions=assumptions
    )
    if _conditional_boolean_status(entire) is True or expr.is_polynomial(var):
        return (
            PropertyResult(Analyticity.ANALYTIC, source="analytic-primitive")
            if return_result
            else Analyticity.ANALYTIC
        )
    # Arithmetic closure is valid without needing a primitive function
    # theorem.  Do not apply generic composition closure to an unknown head.
    if expr.is_Add or expr.is_Mul or (isinstance(expr, sp.Pow) and expr.exp.is_integer):
        vals = [
            _function_analytic_impl(a, var, domain=S.Complexes, assumptions=assumptions)
            for a in expr.args
            if a.has(var)
        ]
        if vals and all(v == Analyticity.ANALYTIC for v in vals):
            return (
                PropertyResult(Analyticity.ANALYTIC, source="holomorphic-closure")
                if return_result
                else Analyticity.ANALYTIC
            )
        if any(v == Analyticity.NOT_ANALYTIC for v in vals):
            return (
                PropertyResult(Analyticity.NOT_ANALYTIC, source="holomorphic-closure")
                if return_result
                else Analyticity.NOT_ANALYTIC
            )
    return (
        PropertyResult(Analyticity.UNKNOWN, source="insufficient-evidence")
        if return_result
        else Analyticity.UNKNOWN
    )


def _real_analyticity_stage(
    expr, var, *, domain, assumptions, return_result, domain_ready
):
    domain_cond = (
        sp.sympify(domain)
        if domain_ready
        else (
            _real_domain_condition(expr, var)
            if domain in (S.Reals, "real", "reals", "R", "RR")
            else (
                domain.as_relational(var)
                if isinstance(domain, sp.Set)
                else sp.sympify(domain)
            )
        )
    )
    if not domain_ready and domain_cond is not None:
        split = evaluate_parameter_strata(
            domain_cond,
            var,
            lambda branch_domain: _function_analytic_impl(
                expr,
                var,
                domain=branch_domain,
                assumptions=assumptions,
                return_result=return_result,
                domain_ready=True,
            ),
            source="analyticity-parameter-stratification",
        )
        if split is not None:
            return split
    if expr.func in (sp.Abs, sp.sign, sp.floor, sp.ceiling):
        cuts = definition_cuts(expr, var, mode="real", assumptions=assumptions)
        bad = []
        unknown_cut = False
        for cut in cuts:
            hit = _condition_intersects_domain(
                And(cut.condition, cut.guard), domain_cond, assumptions
            )
            if hit is True:
                bad.append(cut)
            elif hit is None:
                unknown_cut = True
        if bad:
            result = PropertyResult(
                Analyticity.NOT_ANALYTIC, exceptions=tuple(bad), source="nonsmooth-head"
            )
            return result if return_result else result.primary
        if unknown_cut:
            result = PropertyResult(
                Analyticity.UNKNOWN,
                exceptions=cuts,
                source="conditional-nonsmooth-cuts",
            )
            return result if return_result else result.primary
        # Abs/sign/floor/ceiling are locally analytic on every real component
        # that excludes their kink/jump lattice.  This matters for restricted
        # domains such as Abs(x) on x>0.
        result = PropertyResult(Analyticity.ANALYTIC, source="nonsmooth-head-off-cuts")
        return result if return_result else result.primary
    if isinstance(expr, sp.Piecewise):
        boundaries = _piecewise_boundaries(expr)
        bad = []
        unknown_boundary = False
        for b in boundaries:
            hit = _condition_intersects_domain(b, domain_cond, assumptions)
            if hit is True:
                bad.append(b)
            elif hit is None:
                unknown_boundary = True
        if bad:
            exc = tuple(
                Singularity(
                    SingularityKind.DEFINITION_CUT, b, S.true, source="piecewise"
                )
                for b in bad
            )
            result = PropertyResult(
                Analyticity.NOT_ANALYTIC, exceptions=exc, source="piecewise"
            )
            return result if return_result else result.primary
        if unknown_boundary:
            return (
                PropertyResult(Analyticity.UNKNOWN, source="piecewise-boundary-unknown")
                if return_result
                else Analyticity.UNKNOWN
            )
    cuts = definition_cuts(expr, var, mode="real", assumptions=assumptions)
    bad = []
    unknown = False
    for cut in cuts:
        cond = And(cut.condition, cut.guard)
        intersects = _condition_intersects_domain(cond, domain_cond, assumptions)
        if intersects is True:
            bad.append(cut)
        elif intersects is None:
            unknown = True
    if bad:
        result = PropertyResult(
            Analyticity.NOT_ANALYTIC, exceptions=tuple(bad), source="definition-cuts"
        )
        return result if return_result else result.primary
    if unknown:
        return (
            PropertyResult(
                Analyticity.UNKNOWN,
                exceptions=cuts,
                source="conditional-definition-cuts",
            )
            if return_result
            else Analyticity.UNKNOWN
        )

    # Composition theorem: analytic primitives composed with analytic arguments.
    for arg in getattr(expr, "args", ()):
        if arg.has(var) and arg != var:
            a = _function_analytic_impl(
                arg, var, domain=S.Reals, assumptions=assumptions
            )
            if a == Analyticity.NOT_ANALYTIC:
                return (
                    PropertyResult(Analyticity.NOT_ANALYTIC, source="composition")
                    if return_result
                    else Analyticity.NOT_ANALYTIC
                )
            if a == Analyticity.UNKNOWN:
                unknown = True
    entire = evaluate_property(
        expr, "Entire", mode="complex", var=var, assumptions=assumptions
    )
    local = evaluate_property(
        expr, "LocallyAnalytic", mode="real", var=var, assumptions=assumptions
    )
    if expr.is_polynomial(var) or entire.primary is S.true or local.primary is S.true:
        val = Analyticity.UNKNOWN if unknown else Analyticity.ANALYTIC
        return (
            PropertyResult(val, source="registry+composition") if return_result else val
        )
    if expr.is_rational_function(var) and not cuts:
        return (
            PropertyResult(Analyticity.ANALYTIC, source="rational")
            if return_result
            else Analyticity.ANALYTIC
        )
    if expr.is_Add or expr.is_Mul or (isinstance(expr, sp.Pow) and expr.exp.is_integer):
        vals = [
            _function_analytic_impl(a, var, domain=S.Reals, assumptions=assumptions)
            for a in expr.args
            if a.has(var)
        ]
        if vals and all(v == Analyticity.ANALYTIC for v in vals):
            return (
                PropertyResult(Analyticity.ANALYTIC, source="analytic-closure")
                if return_result
                else Analyticity.ANALYTIC
            )
        if any(v == Analyticity.NOT_ANALYTIC for v in vals):
            return (
                PropertyResult(Analyticity.NOT_ANALYTIC, source="analytic-closure")
                if return_result
                else Analyticity.NOT_ANALYTIC
            )
    return (
        PropertyResult(Analyticity.UNKNOWN, source="insufficient-evidence")
        if return_result
        else Analyticity.UNKNOWN
    )


def _function_analytic_impl(
    expr,
    var,
    *,
    domain=S.Reals,
    assumptions=S.true,
    return_result=False,
    _domain_ready=False,
):
    """Classify real analyticity by default, holomorphy for Complexes."""
    expr = sp.sympify(expr)
    vv = _variable_tuple(var)
    is_complex = domain in (S.Complexes, "complex", "complexes", "C", "CC")
    if len(vv) > 1:
        val = _multivariate_structural_analytic(expr, vv, complex_mode=is_complex)
        result = PropertyResult(val, source="multivariate-structural-analyticity")
        return result if return_result else val
    var = vv[0]
    if (
        not is_complex
        and domain not in (S.Reals, "real", "reals", "R", "RR")
        and not isinstance(domain, (sp.Set, Boolean, Relational))
    ):
        raise ValueError(
            "domain must be Reals, Complexes, a set, or a Boolean domain condition"
        )

    # Algebraic expressions are analytic on their natural domain except at
    # branch/definition boundaries, handled below.
    if is_complex:
        return _complex_analyticity_stage(
            expr, var, assumptions=assumptions, return_result=return_result
        )

    return _real_analyticity_stage(
        expr,
        var,
        domain=domain,
        assumptions=assumptions,
        return_result=return_result,
        domain_ready=_domain_ready,
    )


def meromorphic(expr, var, *, assumptions=S.true, return_result=False):
    """Classify meromorphicity on the complex plane/natural complex domain."""
    expr = sp.sympify(expr)
    var = sp.sympify(var)
    if expr.is_rational_function(var):
        return (
            PropertyResult(Meromorphicity.MEROMORPHIC, source="rational")
            if return_result
            else Meromorphicity.MEROMORPHIC
        )
    branch = _function_branchpoint_records(expr, var, assumptions=assumptions)
    cuts = branchcuts(expr, var, assumptions=assumptions, mode="records")
    essential = essential_singularities(expr, var, assumptions=assumptions)
    bad = branch + cuts + essential
    definite = []
    conditional = []
    for item in bad:
        g = simplify_condition(item.guard, assumptions=assumptions)
        if g is S.true:
            definite.append(item)
        elif g is not S.false:
            conditional.append(item)
    if definite:
        result = PropertyResult(
            Meromorphicity.NOT_MEROMORPHIC,
            exceptions=tuple(definite + conditional),
            source="singularity-analysis",
        )
        return result if return_result else result.primary
    if conditional:
        bad_guard = simplify_condition(
            Or(*[item.guard for item in conditional]), assumptions=assumptions
        )
        cases = (
            ConditionalExpression(Meromorphicity.NOT_MEROMORPHIC, bad_guard),
            ConditionalExpression(Meromorphicity.MEROMORPHIC, sp.Not(bad_guard)),
        )
        result = PropertyResult(
            Meromorphicity.UNKNOWN,
            cases=cases,
            exceptions=tuple(conditional),
            source="conditional-singularity-analysis",
        )
        return result
    reg = evaluate_property(
        expr, "Meromorphic", mode="complex", assumptions=assumptions
    )
    if reg.primary is S.true:
        return (
            PropertyResult(Meromorphicity.MEROMORPHIC, source="registry")
            if return_result
            else Meromorphicity.MEROMORPHIC
        )
    if reg.primary is S.false:
        return (
            PropertyResult(Meromorphicity.NOT_MEROMORPHIC, source="registry")
            if return_result
            else Meromorphicity.NOT_MEROMORPHIC
        )
    entire = evaluate_property(expr, "Entire", mode="complex", assumptions=assumptions)
    if entire.primary is S.true:
        return (
            PropertyResult(Meromorphicity.MEROMORPHIC, source="entire")
            if return_result
            else Meromorphicity.MEROMORPHIC
        )
    if (
        _function_analytic_impl(expr, var, domain=S.Complexes, assumptions=assumptions)
        == Analyticity.ANALYTIC
    ):
        # Holomorphic on its natural domain and with no branch/essential
        # obstruction above: any remaining isolated singularities are poles.
        return (
            PropertyResult(Meromorphicity.MEROMORPHIC, source="analytic-with-poles")
            if return_result
            else Meromorphicity.MEROMORPHIC
        )
    # Closure of meromorphic functions under arithmetic.
    if expr.is_Add or expr.is_Mul or (isinstance(expr, sp.Pow) and expr.exp.is_integer):
        vals = []
        for arg in expr.args:
            if arg.has(var):
                vals.append(meromorphic(arg, var, assumptions=assumptions))
        if vals and all(v == Meromorphicity.MEROMORPHIC for v in vals):
            return (
                PropertyResult(Meromorphicity.MEROMORPHIC, source="closure")
                if return_result
                else Meromorphicity.MEROMORPHIC
            )
        if any(v == Meromorphicity.NOT_MEROMORPHIC for v in vals):
            return (
                PropertyResult(Meromorphicity.NOT_MEROMORPHIC, source="closure")
                if return_result
                else Meromorphicity.NOT_MEROMORPHIC
            )
    return (
        PropertyResult(Meromorphicity.UNKNOWN, source="insufficient-evidence")
        if return_result
        else Meromorphicity.UNKNOWN
    )


def _condition_to_set(condition, var):
    """Convert a real condition using the shared structural set converter."""
    if condition is None:
        return None

    # Import lazily to avoid an import cycle through real_analysis.smoothness.
    from ..real_analysis.common import condition_to_set

    converted = condition_to_set(condition, var)
    if converted is not None:
        return converted
    try:
        return sp.solveset(condition, var, domain=S.Reals)
    except SYMBOLIC_QUERY_ERRORS:
        return None


def _function_continuous_impl(
    expr,
    var,
    *,
    domain=S.Reals,
    assumptions=S.true,
    return_result=False,
    _domain_ready=False,
):
    """Test continuity on the function's natural real domain by default."""
    expr = sp.sympify(expr)
    vv = _variable_tuple(var)
    if len(vv) > 1:
        val = _multivariate_structural_continuity(expr, vv)
        result = PropertyResult(val, source="multivariate-structural-continuity")
        return result if return_result else val
    var = vv[0]
    domain_cond = (
        sp.sympify(domain)
        if _domain_ready
        else (
            _real_domain_condition(expr, var)
            if domain in (S.Reals, "real", "reals", "R", "RR")
            else (
                domain.as_relational(var)
                if isinstance(domain, sp.Set)
                else sp.sympify(domain)
            )
        )
    )
    if not _domain_ready and domain_cond is not None:
        split = evaluate_parameter_strata(
            domain_cond,
            var,
            lambda branch_domain: _function_continuous_impl(
                expr,
                var,
                domain=branch_domain,
                assumptions=assumptions,
                return_result=return_result,
                _domain_ready=True,
            ),
            source="continuity-parameter-stratification",
        )
        if split is not None:
            return split
    discontinuities = _function_discontinuity_records(
        expr, var, assumptions=assumptions
    )
    actual = []
    uncertain = []
    for d in discontinuities:
        hit = _condition_intersects_domain(
            And(d.condition, d.guard), domain_cond, assumptions
        )
        if hit is True:
            actual.append(d)
        elif hit is None:
            uncertain.append(d)
    if actual:
        result = PropertyResult(
            Continuity.DISCONTINUOUS,
            exceptions=tuple(actual),
            source="discontinuity-analysis",
        )
        return result if return_result else result.primary
    if isinstance(expr, sp.Piecewise) and not uncertain:
        from ..real_analysis.piecewise import piecewise_branches

        branches = piecewise_branches(expr, var, domain_cond)
        branch_values = [
            _function_continuous_impl(
                branch.expression,
                var,
                domain=branch.condition,
                assumptions=assumptions,
                _domain_ready=True,
            )
            for branch in branches
            if branch.expression.has(var)
        ]
        if all(v == Continuity.CONTINUOUS for v in branch_values):
            return (
                PropertyResult(Continuity.CONTINUOUS, source="piecewise-analysis")
                if return_result
                else Continuity.CONTINUOUS
            )
    natural_set = _condition_to_set(domain_cond, var)
    if natural_set is not None:
        try:
            cset = continuous_domain(expr, var, S.Reals)
            if cset == natural_set:
                return (
                    PropertyResult(
                        Continuity.CONTINUOUS, source="sympy.continuous_domain"
                    )
                    if return_result
                    else Continuity.CONTINUOUS
                )
            missing = natural_set - cset
            if missing is not S.EmptySet and missing != S.EmptySet:
                result = PropertyResult(
                    Continuity.DISCONTINUOUS,
                    exceptions=(missing,),
                    source="sympy.continuous_domain",
                )
                return result if return_result else result.primary
        except NotImplementedError:
            pass
        except SYMBOLIC_QUERY_ERRORS:
            pass

    # Analytic functions are continuous on their real analytic domain.
    analytic_value = _function_analytic_impl(
        expr, var, domain=domain, assumptions=assumptions
    )
    if analytic_value == Analyticity.ANALYTIC:
        return (
            PropertyResult(Continuity.CONTINUOUS, source="analyticity")
            if return_result
            else Continuity.CONTINUOUS
        )
    if expr.func is sp.Abs:
        return (
            PropertyResult(Continuity.CONTINUOUS, source="Abs")
            if return_result
            else Continuity.CONTINUOUS
        )
    return (
        PropertyResult(Continuity.UNKNOWN, source="insufficient-evidence")
        if return_result
        else Continuity.UNKNOWN
    )
