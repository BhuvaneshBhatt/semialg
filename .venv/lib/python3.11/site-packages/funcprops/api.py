"""Public routing layer for function-property queries.

The wrappers in this module also normalize ``ConditionalExpression`` inputs.  The
underlying property engines remain ordinary SymPy-expression analyzers.
"""

from __future__ import annotations

import functools
import inspect

import sympy as sp
from sympy import S

from .assumptions import normalize_assumptions
from .conditional import (
    combine_assumptions,
    normalize_conditional_input,
    restrict_domain,
)
from .geometry import (
    boundary_behavior as _function_boundary_behavior,
)
from .geometry import (
    geometry as _function_geometry,
)
from .geometry import (
    level_boundedness as _function_level_boundedness,
)
from .geometry import (
    mode_existence as _function_mode_existence,
)
from .real_analysis.boundedness import boundedness as _function_boundedness
from .real_analysis.complete_monotonicity import (
    complete_monotonicity as _function_complete_monotonicity,
)
from .real_analysis.convexity import convexity as _function_convexity
from .real_analysis.differentiability import (
    differentiability as _function_differentiability,
)
from .real_analysis.domain import domain as _function_domain
from .real_analysis.integrability import integrability as _function_integrability
from .real_analysis.mapping import (
    bijective as _function_bijective,
)
from .real_analysis.mapping import (
    injective as _function_injective,
)
from .real_analysis.mapping import (
    surjective as _function_surjective,
)
from .real_analysis.monotonicity import monotonicity as _function_monotonicity
from .real_analysis.oscillation import oscillation as _function_oscillation
from .real_analysis.range import range as _function_range
from .real_analysis.regularity_extensions import (
    holder_regularity as _function_holder_regularity,
)
from .real_analysis.regularity_extensions import (
    lipschitz as _function_lipschitz,
)
from .real_analysis.regularity_extensions import (
    uniformly_continuous as _function_uniformly_continuous,
)
from .real_analysis.sign import sign as _function_sign
from .real_analysis.smoothness import (
    smooth as _smooth,
)
from .real_analysis.smoothness import (
    smoothness as _function_smoothness,
)
from .real_analysis.structural_properties import (
    homogeneity as _function_homogeneity,
)
from .real_analysis.structural_properties import (
    periodicity as _function_periodicity,
)
from .real_analysis.structural_properties import (
    symmetry as _function_symmetry,
)
from .real_analysis.variation import variation as _function_variation
from .real_analysis.zeros import zeros as _function_zeros
from .results import PropertyResult, as_detailed_result
from .singularities import (
    analytic as _function_analytic,
)
from .singularities import (
    branchcuts as _function_branchcuts,
)
from .singularities import (
    branchpoints as _function_branchpoints,
)
from .singularities import (
    continuous as _function_continuous,
)
from .singularities import (
    discontinuities as _function_discontinuities,
)
from .singularities import (
    meromorphic as _function_meromorphic,
)
from .singularities import (
    poles as _function_poles,
)
from .singularities import (
    singularities as _function_singularities,
)
from .singularities import (
    singularity_analysis as _function_singularity_analysis,
)
from .stratification import resolve_result_assumptions
from .transcendental.complex_roots import (
    complex_zeros as _function_complex_zeros,
)
from .transformation import (
    change_of_variables as _function_change_of_variables,
)
from .transformation import (
    determinant_nonvanishing as _function_determinant_nonvanishing,
)
from .transformation import (
    jacobian as _function_jacobian,
)
from .transformation import (
    local_invertibility as _function_local_invertibility,
)
from .transformation import (
    transformation_geometry as _function_transformation_geometry,
)


def _public_signature(sig):
    """Return the uniform root-API signature for an internal property engine."""
    params = []
    for param in sig.parameters.values():
        if param.name.startswith("_"):
            continue
        name = {"var": "variable", "vars": "variables"}.get(param.name, param.name)
        if name == "return_result":
            continue
        params.append(param.replace(name=name))

    insertion = len(params)
    for i, param in enumerate(params):
        if param.kind is inspect.Parameter.VAR_KEYWORD:
            insertion = i
            break
    names = {param.name for param in params}
    additions = []
    if "assumptions" not in names:
        additions.append(
            inspect.Parameter(
                "assumptions", inspect.Parameter.KEYWORD_ONLY, default=True
            )
        )
    additions.append(
        inspect.Parameter(
            "return_result", inspect.Parameter.KEYWORD_ONLY, default=False
        )
    )
    for param in reversed(additions):
        params.insert(insertion, param)
    return sig.replace(parameters=params)


def _conditional_wrapper(fn, *, domain_mode="ordinary"):
    sig = inspect.signature(fn)
    has_native_assumptions = "assumptions" in sig.parameters
    has_native_return_result = "return_result" in sig.parameters

    @functools.wraps(fn)
    def wrapped(expr, *args, **kwargs):
        private = [name for name in kwargs if name.startswith("_")]
        if private:
            raise TypeError(
                f"internal keyword is not part of the public API: {private[0]}"
            )
        if "variable" in kwargs and "var" in sig.parameters:
            kwargs["var"] = kwargs.pop("variable")
        if "variables" in kwargs and "vars" in sig.parameters:
            kwargs["vars"] = kwargs.pop("variables")
        return_result = bool(kwargs.pop("return_result", False))
        normalized = normalize_conditional_input(expr)
        public_assumptions = normalize_assumptions(kwargs.pop("assumptions", True))
        combined = normalize_assumptions(
            combine_assumptions(public_assumptions, normalized.condition)
        )

        bound = sig.bind_partial(normalized.expr, *args, **kwargs)
        var = bound.arguments.get(
            "var", bound.arguments.get("variable", bound.arguments.get("variables"))
        )
        if var is None and args:
            var = args[0]
        if var is not None:
            try:
                var = sp.sympify(var)
            except (TypeError, ValueError) as exc:
                raise TypeError(
                    "analysis variable(s) must be SymPy-compatible symbols or tuples"
                ) from exc
            scalar_parameter = "var" in sig.parameters or "variable" in sig.parameters
            tuple_capable = fn.__name__ in {
                "injective",
                "surjective",
                "bijective",
                "convexity",
                "analytic",
                "continuous",
            }
            if scalar_parameter and not tuple_capable and isinstance(var, sp.Tuple):
                raise TypeError(
                    "expected one analysis variable, but received multiple variables; "
                    "use a multivariate property such as jacobian(), geometry(), or "
                    "homogeneity() when analyzing several variables"
                )
            if (
                scalar_parameter
                and not tuple_capable
                and not isinstance(var, sp.Symbol)
            ):
                raise TypeError(
                    f"analysis variable must be a Symbol, not {type(var).__name__}"
                )

        if has_native_assumptions:
            kwargs["assumptions"] = combined
        elif "domain" in sig.parameters and domain_mode == "ordinary":
            domain = bound.arguments.get("domain", sig.parameters["domain"].default)
            if normalized.condition is not S.true:
                kwargs["domain"] = restrict_domain(domain, normalized.condition, var)

        if (
            "domain" in sig.parameters
            and domain_mode == "ordinary"
            and has_native_assumptions
        ):
            domain = bound.arguments.get("domain", sig.parameters["domain"].default)
            if normalized.condition is not S.true:
                kwargs["domain"] = restrict_domain(domain, normalized.condition, var)

        if has_native_return_result:
            kwargs["return_result"] = return_result
        result = fn(normalized.expr, *args, **kwargs)
        result = resolve_result_assumptions(result, combined)
        if isinstance(result, PropertyResult):
            if return_result:
                return result
            result = result.as_piecewise()
        if domain_mode == "natural-domain" and normalized.condition is not S.true:
            restriction = normalized.condition
            if isinstance(result, sp.Set) and not isinstance(var, sp.Tuple):
                from .real_analysis.common import condition_to_set

                restricted_set = condition_to_set(restriction, var)
                if restricted_set is not None:
                    result = sp.Intersection(result, restricted_set)
                else:
                    result = sp.And(sp.Contains(var, result), restriction)
            elif isinstance(result, sp.Basic):
                result = sp.And(result, restriction)
        if return_result and not has_native_return_result:
            domain = bound.arguments.get("domain")
            return as_detailed_result(result, domain=domain, source=fn.__name__)
        return result

    wrapped.__signature__ = _public_signature(sig)
    return wrapped


boundedness = _conditional_wrapper(_function_boundedness)
convexity = _conditional_wrapper(_function_convexity)
complete_monotonicity = _conditional_wrapper(_function_complete_monotonicity)
differentiability = _conditional_wrapper(_function_differentiability)
domain = _conditional_wrapper(_function_domain, domain_mode="natural-domain")
bijective = _conditional_wrapper(_function_bijective)
injective = _conditional_wrapper(_function_injective)
surjective = _conditional_wrapper(_function_surjective)
monotonicity = _conditional_wrapper(_function_monotonicity)
oscillation = _conditional_wrapper(_function_oscillation)
range = _conditional_wrapper(_function_range)
sign = _conditional_wrapper(_function_sign)
smoothness = _conditional_wrapper(_function_smoothness)
smooth = _conditional_wrapper(_smooth)
lipschitz = _conditional_wrapper(_function_lipschitz)
uniformly_continuous = _conditional_wrapper(_function_uniformly_continuous)
holder_regularity = _conditional_wrapper(_function_holder_regularity)
variation = _conditional_wrapper(_function_variation)
integrability = _conditional_wrapper(_function_integrability)
homogeneity = _conditional_wrapper(_function_homogeneity)
periodicity = _conditional_wrapper(_function_periodicity)
symmetry = _conditional_wrapper(_function_symmetry)
zeros = _conditional_wrapper(_function_zeros)
analytic = _conditional_wrapper(_function_analytic)


def _restrict_condition_locus(locus, condition, var):
    combined = sp.And(sp.sympify(locus), sp.sympify(condition))
    # Branch-cut conditions commonly carry Contains(var, Reals).  Under that
    # guard, refine re(var)/im(var) before reducing inequalities.
    if var is not None and combined.has(sp.Contains(var, S.Reals)):
        combined = sp.refine(combined, sp.Q.real(var))
        combined = combined.xreplace({sp.re(var): var, sp.im(var): S.Zero})
    try:
        reduced = (
            sp.reduce_inequalities(combined, var)
            if var is not None
            else sp.simplify(combined)
        )
    except (NotImplementedError, TypeError, ValueError):
        reduced = sp.simplify(combined)
    return reduced


def branch_cuts(
    expr, variable, *, assumptions=True, mode="conditions", return_result=False
):
    var = variable
    normalized = normalize_conditional_input(expr)
    assumptions = normalize_assumptions(
        combine_assumptions(assumptions, normalized.condition)
    )
    raw = _function_branchcuts(normalized.expr, var, assumptions=assumptions, mode=mode)
    if normalized.condition is S.true or mode == "jumps":
        return as_detailed_result(raw, source="branch_cuts") if return_result else raw
    if mode == "conditions":
        out = []
        for locus in raw:
            restricted = _restrict_condition_locus(
                locus, normalized.condition, sp.sympify(var)
            )
            if restricted is not S.false and restricted != S.false:
                out.append(restricted)
        result = tuple(out)
        return (
            as_detailed_result(result, source="branch_cuts")
            if return_result
            else result
        )
    if mode == "union":
        result = _restrict_condition_locus(raw, normalized.condition, sp.sympify(var))
        return (
            as_detailed_result(result, source="branch_cuts")
            if return_result
            else result
        )
    if mode == "records":
        from .results import Singularity

        out = []
        for record in raw:
            locus = _restrict_condition_locus(
                record.condition, normalized.condition, sp.sympify(var)
            )
            if locus is S.false or locus == S.false:
                continue
            out.append(
                Singularity(
                    record.kind, locus, record.guard, record.multiplicity, record.source
                )
            )
        result = tuple(out)
        return (
            as_detailed_result(result, source="branch_cuts")
            if return_result
            else result
        )
    return as_detailed_result(raw, source="branch_cuts") if return_result else raw


branch_points = _conditional_wrapper(_function_branchpoints)
continuous = _conditional_wrapper(_function_continuous)
discontinuities = _conditional_wrapper(_function_discontinuities)
meromorphic = _conditional_wrapper(_function_meromorphic)
poles = _conditional_wrapper(_function_poles)


def singularities(
    expr,
    variable,
    domain=S.Complexes,
    *,
    assumptions=True,
    return_result=False,
):
    """Return the Boolean singularity locus, preserving function restrictions.

    ``Dummy`` analysis variables are normalized through an ordinary proxy
    symbol.  Several SymPy set/solver paths introduce their own dummies; using
    a proxy prevents those implementation symbols from leaking into the public
    singularity locus and makes the query invariant to ``Symbol`` versus
    ``Dummy`` callers.
    """
    assumptions = normalize_assumptions(assumptions)
    original_var = sp.sympify(variable)
    normalized = normalize_conditional_input(expr)
    concrete_expr = sp.sympify(normalized.expr)
    restriction = sp.sympify(normalized.condition)
    analysis_var = original_var
    if isinstance(original_var, sp.Dummy):
        analysis_var = sp.Symbol("_funcprops_var")
        replacement = {original_var: analysis_var}
        concrete_expr = concrete_expr.xreplace(replacement)
        restriction = restriction.xreplace(replacement)
        assumptions = sp.sympify(assumptions).xreplace(replacement)
    if restriction is not S.true:
        from .results import ConditionalExpression

        concrete_expr = ConditionalExpression(concrete_expr, restriction)
    result = _function_singularities(
        concrete_expr,
        analysis_var,
        domain,
        assumptions=assumptions,
        return_result=return_result,
    )
    if analysis_var != original_var and isinstance(result, sp.Basic):
        result = result.xreplace({analysis_var: original_var})
    return result


singularity_analysis = _conditional_wrapper(_function_singularity_analysis)
boundary_behavior = _conditional_wrapper(_function_boundary_behavior)
level_boundedness = _conditional_wrapper(_function_level_boundedness)
mode_existence = _conditional_wrapper(_function_mode_existence)
geometry = _conditional_wrapper(_function_geometry)


def jacobian(expressions, variables, *, assumptions=True, return_result=False):
    """Return a transformation Jacobian matrix or its detailed certificate."""
    return _function_jacobian(expressions, variables, return_result=return_result)


def determinant_nonvanishing(
    expressions, variables, *, domain=S.Reals, assumptions=True, return_result=False
):
    """Return whether the Jacobian determinant is nonzero throughout the domain."""
    return _function_determinant_nonvanishing(
        expressions,
        variables,
        domain=domain,
        assumptions=assumptions,
        return_result=return_result,
    )


def local_invertibility(
    expressions, variables, *, domain=S.Reals, assumptions=True, return_result=False
):
    """Return whether local invertibility is certified throughout the domain."""
    return _function_local_invertibility(
        expressions,
        variables,
        domain=domain,
        assumptions=assumptions,
        return_result=return_result,
    )


def change_of_variables(
    expressions, variables, *, domain=S.Reals, assumptions=True, return_result=False
):
    return _function_change_of_variables(
        expressions, variables, domain=domain, assumptions=assumptions
    )


def transformation_geometry(
    expressions, variables, *, domain=S.Reals, assumptions=True, return_result=False
):
    return _function_transformation_geometry(
        expressions, variables, domain=domain, assumptions=assumptions
    )


# Complex rectangle isolation cannot in general intersect an arbitrary Boolean
# condition with a rectangle without a certified planar-set decomposition.  We
# nevertheless accept unconditional ConditionalExpression wrappers and reject only
# genuinely restricted inputs explicitly.
def complex_zeros(expr, *args, assumptions=True, return_result=False, **kwargs):
    normalized = normalize_conditional_input(expr)
    assumptions = normalize_assumptions(
        combine_assumptions(assumptions, normalized.condition)
    )
    if assumptions is not S.true:
        raise NotImplementedError(
            "assumption-restricted complex rectangle root isolation requires "
            "certified planar-domain transport"
        )
    return _function_complex_zeros(normalized.expr, *args, **kwargs)


_PROPERTY_CATALOG = {
    "domain": ("sets", "Natural or restricted domain."),
    "range": ("sets", "Exact real range condition."),
    "zeros": ("zeros", "Real zero set or certified zero geometry."),
    "complex_zeros": ("zeros", "Certified zeros in a complex rectangle."),
    "sign": ("real analysis", "Sign classification on a domain."),
    "boundedness": ("real analysis", "Boundedness classification."),
    "boundary_behavior": (
        "real analysis",
        "Certified boundary limits and escape behavior.",
    ),
    "monotonicity": ("real analysis", "Monotonicity classification."),
    "complete_monotonicity": ("real analysis", "Complete-monotonicity classification."),
    "convexity": ("real analysis", "Convexity and concavity classification."),
    "oscillation": ("real analysis", "Oscillation classification."),
    "continuous": ("regularity", "Continuity classification."),
    "differentiability": ("regularity", "Differentiability classification."),
    "smoothness": ("regularity", "Smoothness classification."),
    "lipschitz": ("regularity", "Lipschitz regularity and certified constant."),
    "uniformly_continuous": ("regularity", "Uniform-continuity classification."),
    "holder_regularity": ("regularity", "Certified Hölder exponent."),
    "variation": ("real analysis", "Bounded variation and total variation."),
    "integrability": ("real analysis", "Structured integrability classification."),
    "analytic": ("regularity", "Real-analyticity classification."),
    "meromorphic": ("complex analysis", "Meromorphicity classification."),
    "periodicity": ("structure", "Periodicity and fundamental-period information."),
    "symmetry": ("structure", "Even, odd, central, or other symmetry."),
    "homogeneity": ("structure", "Homogeneity classification and degree."),
    "injective": ("mapping", "Injectivity classification."),
    "surjective": ("mapping", "Surjectivity classification."),
    "bijective": ("mapping", "Bijectivity classification."),
    "singularities": ("complex analysis", "Singularity locus."),
    "poles": ("complex analysis", "Poles with multiplicities."),
    "branch_points": ("complex analysis", "Branch-point locus."),
    "branch_cuts": ("complex analysis", "Principal branch-cut geometry."),
    "discontinuities": ("regularity", "Discontinuity locus."),
    "singularity_analysis": (
        "complex analysis",
        "Detailed singularity classification.",
    ),
    "jacobian": ("transformations", "Exact Jacobian matrix."),
    "determinant_nonvanishing": (
        "transformations",
        "Domain-wide Jacobian nonvanishing.",
    ),
    "local_invertibility": ("transformations", "Local invertibility classification."),
    "change_of_variables": ("transformations", "Change-of-variables geometry."),
    "transformation_geometry": ("transformations", "Combined transformation geometry."),
    "level_boundedness": ("optimization", "Level-boundedness classification."),
    "mode_existence": (
        "optimization",
        "Structural extremum-attainment classification.",
    ),
    "geometry": ("optimization", "Combined function geometry."),
}


def properties(*, category=None):
    """List discoverable root-level mathematical properties.

    Pass ``category=...`` to filter the catalog.  The return value is a tuple of
    ``(name, category, summary)`` triples so it is useful in terminals, notebooks,
    and documentation tooling without introducing another public result type.
    """
    requested = category.casefold() if category is not None else None
    rows = []
    for name, (group, summary) in sorted(_PROPERTY_CATALOG.items()):
        if requested is None or group.casefold() == requested:
            rows.append((name, group, summary))
    if requested is not None and not rows:
        choices = ", ".join(sorted({group for group, _ in _PROPERTY_CATALOG.values()}))
        raise ValueError(
            f"unknown property category {category!r}; choose from: {choices}"
        )
    return tuple(rows)


__all__ = [
    "analytic",
    "bijective",
    "boundary_behavior",
    "boundedness",
    "branch_cuts",
    "branch_points",
    "change_of_variables",
    "complete_monotonicity",
    "complex_zeros",
    "continuous",
    "convexity",
    "determinant_nonvanishing",
    "differentiability",
    "discontinuities",
    "domain",
    "geometry",
    "homogeneity",
    "injective",
    "jacobian",
    "level_boundedness",
    "local_invertibility",
    "meromorphic",
    "mode_existence",
    "monotonicity",
    "oscillation",
    "periodicity",
    "poles",
    "range",
    "sign",
    "singularities",
    "singularity_analysis",
    "smoothness",
    "surjective",
    "symmetry",
    "transformation_geometry",
    "zeros",
]
