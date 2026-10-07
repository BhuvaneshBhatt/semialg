"""Adapter to semialg's public exact real-function API.

funcprops owns the public property vocabulary and evidence objects. semialg
provides the exact backend for semialgebraic functions on a semialgebraic real
domains, including domain, range, sign, monotonicity, convexity, smoothness,
and mapping-property calculations. Imports remain lazy so catalogue and
complex-analysis queries do not load the backend until a real-domain proof
needs it.
"""

from __future__ import annotations

from collections.abc import Callable
from functools import cache

import sympy as sp


class SemialgUnavailable(ImportError):
    """Raised when the required semialg public API cannot be loaded."""


@cache
def _public_function(name: str) -> Callable:
    """Return one function from semialg's documented root API."""
    try:
        import semialg
    except ImportError as exc:
        raise SemialgUnavailable(
            "funcprops real-analysis queries require semialg>=1.3.0"
        ) from exc

    try:
        function = getattr(semialg, name)
    except AttributeError as exc:
        raise SemialgUnavailable(
            f"semialg>=1.3.0 does not provide the required public function {name!r}"
        ) from exc
    if not callable(function):
        raise SemialgUnavailable(f"semialg public attribute {name!r} is not callable")
    return function


def domain(expr, variable):
    """Return semialg's recognized real-domain condition."""
    function = _public_function("function_domain")
    return function(sp.sympify(expr), [variable])


def range(expr, variable, value_symbol, *, constraints=True, method="qe", **kwargs):
    """Return an exact real-range condition from semialg."""
    function = _public_function("function_range")
    return function(
        sp.sympify(expr),
        constraints,
        [variable],
        value_symbol=value_symbol,
        method=method,
        **kwargs,
    )


def convexity(expr, variable, *, domain=True, strict_inequalities=True, **kwargs):
    """Return semialg's strongest certified convexity classification."""
    properties = kwargs.pop("properties", "primary")
    return_result = kwargs.pop("return_result", False)
    parameters = kwargs.pop("parameters", None)
    if kwargs:
        names = ", ".join(sorted(kwargs))
        raise TypeError(f"unsupported semialg convexity options: {names}")

    # semialg computes the strongest exact class. funcprops weakens the class
    # afterwards when callers request non-strict semantics.
    del strict_inequalities
    function = _public_function("function_convexity")
    return function(
        sp.sympify(expr),
        [variable],
        domain=domain,
        parameters=parameters,
        return_result=return_result,
        properties=properties,
    )


def monotonicity(expr, variable, *, domain=True, return_result=False):
    """Return semialg's exact monotonicity classification."""
    function = _public_function("function_monotonicity")
    return function(
        sp.sympify(expr),
        variable,
        domain=domain,
        return_result=return_result,
    )


def smoothness(expr, variable, *, domain=True, max_order=6):
    """Return semialg's exact regularity result on a semialgebraic domain."""
    function = _public_function("function_smoothness")
    return function(
        sp.sympify(expr),
        variable,
        domain=domain,
        max_order=max_order,
    )


def function_mapping_properties(
    expr, variable, *, domain=True, codomain=True, image_variable=None
):
    """Return semialg's exact mapping properties for a semialgebraic map."""
    function = _public_function("function_mapping_properties")
    image_variables = None if image_variable is None else (image_variable,)
    return function(
        expr,
        variable,
        domain=domain,
        codomain=codomain,
        image_variables=image_variables,
    )


def eligible_problem(expr, variable, *, domain=True) -> bool:
    """Return whether an expression/domain pair is safe for semialg delegation."""
    from .real_analysis.common import (
        is_semialgebraic_condition,
        is_semialgebraic_expression,
    )

    return is_semialgebraic_expression(expr, variable) and is_semialgebraic_condition(
        domain, variable
    )


def is_satisfiable(condition, variables):
    """Decide satisfiability of a semialgebraic parameter guard exactly."""
    function = _public_function("is_satisfiable")
    normalized = tuple(sp.sympify(v) for v in variables)
    return bool(function(sp.sympify(condition), normalized))


def sign(expr, variable, *, assumptions=True, return_result=False):
    """Return semialg's canonical exact sign classification."""
    function = _public_function("function_sign")
    return function(
        sp.sympify(expr),
        [variable],
        assumptions=assumptions,
        return_result=return_result,
    )


def sign_partition(expr, variable, *, domain=True):
    """Return semialg's exact univariate sign partition."""
    function = _public_function("function_sign_partition")
    return function(sp.sympify(expr), variable, domain=domain)


def normalize_sign(value):
    """Translate a semialg sign classification into funcprops' Sign enum."""
    from .results import Sign

    classification = getattr(value, "classification", value)
    return {
        "positive": Sign.POSITIVE,
        "negative": Sign.NEGATIVE,
        "zero": Sign.ZERO,
        "nonnegative": Sign.NONNEGATIVE,
        "nonpositive": Sign.NONPOSITIVE,
        "mixed": Sign.SIGN_CHANGING,
    }.get(str(classification), Sign.UNKNOWN)


def normalize_convexity(value):
    """Translate semialg convexity vocabulary into funcprops' enum."""
    from .results import Convexity

    if isinstance(value, Convexity):
        return value
    if hasattr(value, "primary"):
        value = value.primary
    if isinstance(value, int):
        return {
            1: Convexity.STRICTLY_CONVEX,
            2: Convexity.CONVEX,
            -1: Convexity.STRICTLY_CONCAVE,
            -2: Convexity.CONCAVE,
        }.get(int(value), Convexity.UNKNOWN)
    key = str(getattr(value, "name", value)).lower().replace("_", "").replace(" ", "")
    return {
        "stronglyconvex": Convexity.STRONGLY_CONVEX,
        "strictlyconvex": Convexity.STRICTLY_CONVEX,
        "convex": Convexity.CONVEX,
        "affine": Convexity.AFFINE,
        "strictlyconcave": Convexity.STRICTLY_CONCAVE,
        "concave": Convexity.CONCAVE,
        "stronglyconcave": Convexity.STRONGLY_CONCAVE,
        "neither": Convexity.NEITHER,
    }.get(key, Convexity.UNKNOWN)


def normalize_monotonicity(value):
    """Translate semialg monotonicity vocabulary into funcprops' enum."""
    from .results import Monotonicity

    classification = getattr(value, "classification", value)
    return {
        "strictly_increasing": Monotonicity.STRICTLY_INCREASING,
        "increasing": Monotonicity.NONDECREASING,
        "constant": Monotonicity.CONSTANT,
        "decreasing": Monotonicity.NONINCREASING,
        "strictly_decreasing": Monotonicity.STRICTLY_DECREASING,
        "nonmonotonic": Monotonicity.NOT_MONOTONE,
    }.get(str(classification), Monotonicity.UNKNOWN)


def normalize_truth(value):
    """Translate a backend boolean into funcprops' three-valued truth type."""
    from .results import TruthValue

    if value is True:
        return TruthValue.TRUE
    if value is False:
        return TruthValue.FALSE
    return TruthValue.UNKNOWN


def prove_sign(expr, variable, relation, *, assumptions=True):
    """Prove one exact sign relation with semialg's public sign provers."""
    prover_names = {
        "positive": "prove_positive",
        "negative": "prove_negative",
        "nonnegative": "prove_nonnegative",
        "nonpositive": "prove_nonpositive",
    }
    try:
        prover_name = prover_names[relation]
    except KeyError as exc:
        raise ValueError(f"unknown sign relation: {relation!r}") from exc

    prover = _public_function(prover_name)
    return prover(sp.sympify(expr), [variable], assumptions=assumptions)
