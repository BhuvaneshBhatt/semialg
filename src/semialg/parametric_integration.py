from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass

import sympy as sp

from .standard_regions import ParametricRegion


@dataclass(frozen=True)
class ParametricIntegralResult:
    value: sp.Expr
    integrand: sp.Expr
    ambient_variables: tuple[sp.Symbol, ...]
    region: ParametricRegion
    jacobian_factor: sp.Expr
    exact: bool = True
    method: str = "parametric_metric_jacobian"


def metric_jacobian_factor(mapping: Sequence[sp.Expr], parameters: Sequence[sp.Symbol]) -> sp.Expr:
    """Return the metric Jacobian for a parametrized region.

    For full-dimensional square parametrizations this is ``Abs(det(J))``. For
    lower-dimensional parametrizations it is ``sqrt(det(J.T*J))``.
    """

    mapping_vec = sp.Matrix(mapping)
    params = tuple(parameters)
    if not params:
        return sp.Integer(1)
    jac = mapping_vec.jacobian(params)
    if jac.rows == jac.cols:
        return sp.Abs(sp.det(jac))
    gram = jac.T * jac
    return sp.sqrt(sp.simplify(sp.det(gram)))


def reduce_parametric_integral(
    integrand: object,
    ambient_variables: Sequence[sp.Symbol | str],
    region: ParametricRegion,
) -> tuple[sp.Expr, tuple[tuple[sp.Symbol, sp.Expr, sp.Expr], ...], sp.Expr]:
    """Reduce an integral over a parametrized region to an iterated integral.

    Returns ``(transformed_integrand, limits, jacobian_factor)``.
    """

    vars_ = tuple(sp.Symbol(v, real=True) if isinstance(v, str) else v for v in ambient_variables)
    if len(vars_) != len(region.mapping):
        raise ValueError("ambient variable count must match the parametrization mapping dimension")
    expr = sp.sympify(integrand)
    jac_factor = metric_jacobian_factor(region.mapping, region.parameters)
    transformed = sp.simplify(
        expr.subs(dict(zip(vars_, region.mapping, strict=True))) * jac_factor / region.multiplicity
    )
    return transformed, region.limits, jac_factor


def integrate_over_parametric_region(
    integrand: object,
    ambient_variables: Sequence[sp.Symbol | str],
    region: ParametricRegion,
    *,
    method: str = "symbolic",
    precision: int = 50,
    return_result: bool = False,
) -> sp.Expr | ParametricIntegralResult:
    """Integrate over a parametrized region using the metric Jacobian.

    The parametrization need not be full-dimensional. For curves and surfaces,
    the metric Jacobian gives intrinsic length/surface measure.
    """

    transformed, limits, jac_factor = reduce_parametric_integral(
        integrand, ambient_variables, region
    )
    if method not in {"symbolic", "numeric", "auto"}:
        raise ValueError('method must be one of "symbolic", "numeric", or "auto"')
    value = sp.integrate(transformed, *limits)
    exact = True
    if isinstance(value, sp.Integral) or value.has(sp.Integral):
        if method == "symbolic":
            raise NotImplementedError(
                "SymPy could not evaluate the parametric integral symbolically"
            )
        value = sp.Integral(transformed, *limits).evalf(precision)
        exact = False
    elif method == "numeric":
        value = sp.N(value, precision)
        exact = False
    result = ParametricIntegralResult(
        value=sp.simplify(value) if exact else value,
        integrand=sp.sympify(integrand),
        ambient_variables=tuple(
            sp.Symbol(v, real=True) if isinstance(v, str) else v for v in ambient_variables
        ),
        region=region,
        jacobian_factor=jac_factor,
        exact=exact,
    )
    return result if return_result else result.value


__all__ = [
    "ParametricIntegralResult",
    "metric_jacobian_factor",
    "reduce_parametric_integral",
    "integrate_over_parametric_region",
]
