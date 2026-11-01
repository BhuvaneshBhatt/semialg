"""Support-function certificates for affine objectives on exact ellipsoids."""

import sympy as sp

from .exact_arithmetic import exact_sign
from .normalization import require_conjunction
from .optimization_results import OptimizationResult
from .quadratic_forms import positive_quadratic_form
from .relations import split_relation


def linear_ellipsoid_optimum(objective, condition, variables, *, kind):
    """Use Cauchy–Schwarz in the ellipsoid metric, retaining open-set attainment.

    One quadratic relation defines the entire domain. Extra constraints remain
    the responsibility of the general optimizer. Positive definiteness and all
    coefficients are checked exactly, so singular quadrics and parameters do
    not receive a positive-definite support certificate.
    """
    atoms = require_conjunction(condition, message="ellipsoid optimization expects a conjunction")
    if len(atoms) != 1 or not isinstance(atoms[0], sp.core.relational.Relational):
        return None
    try:
        affine = sp.Poly(objective, *variables)
    except sp.PolynomialError:
        return None
    if affine.total_degree() > 1 or any(
        value.free_symbols or value.is_real is not True or value.is_algebraic is not True
        for value in affine.coeffs()
    ):
        return None
    expression, relation = split_relation(atoms[0])
    if relation not in {"<", "<=", ">", ">=", "=="}:
        return None
    if relation in {">", ">="}:
        expression = -expression
    data = positive_quadratic_form(expression, variables)
    if data is None and relation == "==":
        data = positive_quadratic_form(-expression, variables)
    if data is None:
        return None
    form, center, level = data
    try:
        level_sign = exact_sign(level)
    except (ValueError, NotImplementedError):
        return None
    strict = relation in {"<", ">"}
    if level_sign < 0 or (level_sign == 0 and strict):
        raise ValueError("ellipsoid domain is empty")
    coefficients = sp.ImmutableMatrix([affine.coeff_monomial(variable) for variable in variables])
    dual = form.matrix.inv() * coefficients
    norm_squared = sp.cancel((coefficients.T * dual)[0])
    center_value = sp.cancel(objective.subs(dict(zip(variables, center, strict=True))))
    try:
        norm_sign = exact_sign(norm_squared)
    except (ValueError, NotImplementedError):
        return None
    if norm_sign == 0:
        point = sp.MutableDenseMatrix(center)
        if relation == "==" and level_sign > 0:
            point[0] += sp.sqrt(level / form.matrix[0, 0])
        value = center_value
        attained = True
    else:
        direction = -1 if kind == "min" else 1
        point = center + direction * sp.sqrt(level / norm_squared) * dual
        value = sp.simplify(center_value + direction * sp.sqrt(level * norm_squared))
        attained = not strict
    points = (dict(zip(variables, map(sp.simplify, point), strict=True)),) if attained else ()
    return OptimizationResult(
        objective,
        variables,
        value,
        points,
        attained,
        kind,
        "exact_ellipsoid_support",
        {
            "constraint_formula": condition,
            "global_certificate": "positive_definite_ellipsoid_support",
            "matrix": form.matrix,
            "center": center,
            "level": level,
            "dual_norm_squared": norm_squared,
        },
        True,
    )
