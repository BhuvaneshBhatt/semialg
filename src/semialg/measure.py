from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass

import sympy as sp


@dataclass(frozen=True)
class MeasureResult:
    """Exact measure result for a semialgebraic set.

    The ``value`` field is the Lebesgue measure in the ambient variables used
    by the call. ``method`` records the reconstruction strategy that produced
    the result.
    """

    value: sp.Expr
    variables: tuple[sp.Symbol, ...]
    method: str
    diagnostics: Mapping[str, object] | None = None


def semialgebraic_measure(
    condition: object,
    variables: Sequence[sp.Symbol | str],
    *,
    bounds: Sequence[tuple[sp.Symbol | str, object, object]]
    | Mapping[sp.Symbol | str, tuple[object, object]]
    | None = None,
    measure_dimension: object = "ambient",
    return_result: bool = False,
    parameters: Sequence[sp.Symbol | str] | None = None,
    return_stratified: bool = False,
) -> sp.Expr | MeasureResult | object:
    """Return the exact measure of a supported semialgebraic set.

    The measure implementation delegates to the structural region-integral
    reducer with integrand ``1``. This keeps ``semialgebraic_measure`` aligned
    with all standard-shape recognition supported by ``integrate_over_region``:
    one-dimensional cell intervals, axis-aligned boxes, the 2D unit simplex,
    axis-aligned ellipses, origin-centered disks/annuli, simple vertical slices,
    and the first intrinsic-dimensional cases.
    """

    from .region_integrate import integrate_over_region

    result = integrate_over_region(
        1,
        condition,
        variables,
        bounds=bounds,
        measure_dimension=measure_dimension,
        return_result=True,
        parameters=parameters,
        return_stratified=return_stratified,
    )
    if return_stratified:
        return result
    if not hasattr(result, "method"):
        raise TypeError("integration backend returned an unexpected result type")
    method_name = result.method
    if method_name == "one_dimensional_cell_integration":
        method_name = "one_dimensional_cell_sampling"
    measure_result = MeasureResult(
        result.value,
        result.variables,
        method_name,
        {**(result.diagnostics or {}), "delegated_to": "integrate_over_region"},
    )
    return measure_result if return_result else measure_result.value


__all__ = ["MeasureResult", "semialgebraic_measure"]
