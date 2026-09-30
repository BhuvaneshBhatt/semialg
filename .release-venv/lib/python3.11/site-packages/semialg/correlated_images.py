"""Certified correlated images of polynomial/rational maps."""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass, field

import sympy as sp

from .geometry_queries import _mapping_tuple, semialgebraic_image
from .normalization import (
    normalize_formula,
    normalize_problem_variables,
    normalize_symbol_sequence,
)
from .parameterization_analysis import parameterization_geometry


@dataclass(frozen=True)
class CorrelatedMapImageResult:
    mapping: tuple[sp.Expr, ...]
    variables: tuple[sp.Symbol, ...]
    domain: sp.Expr
    image_variables: tuple[sp.Symbol, ...]
    formula: sp.Expr
    method: str = "complete_qe_correlated_image"
    certified: bool = True
    diagnostics: Mapping[str, object] = field(default_factory=dict)

    @property
    def proof_trace(self) -> tuple[str, ...]:
        steps = ["exact_graph", "existential_projection", self.method]
        if self.diagnostics.get("sphere_constraint"):
            steps.insert(0, "sphere_constraint")
        if self.diagnostics.get("critical_value_set") not in (None, sp.false, "False"):
            steps.append("critical_value_geometry")
        return tuple(steps)


def _targets(count: int, forbidden: set[sp.Symbol]) -> tuple[sp.Symbol, ...]:
    out: list[sp.Symbol] = []
    i = 0
    while len(out) < count:
        target = sp.Symbol(f"_image_{i}", real=True)
        i += 1
        if target not in forbidden:
            out.append(target)
    return tuple(out)


def correlated_map_image(mapping, domain=sp.true, variables=None, *, image_variables=None):
    """Return the exact *joint* image of a polynomial/rational map.

    Unlike coordinate-wise range calls, this preserves correlations between
    output coordinates.  Completeness follows from exact graph projection.
    """
    maps = _mapping_tuple(mapping)
    condition = normalize_formula(domain)
    vars_ = normalize_problem_variables(variables, sp.Tuple(*maps), condition)
    if image_variables is None:
        mapping_symbols = set().union(*(expr.free_symbols for expr in maps))
        targets = _targets(len(maps), set(vars_) | condition.free_symbols | mapping_symbols)
    else:
        targets = normalize_symbol_sequence(image_variables)
    if len(targets) != len(maps):
        raise ValueError("image_variables must have the same length as mapping")
    if set(targets) & set(vars_):
        raise ValueError("image variables must be distinct from source variables")
    formula = semialgebraic_image(maps, condition, vars_, image_variables=targets)
    diagnostics = {
        "output_dimension": len(maps),
        "source_dimension": len(vars_),
        "correlated": True,
    }
    try:
        geom = parameterization_geometry(maps, condition, vars_, image_variables=targets)
        diagnostics.update(
            {
                "generic_rank": geom.generic_rank,
                "image_dimension": geom.image_dimension,
                "critical_value_set": geom.critical_value_set,
            }
        )
    except (ArithmeticError, TypeError, ValueError, NotImplementedError, sp.PolynomialError):
        pass
    return CorrelatedMapImageResult(
        maps, vars_, condition, targets, sp.simplify(formula), diagnostics=diagnostics
    )


def angular_map_image(mapping, variables, domain=sp.true, *, image_variables=None, sphere=True):
    """Return the exact correlated image of a map on an angular domain.

    With ``sphere=True`` the unit-sphere equation is added exactly.  Rational
    map denominators are handled by :func:`semialgebraic_image`'s graph logic.
    """
    vars_ = tuple(variables)
    condition = normalize_formula(domain)
    if sphere:
        condition = sp.And(condition, sp.Eq(sum(v**2 for v in vars_), 1))
    result = correlated_map_image(mapping, condition, vars_, image_variables=image_variables)
    diag = dict(result.diagnostics)
    diag["sphere_constraint"] = bool(sphere)
    return CorrelatedMapImageResult(
        result.mapping,
        result.variables,
        result.domain,
        result.image_variables,
        result.formula,
        "complete_qe_angular_map_image",
        True,
        diag,
    )


__all__ = ["CorrelatedMapImageResult", "correlated_map_image", "angular_map_image"]
