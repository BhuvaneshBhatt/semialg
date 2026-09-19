# ruff: noqa: F403, F405
from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass

import sympy as sp

from ._standard_region_geometry import (
    _effective_parametric_data,
    _independent_vectors,
    _PointData,
    _provably_zero,
    _validate_interval,
)
from ._standard_regions_base import *
from ._standard_regions_curved import *
from ._standard_regions_polyhedral import *
from .exact_arithmetic import compare_exact_reals
from .normalization import normalize_symbol_sequence


@dataclass(frozen=True)
class ParametricRegion(StandardRegion):
    parameters: tuple[sp.Symbol, ...]
    limits: tuple[tuple[sp.Symbol, sp.Expr, sp.Expr], ...]
    mapping: tuple[sp.Expr, ...]
    multiplicity: sp.Expr = sp.Integer(1)
    assumptions: sp.Expr = sp.true

    def __init__(
        self,
        parameters: Sequence[sp.Symbol | str],
        limits: Sequence[tuple[sp.Symbol | str, object, object]],
        mapping: Sequence[object],
        *,
        multiplicity: object = 1,
        assumptions: object = True,
    ):
        """Validate the parameter map, limits, multiplicity, and assumptions before freezing the region."""
        params = normalize_symbol_sequence(parameters)
        if len(set(params)) != len(params):
            raise ValueError("parametric region parameters must be unique")
        by_name = {param.name: param for param in params}
        if len(by_name) != len(params):
            raise ValueError("parametric region parameter names must be unique")
        sym_limits: list[tuple[sp.Symbol, sp.Expr, sp.Expr]] = []
        seen: set[sp.Symbol] = set()
        for raw, lo, hi in limits:
            if isinstance(raw, str):
                param = by_name.get(raw)
                if param is None:
                    raise ValueError(f"limit variable {raw!r} is not a declared parameter")
            else:
                param = raw
                if param not in params:
                    raise ValueError(f"limit variable {param!r} is not a declared parameter")
            if param in seen:
                raise ValueError(f"duplicate integration limit for parameter {param!r}")
            seen.add(param)
            lower = sp.sympify(lo)
            upper = sp.sympify(hi)
            _validate_interval(lower, upper, label=f"limits for {param}")
            sym_limits.append((param, lower, upper))
        if seen != set(params):
            missing = tuple(param for param in params if param not in seen)
            raise ValueError(f"missing integration limits for parameters {missing!r}")
        mult = sp.sympify(multiplicity)
        try:
            mult_cmp = compare_exact_reals(mult, sp.Integer(0))
        except (TypeError, ValueError, NotImplementedError):
            if mult.is_positive is True:
                mult_cmp = 1
            elif mult.is_nonpositive is True:
                mult_cmp = 0
            else:
                raise ValueError("parametrization multiplicity must be provably positive") from None
        if mult_cmp <= 0:
            raise ValueError("parametrization multiplicity must be positive")
        object.__setattr__(self, "parameters", params)
        object.__setattr__(self, "limits", tuple(sym_limits))
        object.__setattr__(self, "mapping", tuple(sp.sympify(expr) for expr in mapping))
        object.__setattr__(self, "multiplicity", mult)
        object.__setattr__(self, "assumptions", sp.sympify(assumptions))

    def dimension(self) -> int:
        free, _limits, mapping, assumptions = _effective_parametric_data(
            self.parameters, self.limits, self.mapping, self.assumptions
        )
        if assumptions not in (True, sp.true):
            from .regions.operations import region_dimension
            from .symbolic_regions import as_semialgebraic_region

            output = tuple(sp.Dummy(f"image{i + 1}", real=True) for i in range(len(mapping)))
            return region_dimension(as_semialgebraic_region(self, output), output)
        if not free:
            return 0
        return int(sp.Matrix(mapping).jacobian(free).rank())

    def ambient_dimension(self) -> int:
        return len(self.mapping)

    def generic_map_degree(self):
        """Return the generic algebraic fiber degree of this parametrization."""

        from .map_degree import parametric_map_degree

        return parametric_map_degree(self.mapping, self.parameters)


@dataclass(frozen=True)
class TransformedRegion(StandardRegion):
    base: StandardRegion
    mapping: tuple[sp.Expr, ...]
    base_variables: tuple[sp.Symbol, ...]

    def __init__(
        self,
        base: StandardRegion,
        mapping: Sequence[object],
        base_variables: Sequence[sp.Symbol | str],
    ):
        if not isinstance(base, StandardRegion):
            raise TypeError("TransformedRegion base must be a StandardRegion")
        variables = normalize_symbol_sequence(base_variables)
        if len(variables) != base.ambient_dimension():
            raise ValueError("base variable count must match the base ambient dimension")
        object.__setattr__(self, "base", base)
        object.__setattr__(self, "mapping", tuple(sp.sympify(e) for e in mapping))
        object.__setattr__(self, "base_variables", variables)

    def dimension(self) -> int:
        """Compute the exact image dimension from the base tangent space and mapping Jacobian."""
        base_dim = self.base.dimension()
        if base_dim <= 0:
            return max(base_dim, 0)

        substitutions: dict[sp.Symbol, sp.Expr] = {}
        tangent_vectors: tuple[_PointData, ...] | None = None
        if isinstance(self.base, IntervalRegion):
            if self.base.dimension() == 0:
                return 0
            tangent_vectors = ((sp.Integer(1),),)
        elif isinstance(self.base, BoxRegion):
            vectors: list[_PointData] = []
            for index, (variable, (lower, upper)) in enumerate(
                zip(self.base_variables, self.base.bounds, strict=True)
            ):
                if _provably_zero(upper - lower):
                    substitutions[variable] = lower
                else:
                    vector = tuple(
                        sp.Integer(1 if i == index else 0) for i in range(len(self.base_variables))
                    )
                    vectors.append(vector)
            tangent_vectors = tuple(vectors)
        elif isinstance(self.base, SimplexRegion):
            anchor = self.base.vertices[0]
            substitutions = dict(zip(self.base_variables, anchor, strict=True))
            differences = tuple(
                tuple(sp.Matrix(vertex) - sp.Matrix(anchor)) for vertex in self.base.vertices[1:]
            )
            tangent_vectors = _independent_vectors(differences)
        elif isinstance(self.base, (ParallelogramRegion, ParallelepipedRegion)):
            substitutions = dict(zip(self.base_variables, self.base.origin, strict=True))
            tangent_vectors = _independent_vectors(self.base.vectors)
        elif isinstance(self.base, PointRegion):
            return 0

        if tangent_vectors is not None:
            if not tangent_vectors:
                return 0
            params = tuple(sp.Dummy(f"u{i + 1}", real=True) for i in range(len(tangent_vectors)))
            base_point = tuple(substitutions.get(v, v) for v in self.base_variables)
            restricted = []
            for expr in self.mapping:
                replacement = {}
                for coord, variable in enumerate(self.base_variables):
                    value = base_point[coord] + sum(
                        param * vector[coord]
                        for param, vector in zip(params, tangent_vectors, strict=True)
                    )
                    replacement[variable] = value
                restricted.append(sp.simplify(expr.subs(replacement)))
            return int(sp.Matrix(restricted).jacobian(params).rank())

        from .regions.operations import region_dimension
        from .symbolic_regions import as_semialgebraic_region

        output = tuple(sp.Dummy(f"image{i + 1}", real=True) for i in range(len(self.mapping)))
        return region_dimension(as_semialgebraic_region(self, output), output)

    def ambient_dimension(self) -> int:
        return len(self.mapping)


@dataclass(frozen=True)
class BooleanRegion(StandardRegion):
    op: str
    regions: tuple[StandardRegion, ...]
    assume_disjoint: bool = False

    def __init__(
        self, op: str, regions: Sequence[StandardRegion], *, assume_disjoint: bool = False
    ):
        if op not in {"union", "intersection", "difference", "symmetric_difference", "complement"}:
            raise ValueError("unsupported BooleanRegion op")
        region_tuple = tuple(regions)
        if any(not isinstance(region, StandardRegion) for region in region_tuple):
            raise TypeError("BooleanRegion members must be StandardRegion objects")
        if op in {"difference", "symmetric_difference"} and len(region_tuple) != 2:
            raise ValueError(f"{op} requires exactly two regions")
        if op == "complement" and len(region_tuple) != 1:
            raise ValueError("complement requires exactly one region")
        ambient_dims = {region.ambient_dimension() for region in region_tuple}
        if len(ambient_dims) > 1:
            raise ValueError("Boolean regions must share one ambient dimension")
        object.__setattr__(self, "op", op)
        object.__setattr__(self, "regions", region_tuple)
        object.__setattr__(self, "assume_disjoint", assume_disjoint)

    def dimension(self) -> int:
        if not self.regions:
            return -1
        if self.op == "union":
            return max(r.dimension() for r in self.regions)
        if self.op == "intersection" and all(
            isinstance(region, IntervalRegion) for region in self.regions
        ):
            intervals = tuple(
                region for region in self.regions if isinstance(region, IntervalRegion)
            )
            lower = intervals[0].lower
            upper = intervals[0].upper
            try:
                for interval in intervals[1:]:
                    if compare_exact_reals(interval.lower, lower) > 0:
                        lower = interval.lower
                    if compare_exact_reals(interval.upper, upper) < 0:
                        upper = interval.upper
                relation = compare_exact_reals(lower, upper)
            except (TypeError, ValueError, NotImplementedError) as exc:
                raise NotImplementedError(
                    "interval-intersection dimension requires exactly comparable endpoints"
                ) from exc
            if relation > 0:
                return -1
            if relation < 0:
                return 1
            for interval in intervals:
                try:
                    at_lower = compare_exact_reals(lower, interval.lower) == 0
                    at_upper = compare_exact_reals(upper, interval.upper) == 0
                except (TypeError, ValueError, NotImplementedError) as exc:
                    raise NotImplementedError(
                        "interval-intersection dimension requires exactly comparable endpoints"
                    ) from exc
                if (at_lower and not interval.lower_closed) or (
                    at_upper and not interval.upper_closed
                ):
                    return -1
            return 0
        raise NotImplementedError(
            "exact dimension is not structurally determined for this BooleanRegion operation; "
            "convert it to SemialgebraicRegion and use region_dimension()"
        )

    def ambient_dimension(self) -> int:
        return self.regions[0].ambient_dimension() if self.regions else 0


def RegionUnion(*regions: StandardRegion, assume_disjoint: bool = False) -> BooleanRegion:
    """Return the Boolean union of standard regions."""
    return BooleanRegion("union", regions, assume_disjoint=assume_disjoint)


def RegionIntersection(*regions: StandardRegion) -> BooleanRegion:
    """Return the Boolean intersection of standard regions."""
    return BooleanRegion("intersection", regions)


def RegionDifference(a: StandardRegion, b: StandardRegion) -> BooleanRegion:
    """Return the Boolean difference of two standard regions."""
    return BooleanRegion("difference", (a, b))


def RegionSymmetricDifference(a: StandardRegion, b: StandardRegion) -> BooleanRegion:
    """Return the Boolean symmetric difference of two standard regions."""
    return BooleanRegion("symmetric_difference", (a, b))


def is_standard_region(obj: object) -> bool:
    return isinstance(obj, StandardRegion)
