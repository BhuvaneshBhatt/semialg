from __future__ import annotations

from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass, field
from itertools import combinations

import sympy as sp


@dataclass(frozen=True)
class ProjectionPolynomial:
    """One polynomial in a projection tower with auditable provenance."""

    poly: sp.Poly
    level: int
    source: str
    parents: tuple[str, ...] = ()
    operation_variable: sp.Symbol | None = None
    expression: sp.Expr | None = None

    @property
    def key(self) -> str:
        return _poly_sort_key(self.poly)[0]


@dataclass(frozen=True)
class ProjectionLevel:
    """Projection data for one CAD level."""

    level: int
    variable: sp.Symbol | None
    polynomials: tuple[sp.Poly, ...]
    entries: tuple[ProjectionPolynomial, ...] = ()

    def __post_init__(self) -> None:
        if not self.entries:
            entries = tuple(
                ProjectionPolynomial(poly=poly, level=self.level, source="unspecified")
                for poly in self.polynomials
            )
            object.__setattr__(self, "entries", entries)

    @property
    def provenance_by_key(self) -> dict[str, ProjectionPolynomial]:
        return {entry.key: entry for entry in self.entries}


@dataclass(frozen=True)
class ProjectionTower:
    """Complete Collins projection tower."""

    variables: tuple[sp.Symbol, ...]
    levels: tuple[ProjectionLevel, ...]
    original_polynomials: tuple[sp.Poly, ...]
    metadata: Mapping[str, object] = field(default_factory=dict)

    def level(self, level: int) -> ProjectionLevel:
        return self.levels[level - 1]

    @property
    def by_level(self) -> dict[int, tuple[sp.Poly, ...]]:
        return {item.level: item.polynomials for item in self.levels}

    def entries_by_level(self) -> dict[int, tuple[ProjectionPolynomial, ...]]:
        return {item.level: item.entries for item in self.levels}

    def poly_count_by_level(self) -> dict[int, int]:
        return {item.level: len(item.polynomials) for item in self.levels}


def _poly_sort_key(poly: sp.Poly) -> tuple[str, tuple[str, ...]]:
    return (sp.sstr(sp.expand(poly.as_expr())), tuple(sp.sstr(gen) for gen in poly.gens))


def normalize_poly(poly: sp.Poly) -> sp.Poly | None:
    if poly.is_zero:
        return None
    primitive = poly.primitive()[1]
    if primitive.is_zero:
        return None
    if primitive.LC().could_extract_minus_sign():
        primitive = -primitive
    try:
        factors = sp.factor_list(primitive.as_expr(), *primitive.gens)[1]
    except Exception:
        factors = [(primitive.as_expr(), 1)]
    pieces: list[sp.Poly] = []
    for factor, _mult in factors:
        factor_poly = sp.Poly(factor, *primitive.gens)
        if factor_poly.total_degree() > 0:
            if factor_poly.LC().could_extract_minus_sign():
                factor_poly = -factor_poly
            pieces.append(factor_poly)
    if not pieces:
        return None
    result = pieces[0]
    for piece in pieces[1:]:
        result *= piece
    result = result.primitive()[1]
    if result.LC().could_extract_minus_sign():
        result = -result
    return result


def _squarefree_factors(poly: sp.Poly) -> tuple[sp.Poly, ...]:
    normalized = normalize_poly(poly)
    if normalized is None or normalized.total_degree() == 0:
        return tuple()
    factors: list[sp.Poly] = [normalized]
    try:
        raw_factors = sp.factor_list(normalized.as_expr(), *normalized.gens)[1]
    except Exception:
        raw_factors = []
    for factor_expr, _mult in raw_factors:
        factor = normalize_poly(sp.Poly(factor_expr, *normalized.gens))
        if factor is not None and factor.total_degree() > 0:
            factors.append(factor)
    return tuple(factors)


def squarefree_basis(polys: Iterable[sp.Poly]) -> tuple[sp.Poly, ...]:
    seen: set[tuple[str, tuple[str, ...]]] = set()
    out: list[sp.Poly] = []
    for poly in polys:
        for factor in _squarefree_factors(poly):
            key = _poly_sort_key(factor)
            if key not in seen:
                seen.add(key)
                out.append(factor)
    return tuple(sorted(out, key=_poly_sort_key))


def _as_poly(expr: sp.Expr, gens: Sequence[sp.Symbol]) -> sp.Poly | None:
    expr = sp.expand(expr)
    if not gens:
        return None
    try:
        return normalize_poly(sp.Poly(expr, *gens))
    except Exception:
        return None


def _entry(
    poly: sp.Poly | None,
    *,
    level: int,
    source: str,
    parents: tuple[str, ...],
    var: sp.Symbol | None,
    expr: sp.Expr | None = None,
) -> ProjectionPolynomial | None:
    if poly is None:
        return None
    return ProjectionPolynomial(
        poly=poly,
        level=level,
        source=source,
        parents=parents,
        operation_variable=var,
        expression=expr,
    )


def _dedupe_entries(entries: Iterable[ProjectionPolynomial]) -> tuple[ProjectionPolynomial, ...]:
    by_key: dict[tuple[str, tuple[str, ...]], ProjectionPolynomial] = {}
    for item in entries:
        key = _poly_sort_key(item.poly)
        if key not in by_key:
            by_key[key] = item
            continue
        old = by_key[key]
        source = old.source if item.source == old.source else f"{old.source}|{item.source}"
        parents = tuple(dict.fromkeys((*old.parents, *item.parents)))
        by_key[key] = ProjectionPolynomial(
            poly=old.poly,
            level=old.level,
            source=source,
            parents=parents,
            operation_variable=old.operation_variable,
            expression=old.expression,
        )
    return tuple(sorted(by_key.values(), key=lambda entry: _poly_sort_key(entry.poly)))


def _content(
    poly: sp.Poly, var: sp.Symbol, lower_gens: Sequence[sp.Symbol], level: int
) -> list[ProjectionPolynomial]:
    content = sp.Poly(poly.as_expr(), var).content()
    item = _entry(
        _as_poly(content, lower_gens),
        level=level,
        source="content",
        parents=(_poly_sort_key(poly)[0],),
        var=var,
        expr=content,
    )
    return [] if item is None else [item]


def _coefficients(
    poly: sp.Poly, var: sp.Symbol, lower_gens: Sequence[sp.Symbol], level: int
) -> list[ProjectionPolynomial]:
    out: list[ProjectionPolynomial] = []
    parent = _poly_sort_key(poly)[0]
    for coeff in sp.Poly(poly.as_expr(), var).all_coeffs():
        projected = _as_poly(coeff, lower_gens)
        item = _entry(
            projected, level=level, source="coefficient", parents=(parent,), var=var, expr=coeff
        )
        if item is not None:
            out.append(item)
    return out


def _discriminant(
    poly: sp.Poly, var: sp.Symbol, lower_gens: Sequence[sp.Symbol], level: int
) -> list[ProjectionPolynomial]:
    if poly.degree(var) <= 1:
        return []
    expr = sp.discriminant(poly.as_expr(), var)
    item = _entry(
        _as_poly(expr, lower_gens),
        level=level,
        source="discriminant",
        parents=(_poly_sort_key(poly)[0],),
        var=var,
        expr=expr,
    )
    return [] if item is None else [item]


def _resultant(
    left: sp.Poly, right: sp.Poly, var: sp.Symbol, lower_gens: Sequence[sp.Symbol], level: int
) -> list[ProjectionPolynomial]:
    expr = sp.resultant(left.as_expr(), right.as_expr(), var)
    item = _entry(
        _as_poly(expr, lower_gens),
        level=level,
        source="resultant",
        parents=(_poly_sort_key(left)[0], _poly_sort_key(right)[0]),
        var=var,
        expr=expr,
    )
    return [] if item is None else [item]


def collins_proj_entries(
    polys: Sequence[sp.Poly], var: sp.Symbol, lower_gens: Sequence[sp.Symbol], level: int
) -> tuple[ProjectionPolynomial, ...]:
    basis = squarefree_basis(polys)
    active = [poly for poly in basis if poly.degree(var) > 0]
    inactive = [poly for poly in basis if poly.degree(var) == 0]
    projected: list[ProjectionPolynomial] = []
    # Polynomials that do not contain the eliminated variable remain constraints
    # on the lower-dimensional base and must be carried downward. Dropping
    # them loses parameter-only bounds such as x <= 4 in exists y. y^2 = x.
    if lower_gens:
        for poly in inactive:
            item = _entry(
                _as_poly(poly.as_expr(), lower_gens),
                level=level,
                source="inactive",
                parents=(_poly_sort_key(poly)[0],),
                var=var,
                expr=poly.as_expr(),
            )
            if item is not None:
                projected.append(item)
    for poly in active:
        projected.extend(_content(poly, var, lower_gens, level))
        projected.extend(_coefficients(poly, var, lower_gens, level))
        projected.extend(_discriminant(poly, var, lower_gens, level))
    for left, right in combinations(active, 2):
        projected.extend(_resultant(left, right, var, lower_gens, level))
    return _dedupe_entries(projected)


def collins_projection_step(
    polys: Sequence[sp.Poly], var: sp.Symbol, lower_gens: Sequence[sp.Symbol]
) -> tuple[sp.Poly, ...]:
    return tuple(
        entry.poly for entry in collins_proj_entries(polys, var, lower_gens, len(lower_gens))
    )


def _input_entries(polys: Sequence[sp.Poly], level: int) -> tuple[ProjectionPolynomial, ...]:
    return _dedupe_entries(
        ProjectionPolynomial(poly=poly, level=level, source="input")
        for poly in squarefree_basis(polys)
    )


def build_collins_proj_set(
    polys: Sequence[sp.Expr | sp.Poly], variables: Sequence[sp.Symbol]
) -> ProjectionTower:
    vars_tuple = tuple(variables)
    if not vars_tuple:
        raise ValueError("at least one variable is required")
    current = squarefree_basis(
        poly if isinstance(poly, sp.Poly) else sp.Poly(sp.expand(poly), *vars_tuple)
        for poly in polys
    )
    poly_levels: dict[int, tuple[sp.Poly, ...]] = {len(vars_tuple): current}
    entry_levels: dict[int, tuple[ProjectionPolynomial, ...]] = {
        len(vars_tuple): _input_entries(current, len(vars_tuple))
    }
    for level in range(len(vars_tuple), 1, -1):
        entries = collins_proj_entries(
            current, vars_tuple[level - 1], vars_tuple[: level - 1], level - 1
        )
        current = tuple(entry.poly for entry in entries)
        poly_levels[level - 1] = current
        entry_levels[level - 1] = entries
    levels = tuple(
        ProjectionLevel(
            i,
            vars_tuple[i - 1],
            poly_levels.get(i, ()),
            entry_levels.get(i, ()),
        )
        for i in range(1, len(vars_tuple) + 1)
    )
    return ProjectionTower(
        variables=vars_tuple,
        levels=levels,
        original_polynomials=poly_levels[len(vars_tuple)],
        metadata={"projection": "collins", "complete": True, "provenance": True},
    )
