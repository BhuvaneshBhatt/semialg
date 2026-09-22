"""Symbolic membership and relation predicates for semialgebraic regions."""

from __future__ import annotations

from collections.abc import Sequence

import sympy as sp
from sympy.logic.boolalg import Boolean, BooleanFunction, as_Boolean

from .normalization import normalize_formula, normalize_variables
from .quantifiers import Exists, ForAll
from .structural_keys import ordered_symbols
from .symbolic_regions import (
    SemialgebraicRegion,
    _fresh_symbols,
    _point_tuple,
    _require_same_ambient,
)


def as_semialgebraic_region(region, variables=None):
    """Coerce a region lazily so direct module imports remain acyclic."""

    from .region_coercion import as_semialgebraic_region as coerce_region

    return coerce_region(region, variables)


class RegionElement(Boolean):
    """Symbolic assertion that ``point`` belongs to ``region``."""

    def __new__(cls, point: object, region: object):
        reg = as_semialgebraic_region(region)
        pt = _point_tuple(point)
        if len(pt) != reg.ambient_dimension:
            raise ValueError("point dimension does not match region ambient dimension")
        return sp.Basic.__new__(cls, sp.Tuple(*pt), reg)

    @property
    def point(self) -> tuple[sp.Expr, ...]:
        return tuple(self.args[0])

    @property
    def region(self) -> SemialgebraicRegion:
        return self.args[1]

    def as_formula(self, *, eliminate: bool = False) -> sp.Expr:
        return self.region.membership_formula(self.point, eliminate=eliminate)

    def evaluate(self, *, strategy: str = "auto") -> bool:
        return self.region.contains(self.point, strategy=strategy)


class RegionNotElement(Boolean):
    """Symbolic assertion that ``point`` does not belong to ``region``."""

    def __new__(cls, point: object, region: object):
        element = RegionElement(point, region)
        return sp.Basic.__new__(cls, element.args[0], element.args[1])

    @property
    def point(self) -> tuple[sp.Expr, ...]:
        return tuple(self.args[0])

    @property
    def region(self) -> SemialgebraicRegion:
        return self.args[1]

    def as_formula(self, *, eliminate: bool = False) -> sp.Expr:
        return normalize_formula(
            sp.Not(self.region.membership_formula(self.point, eliminate=eliminate))
        )

    def evaluate(self, *, strategy: str = "auto") -> bool:
        return not self.region.contains(self.point, strategy=strategy)


class _RegionRelation(Boolean):
    __slots__ = ()

    @property
    def regions(self) -> tuple[SemialgebraicRegion, ...]:
        return tuple(self.args)

    @property
    def ambient_dimension(self) -> int:
        return self.regions[0].ambient_dimension if self.regions else 0

    def evaluate(self, *, strategy: str | None = None) -> bool:
        from .reasoning_regions import region_disjoint, region_equal, region_subset

        regs = self.regions
        for region in regs:
            parameters = region.quantifier_free_formula().free_symbols - set(region.variables)
            if parameters:
                raise ValueError(
                    "parameter-dependent region relation is symbolic; use "
                    "region_relation_conditions(..., eliminate=True)"
                )
        if isinstance(self, RegionSubset):
            return bool(
                region_subset(
                    regs[0].quantifier_free_formula(),
                    regs[1].quantifier_free_formula(),
                    regs[0].variables,
                    strategy=strategy,
                )
            )
        if isinstance(self, RegionDisjoint):
            for i, left in enumerate(regs):
                for right in regs[i + 1 :]:
                    if not region_disjoint(
                        left.quantifier_free_formula(),
                        right.quantifier_free_formula(),
                        left.variables,
                        strategy=strategy,
                    ):
                        return False
            return True
        if isinstance(self, RegionEqual):
            if len(regs) < 2:
                return True
            first = regs[0]
            return all(
                region_equal(
                    first.quantifier_free_formula(),
                    other.quantifier_free_formula(),
                    first.variables,
                    strategy=strategy,
                )
                for other in regs[1:]
            )
        raise TypeError(type(self).__name__)


class RegionSubset(_RegionRelation):
    """Symbolic region-subset relation."""

    def __new__(cls, lhs: object, rhs: object):
        left = as_semialgebraic_region(lhs)
        right = as_semialgebraic_region(rhs, left.variables)
        _require_same_ambient((left, right))
        return sp.Basic.__new__(cls, left, right)

    def as_formula(self, *, quantifier_free_regions: bool = False) -> sp.Expr:
        left, right = self.regions
        coords = _fresh_symbols("r", left.ambient_dimension)
        a = left.membership_formula(coords, eliminate=quantifier_free_regions)
        b = right.membership_formula(coords, eliminate=quantifier_free_regions)
        return ForAll(coords, sp.Or(sp.Not(a), b))


class RegionDisjoint(_RegionRelation):
    """Symbolic assertion that all supplied regions are pairwise disjoint."""

    def __new__(cls, *regions: object):
        regs = tuple(as_semialgebraic_region(region) for region in regions)
        if regs:
            base_vars = regs[0].variables
            regs = tuple(as_semialgebraic_region(region, base_vars) for region in regs)
            _require_same_ambient(regs)
        return sp.Basic.__new__(cls, *regs)

    def as_formula(self, *, quantifier_free_regions: bool = False) -> sp.Expr:
        if len(self.regions) < 2:
            return sp.true
        coords = _fresh_symbols("r", self.ambient_dimension)
        clauses = []
        for i, left in enumerate(self.regions):
            for right in self.regions[i + 1 :]:
                clauses.append(
                    sp.Not(
                        Exists(
                            coords,
                            sp.And(
                                left.membership_formula(coords, eliminate=quantifier_free_regions),
                                right.membership_formula(coords, eliminate=quantifier_free_regions),
                            ),
                        )
                    )
                )
        return sp.And(*clauses)


class RegionEqual(_RegionRelation):
    """Symbolic assertion that all supplied regions are equal."""

    def __new__(cls, *regions: object):
        regs = tuple(as_semialgebraic_region(region) for region in regions)
        if regs:
            base_vars = regs[0].variables
            regs = tuple(as_semialgebraic_region(region, base_vars) for region in regs)
            _require_same_ambient(regs)
        return sp.Basic.__new__(cls, *regs)

    def as_formula(self, *, quantifier_free_regions: bool = False) -> sp.Expr:
        if len(self.regions) < 2:
            return sp.true
        first = self.regions[0]
        return sp.And(
            *(
                RegionSubset(first, other).as_formula(
                    quantifier_free_regions=quantifier_free_regions
                )
                for other in self.regions[1:]
            ),
            *(
                RegionSubset(other, first).as_formula(
                    quantifier_free_regions=quantifier_free_regions
                )
                for other in self.regions[1:]
            ),
        )


def _eliminate_relation_formula(
    expr: sp.Expr, free_variables: Sequence[sp.Symbol] | None = None
) -> sp.Expr:
    from .formula import parse_quantified_expr
    from .quantifiers import split_quantifiers
    from .solve.reduce import reduce_formula

    order = (
        tuple(free_variables) if free_variables is not None else ordered_symbols(expr.free_symbols)
    )
    prefix, _ = split_quantifiers(expr)
    bound = tuple(variable for _, variable in prefix if variable not in order)
    parsed = parse_quantified_expr(expr, variable_order=order + bound)
    return normalize_formula(reduce_formula(parsed, strategy="auto"))


def _region_parameters(expr: object) -> tuple[sp.Symbol, ...]:
    """Return algebraic-level parameters carried by symbolic region atoms."""

    parameters: set[sp.Symbol] = set()
    for node in sp.preorder_traversal(as_Boolean(expr)):
        regions: tuple[SemialgebraicRegion, ...] = ()
        if isinstance(node, (RegionElement, RegionNotElement)):
            regions = (node.region,)
        elif isinstance(node, _RegionRelation):
            regions = node.regions
        for region in regions:
            parameters.update(region.formula.free_symbols - set(region.variables))
    return ordered_symbols(parameters)


def _real_parameter_formula(parameters: Sequence[sp.Symbol]) -> sp.Expr:
    """Build explicit real-domain membership conditions without auto-evaluation."""

    if not parameters:
        return sp.true
    return sp.And(*(sp.Contains(parameter, sp.S.Reals, evaluate=False) for parameter in parameters))


def _contains_quantifiers(expr: sp.Expr) -> bool:
    return any(isinstance(node, (Exists, ForAll)) for node in sp.preorder_traversal(expr))


def region_element_conditions(
    expr: object,
    *,
    eliminate: bool = False,
    real_parameters: bool = False,
) -> sp.Expr:
    """Lower symbolic region-membership predicates inside a Boolean expression.

    By default this is a structural, non-CAD transformation: membership formulas
    that still contain unresolved quantifiers are retained as symbolic region atoms.
    Set ``eliminate=True`` to request exact quantifier elimination.  Set
    ``real_parameters=True`` to conjoin explicit ``Contains(parameter, Reals)``
    conditions for algebraic-level parameters carried by the referenced regions.
    """

    use_elimination = bool(eliminate)

    def lower(value: object) -> sp.Expr:
        if isinstance(value, RegionElement):
            formula = value.as_formula(eliminate=use_elimination)
            if not eliminate and _contains_quantifiers(formula):
                return value
            return formula
        if isinstance(value, RegionNotElement):
            formula = value.as_formula(eliminate=use_elimination)
            if not eliminate and _contains_quantifiers(formula):
                return value
            return formula
        if isinstance(value, sp.Not) and isinstance(value.args[0], RegionElement):
            formula = normalize_formula(sp.Not(value.args[0].as_formula(eliminate=use_elimination)))
            if not eliminate and _contains_quantifiers(formula):
                return value
            return formula
        if isinstance(value, (Exists, ForAll)):
            return value.func(value.variables, lower(value.formula))
        if isinstance(value, BooleanFunction):
            return value.func(*(lower(arg) for arg in value.args), evaluate=False)
        return as_Boolean(value)

    lowered = lower(expr)
    if real_parameters:
        lowered = normalize_formula(
            sp.And(_real_parameter_formula(_region_parameters(expr)), lowered)
        )
    return lowered


def region_relation_conditions(
    expr: object,
    variables: Sequence[sp.Symbol | str] | None = None,
    *,
    eliminate: bool = False,
    real_parameters: bool = False,
) -> sp.Expr:
    """Lower ``RegionSubset``/``RegionDisjoint``/``RegionEqual`` to quantified conditions.

    With ``eliminate=True`` the generated real quantifiers are eliminated using
    semialg's exact QE stack. ``variables`` specifies the desired free-variable
    order for parameter-dependent relations. ``real_parameters=True`` additionally
    records explicit real-domain membership for algebraic-level region parameters.
    """

    def lower(value: object) -> sp.Expr:
        if isinstance(value, _RegionRelation):
            formula = value.as_formula(quantifier_free_regions=eliminate)
            if eliminate:
                free = (
                    tuple(normalize_variables(variables, formula, append_context_symbols=False))
                    if variables is not None
                    else ordered_symbols(formula.free_symbols)
                )
                return _eliminate_relation_formula(formula, free)
            return formula
        if isinstance(value, (Exists, ForAll)):
            return value.func(value.variables, lower(value.formula))
        if isinstance(value, BooleanFunction):
            return normalize_formula(
                value.func(*(lower(arg) for arg in value.args), evaluate=False)
            )
        return as_Boolean(value)

    lowered = lower(expr)
    if real_parameters:
        lowered = normalize_formula(
            sp.And(_real_parameter_formula(_region_parameters(expr)), lowered)
        )
    return lowered


__all__ = [
    "RegionElement",
    "RegionNotElement",
    "RegionSubset",
    "RegionDisjoint",
    "RegionEqual",
    "region_element_conditions",
    "region_relation_conditions",
]
