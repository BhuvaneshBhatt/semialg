"""Exact structure and simplification for polynomial equality ideals.

This module centralizes Groebner-basis computations used by the exact solver.
It intentionally exposes algebraic structure rather than a second solver API:
callers use :class:`EqualityIdealContext` to inspect dimension, quotient
algebra size, coordinate eliminants, certified substitutions, and radical
consequences of an equality subsystem.
"""

from __future__ import annotations

from collections.abc import Iterable, Sequence
from dataclasses import dataclass

import sympy as sp

from ..internal_symbols import fresh_dummy
from .groebner_utils import compute_groebner_basis
from .rational_univariate._domains import require_exact_rur_domain
from .rational_univariate.quotient import (
    _leading_exponent_grevlex,
    _standard_exponent_count,
    _standard_exponents,
)
from .rational_univariate.representation import RationalUnivariateError


@dataclass(frozen=True)
class LinearVariableElimination:
    """A certified polynomial replacement for one variable modulo an ideal."""

    variable: sp.Symbol
    replacement: sp.Expr
    remaining_generators: tuple[sp.Expr, ...]
    certificate: sp.Expr

    @property
    def remaining_basis(self) -> tuple[sp.Expr, ...]:
        return self.remaining_generators


@dataclass(frozen=True)
class ZeroDimensionalFilterResult:
    """Exact real points from a finite equality variety after constraints."""

    variables: tuple[sp.Symbol, ...]
    points: tuple[tuple[sp.Expr, ...], ...]
    quotient_dimension: int | None
    coordinate_polynomials: tuple[sp.Expr, ...]
    representation: object | None = None

    @property
    def satisfiable(self) -> bool:
        return bool(self.points)


@dataclass(frozen=True)
class EqualityIdealAnalysis:
    """Exact structural invariants of a polynomial equality ideal."""

    variables: tuple[sp.Symbol, ...]
    generators: tuple[sp.Expr, ...]
    groebner_basis: tuple[sp.Expr, ...]
    dimension: int
    quotient_dimension: int | None
    zero_dimensional: bool
    inconsistent: bool


def _normalize_generators(
    generators: Iterable[sp.Expr | sp.Equality | bool],
    variables: Sequence[sp.Symbol],
) -> tuple[tuple[sp.Expr, ...], object]:
    vars_ = tuple(variables)
    raw: list[sp.Expr] = []
    for item in generators:
        if item is True or item is sp.true or item is sp.S.true:
            continue
        if item is False or item is sp.false or item is sp.S.false:
            raw.append(sp.Integer(1))
            continue
        expr = sp.sympify(item)
        if isinstance(expr, sp.Equality):
            expr = sp.expand(expr.lhs - expr.rhs)
        else:
            expr = sp.expand(expr)
        if expr != 0:
            raw.append(expr)
    if not vars_:
        raise RationalUnivariateError("a nonempty variable list is required")
    if not raw:
        return tuple(), sp.QQ
    try:
        polys, options = sp.parallel_poly_from_expr(raw, *vars_, extension=True)
    except (
        sp.PolynomialError,
        sp.polys.polyerrors.CoercionFailed,
        ValueError,
        TypeError,
        NotImplementedError,
    ) as exc:
        raise RationalUnivariateError(
            "equality-ideal analysis requires exact polynomial generators"
        ) from exc
    domain = require_exact_rur_domain(options["domain"])
    normalized = tuple(sp.Poly(poly.as_expr(), *vars_, domain=domain).as_expr() for poly in polys)
    return normalized, domain


def _leading_exponents(basis: sp.polys.polytools.GroebnerBasis) -> tuple[tuple[int, ...], ...]:
    if str(basis.order) != "grevlex":
        raise RationalUnivariateError("leading-exponent analysis requires a grevlex basis")
    return tuple(_leading_exponent_grevlex(poly) for poly in basis.polys)


def _monomial_ideal_dimension(
    leading_exponents: Sequence[Sequence[int]], variable_count: int
) -> int:
    """Return ``dim k[x]/J`` for a monomial ideal ``J`` exactly.

    The height of a monomial ideal is the minimum size of a set of variables
    meeting the support of every monomial generator.  We compute that minimum
    vertex cover exactly; the small variable counts reached by this solver make
    this preferable to probabilistic dimension tests.
    """

    if not leading_exponents:
        return variable_count
    supports = tuple(
        frozenset(index for index, power in enumerate(exponent) if int(power) > 0)
        for exponent in leading_exponents
    )
    if any(not support for support in supports):
        return -1  # unit ideal / zero ring
    # Remove redundant supersets before the exact minimum hitting-set search.
    reduced = tuple(
        support for support in supports if not any(other < support for other in supports)
    )

    memo: dict[tuple[tuple[int, ...], ...], int] = {}

    def minimum_cover(edges: tuple[frozenset[int], ...]) -> int:
        if not edges:
            return 0
        key = tuple(sorted(tuple(sorted(edge)) for edge in edges))
        cached = memo.get(key)
        if cached is not None:
            return cached
        pivot = min(edges, key=len)
        best = variable_count + 1
        for variable in pivot:
            remaining = tuple(edge for edge in edges if variable not in edge)
            best = min(best, 1 + minimum_cover(remaining))
        memo[key] = best
        return best

    return variable_count - minimum_cover(reduced)


def _rebuild_relation(atom: sp.Expr, residual: sp.Expr) -> sp.Expr:
    zero = sp.Integer(0)
    if isinstance(atom, sp.Equality):
        return sp.Eq(residual, zero)
    if isinstance(atom, sp.Unequality):
        return sp.Ne(residual, zero)
    if isinstance(atom, sp.StrictLessThan):
        return residual < zero
    if isinstance(atom, sp.LessThan):
        return residual <= zero
    if isinstance(atom, sp.StrictGreaterThan):
        return residual > zero
    if isinstance(atom, sp.GreaterThan):
        return residual >= zero
    raise RationalUnivariateError(f"unsupported relational atom: {atom!s}")


def _truth_when_zero(atom: sp.Expr) -> bool:
    if isinstance(atom, sp.Equality):
        return True
    if isinstance(atom, sp.Unequality):
        return False
    if isinstance(atom, (sp.StrictLessThan, sp.StrictGreaterThan)):
        return False
    if isinstance(atom, (sp.LessThan, sp.GreaterThan)):
        return True
    raise RationalUnivariateError(f"unsupported relational atom: {atom!s}")


class EqualityIdealContext:
    """Reusable exact Groebner context for an equality subsystem.

    A grevlex basis is computed once.  Dimension comes from its leading
    monomial ideal; for a zero-dimensional ideal the standard monomials give
    the exact quotient-algebra dimension.  More expensive lex/FGLM bases are
    created lazily and cached per coordinate.
    """

    @classmethod
    def build(
        cls,
        generators: Iterable[sp.Expr | sp.Equality | bool],
        variables: Sequence[sp.Symbol],
    ) -> EqualityIdealContext:
        """Construct a context from exact polynomial equality generators."""

        return cls(generators, variables)

    def __init__(
        self,
        generators: Iterable[sp.Expr | sp.Equality | bool],
        variables: Sequence[sp.Symbol],
    ) -> None:
        """Normalize generators and initialize lazy Groebner and normal-form caches for one equality ideal."""
        self.variables = tuple(variables)
        self.generators, self.domain = _normalize_generators(generators, self.variables)
        self._lex_by_first: dict[sp.Symbol, sp.polys.polytools.GroebnerBasis] = {}
        self._fglm_by_last: dict[sp.Symbol, sp.polys.polytools.GroebnerBasis] = {}
        self._fglm_by_order: dict[tuple[sp.Symbol, ...], sp.polys.polytools.GroebnerBasis] = {}
        self._normal_form_cache: dict[sp.Expr, sp.Expr] = {}
        self._radical_cache: dict[sp.Expr, bool] = {}
        self._common_zero_cache: dict[sp.Expr, bool] = {}
        self._coordinate_cache: dict[sp.Symbol, sp.Poly] = {}
        self._standard_cache: tuple[tuple[int, ...], ...] | None = None
        self._leading_exponents: tuple[tuple[int, ...], ...] = tuple()
        self._stats: dict[str, int] = {
            "groebner_basis_count": 0,
            "fglm_conversion_count": 0,
            "normal_form_reductions": 0,
            "radical_basis_count": 0,
            "common_zero_basis_count": 0,
            "standard_materializations": 0,
        }
        if not self.generators:
            self.groebner_basis = None
            self.inconsistent = False
            self.dimension = len(self.variables)
            self.quotient_dimension = None
            self._standard_cache = tuple()
            return
        self.groebner_basis = compute_groebner_basis(
            self.generators, self.variables, order="grevlex", domain=self.domain
        )
        self._stats["groebner_basis_count"] += 1
        unit = len(self.groebner_basis.polys) == 1 and self.groebner_basis.polys[0].is_one
        self.inconsistent = bool(unit)
        if unit:
            self.dimension = -1
            self.quotient_dimension = 0
            self._standard_cache = tuple()
            return
        leading = _leading_exponents(self.groebner_basis)
        self._leading_exponents = leading
        self.dimension = _monomial_ideal_dimension(leading, len(self.variables))
        if self.dimension == 0:
            self.quotient_dimension = _standard_exponent_count(leading, len(self.variables))
        else:
            self._standard_cache = tuple()
            self.quotient_dimension = None

    @property
    def unit_ideal(self) -> bool:
        return self.inconsistent

    @property
    def zero_dimensional(self) -> bool:
        return self.dimension == 0

    @property
    def standard_exponents(self) -> tuple[tuple[int, ...], ...]:
        """Materialize the standard monomial exponents only when explicitly needed."""

        if not self.zero_dimensional:
            return tuple()
        if self._standard_cache is None:
            self._standard_cache = _standard_exponents(self._leading_exponents, len(self.variables))
            self._stats["standard_materializations"] += 1
        return self._standard_cache

    @property
    def analysis(self) -> EqualityIdealAnalysis:
        basis = (
            tuple(poly.as_expr() for poly in self.groebner_basis.polys)
            if self.groebner_basis is not None
            else tuple()
        )
        return EqualityIdealAnalysis(
            variables=self.variables,
            generators=self.generators,
            groebner_basis=basis,
            dimension=self.dimension,
            quotient_dimension=self.quotient_dimension,
            zero_dimensional=self.zero_dimensional,
            inconsistent=self.inconsistent,
        )

    def diagnostics(self) -> dict[str, int]:
        """Return deterministic counters for expensive equality-ideal operations."""

        return dict(self._stats)

    def normal_form(self, expression: sp.Expr) -> sp.Expr:
        """Return the exact normal form modulo the equality ideal."""

        expr = sp.expand(sp.sympify(expression))
        cached = self._normal_form_cache.get(expr)
        if cached is not None:
            return cached
        if self.groebner_basis is None:
            return expr
        try:
            self._stats["normal_form_reductions"] += 1
            _, remainder = self.groebner_basis.reduce(expr)
        except (sp.PolynomialError, TypeError, ValueError) as exc:
            raise RationalUnivariateError(f"Groebner reduction failed: {exc}") from exc
        reduced = sp.expand(remainder)
        self._normal_form_cache[expr] = reduced
        return reduced

    def _lex_basis_with_first(self, variable: sp.Symbol):
        if variable not in self.variables:
            raise ValueError(f"{variable!s} is not an ideal variable")
        cached = self._lex_by_first.get(variable)
        if cached is not None:
            return cached
        order = (variable, *(v for v in self.variables if v != variable))
        basis = compute_groebner_basis(self.generators, order, order="lex", domain=self.domain)
        self._stats["groebner_basis_count"] += 1
        self._lex_by_first[variable] = basis
        return basis

    def fglm_lex_basis(
        self, variable_order: Sequence[sp.Symbol]
    ) -> sp.polys.polytools.GroebnerBasis:
        """Return a cached exact FGLM lex basis for ``variable_order``."""

        if not self.zero_dimensional:
            raise RationalUnivariateError("FGLM conversion requires a zero-dimensional ideal")
        order = tuple(variable_order)
        if set(order) != set(self.variables) or len(order) != len(self.variables):
            raise ValueError("FGLM variable order must be a permutation of the ideal variables")
        cached = self._fglm_by_order.get(order)
        if cached is not None:
            return cached
        source = compute_groebner_basis(self.generators, order, order="grevlex", domain=self.domain)
        self._stats["groebner_basis_count"] += 1
        try:
            basis = source.fglm(order="lex") if len(order) > 1 else source
        except (NotImplementedError, ValueError, TypeError, sp.PolynomialError) as exc:
            raise RationalUnivariateError(f"FGLM conversion failed: {exc}") from exc
        if len(order) > 1:
            self._stats["fglm_conversion_count"] += 1
        self._fglm_by_order[order] = basis
        return basis

    def triangular_basis(self, lifting_order: Sequence[sp.Symbol]) -> tuple[sp.Expr, ...]:
        """Return a lex basis triangular in the supplied CAD lifting order."""

        order = tuple(reversed(tuple(lifting_order)))
        basis = self.fglm_lex_basis(order)
        return tuple(sp.expand(poly.as_expr()) for poly in basis.polys)

    def _fglm_basis_with_last(self, variable: sp.Symbol):
        if not self.zero_dimensional:
            raise RationalUnivariateError("FGLM conversion requires a zero-dimensional ideal")
        if variable not in self.variables:
            raise ValueError(f"{variable!s} is not an ideal variable")
        cached = self._fglm_by_last.get(variable)
        if cached is not None:
            return cached
        order = (*(v for v in self.variables if v != variable), variable)
        basis = self.fglm_lex_basis(order)
        self._fglm_by_last[variable] = basis
        return basis

    def coordinate_polynomial(self, variable: sp.Symbol) -> sp.Poly:
        """Return the minimal univariate coordinate eliminant found by FGLM."""

        cached = self._coordinate_cache.get(variable)
        if cached is not None:
            return cached
        # Reuse any FGLM basis already built for triangular CAD projection if it
        # already exposes this coordinate eliminant.
        candidates: list[sp.Poly] = []
        for cached_basis in self._fglm_by_order.values():
            for item in cached_basis.polys:
                expr = sp.expand(item.as_expr())
                if expr != 0 and expr.free_symbols <= {variable}:
                    poly = sp.Poly(expr, variable, domain=self.domain)
                    if poly.degree() > 0:
                        candidates.append(poly.monic())
        if candidates:
            result = min(candidates, key=lambda poly: (poly.degree(), len(poly.terms())))
            self._coordinate_cache[variable] = result
            return result
        basis = self._fglm_basis_with_last(variable)
        candidates = []
        for item in basis.polys:
            expr = sp.expand(item.as_expr())
            if expr != 0 and expr.free_symbols <= {variable}:
                poly = sp.Poly(expr, variable, domain=self.domain)
                if poly.degree() > 0:
                    candidates.append(poly.monic())
        if not candidates:
            raise RationalUnivariateError(
                f"FGLM basis did not expose a coordinate polynomial for {variable!s}"
            )
        result = min(candidates, key=lambda poly: (poly.degree(), len(poly.terms())))
        self._coordinate_cache[variable] = result
        return result

    def coordinate_polynomials(self) -> dict[sp.Symbol, sp.Poly]:
        """Return exact FGLM coordinate eliminants for every variable."""

        return {variable: self.coordinate_polynomial(variable) for variable in self.variables}

    def fglm_coordinate_polynomials(self) -> tuple[sp.Expr, ...]:
        """Return FGLM coordinate eliminants in the context variable order."""

        polynomials = self.coordinate_polynomials()
        return tuple(polynomials[variable].as_expr() for variable in self.variables)

    def eliminate_linear_variable(self, variable: sp.Symbol) -> LinearVariableElimination | None:
        """Find a global relation ``variable = p(other variables)`` in the ideal.

        A lexicographic basis with ``variable`` first guarantees that such a
        polynomial relation, when present in reduced form with constant leading
        coefficient, is exposed without dividing by a possibly vanishing
        polynomial.  The returned replacement is therefore globally valid on
        the equality variety.
        """

        if self.inconsistent:
            return None
        basis = self._lex_basis_with_first(variable)
        candidates: list[tuple[int, sp.Expr, sp.Expr]] = []
        for item in basis.polys:
            expr = sp.expand(item.as_expr())
            poly = sp.Poly(expr, variable, domain="EX")
            if poly.degree() != 1:
                continue
            coefficient = sp.expand(poly.nth(1))
            constant_part = sp.expand(poly.nth(0))
            if coefficient == 0 or coefficient.free_symbols.intersection(self.variables):
                continue
            replacement = sp.cancel(-constant_part / coefficient)
            if variable in replacement.free_symbols:
                continue
            certificate = sp.expand(variable - replacement)
            try:
                _, remainder = basis.reduce(certificate)
            except (sp.PolynomialError, TypeError, ValueError):
                continue
            if sp.expand(remainder) != 0:
                continue
            score = int(sp.count_ops(replacement))
            candidates.append((score, replacement, certificate))
        if not candidates:
            return None
        _, replacement, certificate = min(candidates, key=lambda item: (item[0], sp.sstr(item[1])))
        remaining = tuple(
            sp.expand(item.as_expr())
            for item in basis.polys
            if variable not in item.as_expr().free_symbols
        )
        return LinearVariableElimination(variable, replacement, remaining, certificate)

    def linear_eliminations(self) -> tuple[LinearVariableElimination, ...]:
        """Return all globally certified single-variable polynomial replacements."""

        eliminations: list[LinearVariableElimination] = []
        for variable in self.variables:
            elimination = self.eliminate_linear_variable(variable)
            if elimination is not None:
                eliminations.append(elimination)
        return tuple(eliminations)

    def radical_contains(self, expression: sp.Expr) -> bool:
        """Return whether ``expression`` lies in the radical of the ideal."""

        return self.in_radical(expression)

    def simplify_relations(self, formula: sp.Expr) -> tuple[sp.Expr, dict[str, object]]:
        """Reduce residual relations modulo the ideal and collapse ideal consequences."""

        atoms = tuple(formula.args) if isinstance(formula, sp.And) else (formula,)
        simplified: list[sp.Expr] = []
        normal_form_changes = 0
        radical_collapses = 0
        for atom in atoms:
            if isinstance(atom, sp.Equality):
                continue
            if not getattr(atom, "is_Relational", False):
                simplified.append(atom)
                continue
            residual = sp.expand(atom.lhs - atom.rhs)
            remainder = self.normal_form(residual)
            if sp.expand(remainder - residual) != 0:
                normal_form_changes += 1
            if remainder == 0 or self.in_radical(remainder):
                radical_collapses += 1
                if not _truth_when_zero(atom):
                    return sp.false, {
                        "normal_form_changes": normal_form_changes,
                        "radical_collapses": radical_collapses,
                    }
                continue
            rebuilt = self._simplify_reduced_relation(atom, remainder)
            if rebuilt is sp.false or rebuilt == sp.false:
                return sp.false, {
                    "normal_form_changes": normal_form_changes,
                    "radical_collapses": radical_collapses,
                }
            if rebuilt is not sp.true and rebuilt != sp.true:
                simplified.append(sp.sympify(rebuilt))
        equations = (
            [sp.Eq(poly.as_expr(), 0) for poly in self.groebner_basis.polys]
            if self.groebner_basis is not None
            else []
        )
        result = sp.And(*equations, *simplified)
        return result, {
            "normal_form_changes": normal_form_changes,
            "radical_collapses": radical_collapses,
        }

    def filter_zero_dimensional(
        self, constraints: sp.Expr | bool = sp.true
    ) -> ZeroDimensionalFilterResult:
        """Enumerate the finite real variety and filter residual constraints exactly."""

        if not self.zero_dimensional:
            raise ValueError("zero-dimensional filtering requires a finite equality variety")
        from .rational_univariate import solve_and_filter_zero_dimensional_system_with_rur

        simplified = self.simplify_constraints(constraints)
        filtered = solve_and_filter_zero_dimensional_system_with_rur(
            self.generators, self.variables, simplified, real=True
        )
        return ZeroDimensionalFilterResult(
            variables=self.variables,
            points=filtered.points,
            quotient_dimension=self.quotient_dimension,
            coordinate_polynomials=self.fglm_coordinate_polynomials(),
            representation=filtered.representation,
        )

    def in_radical(self, expression: sp.Expr) -> bool:
        """Return whether ``expression`` belongs to the radical of this ideal."""

        expr = self.normal_form(sp.expand(sp.sympify(expression)))
        if expr == 0 or self.inconsistent:
            return True
        cached = self._radical_cache.get(expr)
        if cached is not None:
            return cached
        parameter = fresh_dummy("radical_membership")
        variables = (parameter, *self.variables)
        try:
            self._stats["radical_basis_count"] += 1
            basis = compute_groebner_basis(
                (*self.generators, 1 - parameter * expr),
                variables,
                order="grevlex",
                extension=True,
            )
        except (
            sp.PolynomialError,
            sp.polys.polyerrors.CoercionFailed,
            TypeError,
            ValueError,
            NotImplementedError,
        ) as exc:
            raise RationalUnivariateError(f"radical membership test failed: {exc}") from exc
        result = len(basis.polys) == 1 and basis.polys[0].is_one
        self._radical_cache[expr] = result
        return result

    def has_common_zero_with(self, expression: sp.Expr) -> bool:
        """Return whether the ideal together with ``expression == 0`` has a complex zero."""

        if self.inconsistent:
            return False
        expr = self.normal_form(sp.expand(sp.sympify(expression)))
        if expr == 0:
            return True
        cached = self._common_zero_cache.get(expr)
        if cached is not None:
            return cached
        try:
            self._stats["common_zero_basis_count"] += 1
            basis = compute_groebner_basis(
                (*self.generators, expr), self.variables, order="grevlex", extension=True
            )
        except (
            sp.PolynomialError,
            sp.polys.polyerrors.CoercionFailed,
            TypeError,
            ValueError,
            NotImplementedError,
        ) as exc:
            raise RationalUnivariateError(f"common-zero test failed: {exc}") from exc
        result = not (len(basis.polys) == 1 and basis.polys[0].is_one)
        self._common_zero_cache[expr] = result
        return result

    def _simplify_reduced_relation(self, atom: sp.Expr, remainder: sp.Expr) -> sp.Expr | bool:
        """Simplify a relation from a normal form already computed by this context."""

        if not remainder.free_symbols.intersection(self.variables):
            return sp.simplify(_rebuild_relation(atom, remainder))
        if isinstance(atom, sp.Equality) and not self.has_common_zero_with(remainder):
            return sp.false
        if isinstance(atom, sp.Unequality) and not self.has_common_zero_with(remainder):
            return sp.true
        return _rebuild_relation(atom, remainder)

    def simplify_relation(self, atom: sp.Expr) -> sp.Expr | bool:
        """Simplify one polynomial relation using exact ideal consequences."""

        if not getattr(atom, "is_Relational", False):
            raise RationalUnivariateError(f"expected a relational atom, got {atom!s}")
        residual = sp.expand(atom.lhs - atom.rhs)
        remainder = self.normal_form(residual)
        if remainder == 0 or self.in_radical(remainder):
            return sp.true if _truth_when_zero(atom) else sp.false
        return self._simplify_reduced_relation(atom, remainder)

    def simplify_constraints(self, formula: sp.Expr | bool) -> sp.Expr | bool:
        """Simplify a Boolean constraint formula on the equality variety."""

        if formula is True or formula is sp.true or formula is sp.S.true:
            return sp.true
        if formula is False or formula is sp.false or formula is sp.S.false:
            return sp.false
        expr = sp.sympify(formula)
        if isinstance(expr, sp.And):
            return sp.And(*(self.simplify_constraints(arg) for arg in expr.args))
        if isinstance(expr, sp.Or):
            return sp.Or(*(self.simplify_constraints(arg) for arg in expr.args))
        if isinstance(expr, sp.Not):
            return sp.Not(self.simplify_constraints(expr.args[0]))
        return self.simplify_relation(expr)


def analyze_equality_ideal(
    generators: Iterable[sp.Expr | sp.Equality | bool], variables: Sequence[sp.Symbol]
) -> EqualityIdealAnalysis:
    """Return exact ideal dimension and quotient-dimension information."""

    return EqualityIdealContext(generators, variables).analysis


__all__ = [
    "EqualityIdealAnalysis",
    "EqualityIdealContext",
    "LinearVariableElimination",
    "ZeroDimensionalFilterResult",
    "analyze_equality_ideal",
]
