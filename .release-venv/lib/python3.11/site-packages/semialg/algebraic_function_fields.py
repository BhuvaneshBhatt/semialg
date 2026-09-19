"""Exact monogenic algebraic function fields used by triangular decomposition.

The implementation is small and certificate-oriented.  It supports
rational-function base fields and finite monogenic towers ``K[a]/(m(a))``.
Univariate factorization uses a Trager-style norm descent when possible; callers
can fall back to quotient-ring splitting when the norm route is inconclusive.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass

import sympy as sp


class FunctionFieldError(Exception):
    pass


@dataclass(frozen=True)
class CertifiedFactorization:
    """Exact factorization together with independent reconstruction status."""

    unit: object
    factors: tuple[tuple[tuple[object, ...], int], ...]
    reconstruction_verified: bool
    method: str

    @property
    def irreducible(self) -> bool:
        return self.reconstruction_verified and len(self.factors) == 1 and self.factors[0][1] == 1


@dataclass(frozen=True)
class PrimitiveElementCompressionStep:
    """Replayable certificate for one adjacent primitive-element merge."""

    lower_generator: sp.Symbol
    upper_generator: sp.Symbol
    parameter_symbol: sp.Symbol
    parameter: int
    bad_parameter_guard: sp.Expr
    degree_before: int
    degree_after: int
    primitive_expression: sp.Expr


@dataclass(frozen=True)
class PrimitiveElementCompression:
    """Exact compression of a finite monogenic tower to one primitive element."""

    field: object
    generator_images: tuple[tuple[sp.Symbol, object], ...]
    complete: bool
    parameter: int | None
    method: str
    steps: tuple[PrimitiveElementCompressionStep, ...] = ()
    primitive_expression: sp.Expr | None = None
    source_field: object | None = None

    def image(self, generator: sp.Symbol):
        for symbol, value in self.generator_images:
            if symbol == generator:
                return value
        raise KeyError(generator)

    def to_compressed(self, value):
        """Map an element of the original tower into the compressed field."""
        if self.source_field is None:
            raise FunctionFieldError("compression does not record its source field")
        expr = field_element_to_expr(value, self.source_field)
        return evaluate_expression(expr, self.field, dict(self.generator_images))

    def to_source(self, value):
        """Map an element of the compressed field back to the original tower."""
        if self.source_field is None or self.primitive_expression is None:
            raise FunctionFieldError("compression does not record an inverse primitive map")
        expr = field_element_to_expr(value, self.field)
        expr = sp.cancel(expr.subs(self.field.generator, self.primitive_expression))
        return tower_evaluate(expr, self.source_field)


@dataclass(frozen=True)
class RationalFunctionField:
    parameters: tuple[sp.Symbol, ...]

    @property
    def zero(self):
        return sp.Integer(0)

    @property
    def one(self):
        return sp.Integer(1)

    def convert(self, value):
        return self.normalize(value)

    def normalize(self, value):
        domain = sp.QQ.frac_field(*self.parameters) if self.parameters else sp.QQ
        # Domain conversion already performs exact rational-function
        # cancellation.  Running SymPy ``cancel`` before and after it builds
        # large transient expression graphs in long algebraic workloads.
        return domain.to_sympy(domain.convert(sp.sympify(value)))

    def is_zero(self, value) -> bool:
        return self.normalize(value) == 0

    def add(self, a, b):
        return self.normalize(a + b)

    def sub(self, a, b):
        return self.normalize(a - b)

    def mul(self, a, b):
        return self.normalize(a * b)

    def neg(self, a):
        return self.normalize(-a)

    def pow(self, a, n):
        return self.normalize(a**n)

    def inv(self, value):
        value = self.normalize(value)
        if value == 0:
            raise ZeroDivisionError
        return self.normalize(1 / value)

    def squarefree_decomposition(self, coeffs: Sequence[object]):
        return _squarefree_decomposition_with_unit(coeffs, self)

    def factor_univariate(self, coeffs: Sequence[object]):
        normalized = tuple(self.normalize(c) for c in coeffs)
        # Exact QQ leaves use the modular integer factorizer followed by exact
        # rational reconstruction verification.  Parameter-dependent leaves
        # remain on SymPy's exact rational-function factorization path.
        if all(not sp.sympify(c).free_symbols.intersection(self.parameters) for c in normalized):
            from .algebraic.modular import factor_univariate_qq_modular

            modular = factor_univariate_qq_modular(normalized)
            return self.normalize(modular.unit), tuple(
                (tuple(self.normalize(c) for c in factor), multiplicity)
                for factor, multiplicity in modular.factors
            )
        z = sp.Dummy("ff_z")
        expr = sum(self.normalize(c) * z**i for i, c in enumerate(normalized))
        domain = sp.QQ.frac_field(*self.parameters) if self.parameters else sp.QQ
        poly = sp.Poly(expr, z, domain=domain)
        unit, factors = poly.factor_list()
        return self.normalize(unit), tuple(
            (tuple(self.normalize(c) for c in reversed(f.all_coeffs())), m) for f, m in factors
        )

    def factor_univariate_certified(self, coeffs: Sequence[object]):
        return certified_factor_univariate(self, coeffs)


@dataclass(frozen=True)
class MonogenicElement:
    field: MonogenicFunctionField
    coeffs: tuple[object, ...]

    def __add__(self, other):
        return self.field.add(self, other)

    def __radd__(self, other):
        return self.field.add(other, self)

    def __sub__(self, other):
        return self.field.sub(self, other)

    def __rsub__(self, other):
        return self.field.sub(other, self)

    def __neg__(self):
        return self.field.neg(self)

    def __mul__(self, other):
        return self.field.mul(self, other)

    def __rmul__(self, other):
        return self.field.mul(other, self)

    def __truediv__(self, other):
        return self.field.mul(self, self.field.inv(other))

    def __rtruediv__(self, other):
        return self.field.mul(other, self.field.inv(self))

    def __pow__(self, n: int):
        return self.field.pow(self, n)


@dataclass(frozen=True)
class MonogenicFunctionField:
    base: object
    generator: sp.Symbol
    modulus: tuple[object, ...]  # low -> high, monic

    def __post_init__(self):
        if len(self.modulus) < 2:
            raise ValueError("modulus must have positive degree")
        if not self.base.is_zero(self.modulus[-1] - self.base.one):
            raise ValueError("modulus must be monic")

    @property
    def degree(self):
        return len(self.modulus) - 1

    @property
    def zero(self):
        return MonogenicElement(self, (self.base.zero,) * self.degree)

    @property
    def one(self):
        return self.convert(self.base.one)

    @property
    def alpha(self):
        coeffs = [self.base.zero] * self.degree
        if self.degree == 1:
            coeffs[0] = (
                self.base.neg(self.modulus[0]) if hasattr(self.base, "neg") else -self.modulus[0]
            )
        else:
            coeffs[1] = self.base.one
        return MonogenicElement(self, tuple(coeffs))

    def _base_add(self, a, b):
        return (
            self.base.normalize(a + b)
            if isinstance(self.base, RationalFunctionField)
            else self.base.add(a, b)
        )

    def _base_sub(self, a, b):
        return (
            self.base.normalize(a - b)
            if isinstance(self.base, RationalFunctionField)
            else self.base.sub(a, b)
        )

    def _base_mul(self, a, b):
        return (
            self.base.normalize(a * b)
            if isinstance(self.base, RationalFunctionField)
            else self.base.mul(a, b)
        )

    def _base_neg(self, a):
        return (
            self.base.normalize(-a)
            if isinstance(self.base, RationalFunctionField)
            else self.base.neg(a)
        )

    def convert(self, value):
        if isinstance(value, MonogenicElement):
            if value.field == self:
                return value
            if value.field == self.base:
                base_value = value
                return MonogenicElement(self, (base_value,) + (self.base.zero,) * (self.degree - 1))
            raise TypeError("element from different field")
        base_value = (
            self.base.normalize(value)
            if isinstance(self.base, RationalFunctionField)
            else self.base.convert(value)
        )
        return MonogenicElement(self, (base_value,) + (self.base.zero,) * (self.degree - 1))

    def normalize(self, value):
        return self.convert(value)

    def is_zero(self, value):
        value = self.convert(value)
        return all(self.base.is_zero(c) for c in value.coeffs)

    def add(self, a, b):
        a = self.convert(a)
        b = self.convert(b)
        return MonogenicElement(
            self, tuple(self._base_add(x, y) for x, y in zip(a.coeffs, b.coeffs, strict=True))
        )

    def neg(self, a):
        a = self.convert(a)
        return MonogenicElement(self, tuple(self._base_neg(x) for x in a.coeffs))

    def sub(self, a, b):
        return self.add(a, self.neg(b))

    def mul(self, a, b):
        a = self.convert(a)
        b = self.convert(b)
        d = self.degree
        tmp = [self.base.zero] * (2 * d - 1)
        for i, x in enumerate(a.coeffs):
            for j, y in enumerate(b.coeffs):
                tmp[i + j] = self._base_add(tmp[i + j], self._base_mul(x, y))
        for k in range(len(tmp) - 1, d - 1, -1):
            lead = tmp[k]
            if self.base.is_zero(lead):
                continue
            for j in range(d):
                tmp[k - d + j] = self._base_sub(
                    tmp[k - d + j], self._base_mul(lead, self.modulus[j])
                )
        return MonogenicElement(self, tuple(tmp[:d]))

    def pow(self, a, n):
        if n < 0:
            return self.pow(self.inv(a), -n)
        out = self.one
        base = self.convert(a)
        while n:
            if n & 1:
                out = self.mul(out, base)
            base = self.mul(base, base)
            n //= 2
        return out

    def inv(self, a):
        a = self.convert(a)
        if self.is_zero(a):
            raise ZeroDivisionError
        # Extended Euclid in base[t].
        p = list(a.coeffs)
        m = list(self.modulus)
        r0, r1 = m, p
        t0, t1 = [self.base.zero], [self.base.one]
        while _poly_degree(r1, self.base) >= 0:
            q, r2 = _poly_divmod(r0, r1, self.base)
            t2 = _poly_sub(t0, _poly_mul(q, t1, self.base), self.base)
            r0, r1, t0, t1 = r1, r2, t1, t2
        if _poly_degree(r0, self.base) != 0:
            raise FunctionFieldError("noninvertible residue")
        invc = self.base.inv(r0[0])
        coeffs = [self._base_mul(c, invc) for c in t0]
        return self.mul(
            MonogenicElement(self, tuple(coeffs + [self.base.zero] * (self.degree - len(coeffs)))),
            self.one,
        )

    def squarefree_decomposition(self, coeffs: Sequence[object]):
        """Return exact monic squarefree factors and their multiplicities."""
        return _squarefree_decomposition_with_unit(coeffs, self)

    def factor_univariate(self, coeffs: Sequence[object], *, auto_compress: bool = True):
        """Factor an arbitrary nonzero univariate polynomial exactly.

        Characteristic-zero squarefree decomposition is performed directly in
        the algebraic-function-field tower.  Each squarefree component is then
        factored by exact Trager norm descent.

        The returned unit is the original leading coefficient and every factor
        is monic; multiplicities come from the squarefree decomposition.
        """
        coeffs = _poly_trim([self.convert(c) for c in coeffs], self)
        degree = _poly_degree(coeffs, self)
        if degree < 0:
            raise FunctionFieldError("cannot factor the zero polynomial")
        if degree == 0:
            return coeffs[0], ()

        if auto_compress and should_compress_tower(self) and tower_degree(self) * degree <= 6:
            compression = maybe_compress_primitive_element(self)
            if compression is not None:
                try:
                    compressed_coeffs = [compression.to_compressed(c) for c in coeffs]
                    compressed_unit, compressed_factors = compression.field.factor_univariate(
                        compressed_coeffs, auto_compress=False
                    )
                    source_unit = compression.to_source(compressed_unit)
                    source_factors = tuple(
                        (tuple(compression.to_source(c) for c in factor), multiplicity)
                        for factor, multiplicity in compressed_factors
                    )
                    reconstructed = _poly_reconstruct_factorization(
                        source_unit, source_factors, self
                    )
                    if len(reconstructed) == len(coeffs) and all(
                        self.is_zero(a - b) for a, b in zip(reconstructed, coeffs, strict=True)
                    ):
                        return source_unit, source_factors
                except (
                    ArithmeticError,
                    ValueError,
                    TypeError,
                    ZeroDivisionError,
                    FunctionFieldError,
                    sp.PolynomialError,
                    sp.polys.polyerrors.CoercionFailed,
                ):
                    pass

        unit, squarefree_parts = self.squarefree_decomposition(coeffs)
        factors = []
        for squarefree_factor, multiplicity in squarefree_parts:
            _unit, irreducibles = self._factor_squarefree_univariate(
                squarefree_factor,
            )
            for factor, inner_multiplicity in irreducibles:
                factors.append((factor, multiplicity * inner_multiplicity))
        return unit, tuple(factors)

    def factor_univariate_certified(self, coeffs: Sequence[object]):
        return certified_factor_univariate(self, coeffs)

    def _factor_squarefree_univariate(self, coeffs: Sequence[object]):
        """Factor a monic squarefree polynomial by exact Trager norm descent."""
        coeffs = _poly_trim([self.convert(c) for c in coeffs], self)
        if _poly_degree(coeffs, self) <= 1:
            return self.one, ((tuple(_poly_monic(coeffs, self)), 1),)
        if not _poly_is_squarefree(coeffs, self):
            raise FunctionFieldError("internal Trager input is not squarefree")

        shift, norm = _first_good_trager_shift(coeffs, self)
        try:
            _unit, lower_factors = self.base.factor_univariate(norm)
        except (
            ArithmeticError,
            ValueError,
            TypeError,
            ZeroDivisionError,
            FunctionFieldError,
            sp.PolynomialError,
            sp.polys.polyerrors.CoercionFailed,
        ) as exc:
            raise FunctionFieldError("norm factorization inconclusive") from exc
        if not lower_factors or any(m != 1 for _g, m in lower_factors):
            raise FunctionFieldError("certified-good norm did not factor squarefreely")

        pieces = []
        remaining = coeffs
        for lower_factor, _m in lower_factors:
            lifted = [self.convert(c) for c in lower_factor]
            # Undo f(z - s*alpha): recover with g(z + s*alpha).
            lifted = _poly_shift_alpha(lifted, self, -shift)
            g = _poly_gcd(remaining, lifted, self)
            dg = _poly_degree(g, self)
            if 0 < dg < _poly_degree(remaining, self):
                pieces.append((tuple(_poly_monic(g, self)), 1))
                q, r = _poly_divmod(remaining, g, self)
                if _poly_degree(r, self) >= 0:
                    raise FunctionFieldError("inexact factor recovery")
                remaining = q
        if pieces:
            if _poly_degree(remaining, self) > 0:
                pieces.append((tuple(_poly_monic(remaining, self)), 1))
            return self.one, tuple(pieces)

        # For a certified-good Trager shift, an irreducible norm certifies that
        # the original polynomial is irreducible over the extension field.
        if len(lower_factors) == 1 and _poly_degree(lower_factors[0][0], self.base) == _poly_degree(
            norm, self.base
        ):
            return self.one, ((tuple(_poly_monic(coeffs, self)), 1),)
        raise FunctionFieldError("certified norm factors did not recover factors")


def tower_depth(field) -> int:
    """Return the number of monogenic layers above the rational base."""
    _base, levels = _tower_levels(field)
    return len(levels)


def tower_degree(field) -> int:
    """Return the total finite extension degree over the rational base."""
    _base, levels = _tower_levels(field)
    degree = 1
    for level in levels:
        degree *= level.degree
    return degree


def should_compress_tower(
    field, *, min_depth: int = 3, min_degree: int = 8, max_degree: int = 64
) -> bool:
    """Deterministic expression-growth heuristic for automatic compression."""
    if not isinstance(field, MonogenicFunctionField):
        return False
    depth = tower_depth(field)
    degree = tower_degree(field)
    return degree <= max_degree and (depth >= min_depth or (depth >= 2 and degree >= min_degree))


def maybe_compress_primitive_element(field, **heuristic_options):
    """Compress when the deterministic tower-growth heuristic recommends it."""
    if not should_compress_tower(field, **heuristic_options):
        return None
    result = compress_primitive_element(field)
    return result if verify_primitive_element_compression(result) else None


def element_to_expr(value, generator: sp.Symbol):
    if not isinstance(value, MonogenicElement):
        return sp.sympify(value)
    return sp.cancel(
        sum(element_to_expr(c, generator) * generator**i for i, c in enumerate(value.coeffs))
    )


def evaluate_expression(expression: sp.Expr, field, generators: dict[sp.Symbol, object]):
    """Evaluate an exact rational expression in a function-field tower."""
    expression = sp.cancel(sp.sympify(expression))
    if expression.is_Rational:
        return (
            field.convert(expression) if hasattr(field, "convert") else field.normalize(expression)
        )
    if expression.is_Symbol:
        if expression in generators:
            return generators[expression]
        return (
            field.convert(expression) if hasattr(field, "convert") else field.normalize(expression)
        )
    if expression.is_Add:
        out = field.zero
        for arg in expression.args:
            out = out + evaluate_expression(arg, field, generators)
        return out
    if expression.is_Mul:
        out = field.one
        for arg in expression.args:
            out = out * evaluate_expression(arg, field, generators)
        return out
    if expression.is_Pow and expression.exp.is_Integer:
        base = evaluate_expression(expression.base, field, generators)
        n = int(expression.exp)
        if n >= 0:
            return base**n if hasattr(base, "__pow__") else field.normalize(base**n)
        inv = field.inv(base)
        return inv ** (-n) if hasattr(inv, "__pow__") else field.normalize(inv ** (-n))
    raise FunctionFieldError(f"unsupported field expression: {expression}")


def field_element_to_expr(value, field) -> sp.Expr:
    """Convert a tower element back to a SymPy representative."""
    if isinstance(field, RationalFunctionField):
        return sp.cancel(value)
    value = field.convert(value)
    return sp.cancel(
        sum(
            field_element_to_expr(c, field.base) * field.generator**i
            for i, c in enumerate(value.coeffs)
        )
    )


def tower_evaluate(expression: sp.Expr, field):
    """Evaluate using generator names encoded by a nested field tower."""
    expression = sp.cancel(sp.sympify(expression))
    if isinstance(field, RationalFunctionField):
        return field.normalize(expression)
    if expression == field.generator:
        return field.alpha
    if field.generator not in expression.free_symbols:
        return field.convert(tower_evaluate(expression, field.base))
    if expression.is_Add:
        out = field.zero
        for arg in expression.args:
            out = out + tower_evaluate(arg, field)
        return out
    if expression.is_Mul:
        out = field.one
        for arg in expression.args:
            out = out * tower_evaluate(arg, field)
        return out
    if expression.is_Pow and expression.exp.is_Integer:
        base = tower_evaluate(expression.base, field)
        return field.pow(base, int(expression.exp))
    raise FunctionFieldError(f"unsupported tower expression: {expression}")


__all__ = [
    "CertifiedFactorization",
    "FunctionFieldError",
    "MonogenicElement",
    "MonogenicFunctionField",
    "PrimitiveElementCompression",
    "PrimitiveElementCompressionStep",
    "RationalFunctionField",
    "certified_factor_univariate",
    "compress_primitive_element",
    "maybe_compress_primitive_element",
    "should_compress_tower",
    "tower_degree",
    "tower_depth",
    "verify_primitive_element_compression",
    "evaluate_expression",
    "field_element_to_expr",
    "tower_evaluate",
]

# Factorization machinery is implemented out of the representation module.
# Primitive-element realization/certification is kept separate from field arithmetic.
from ._function_field_compression import (  # noqa: E402, F401  # noqa: E402, F401
    _first_integer_off_guard,
    _merge_two_monogenic_fields,
    _primitive_bad_parameter_guard,
    _tower_levels,
    compress_primitive_element,
    verify_primitive_element_compression,
)

# Private implementation imports shared with sibling algebraic modules.
from ._function_field_factorization import (  # noqa: E402, F401  # noqa: E402, F401
    _det_poly,
    _first_good_trager_shift,
    _norm_polynomial,
    _poly_add,
    _poly_degree,
    _poly_derivative,
    _poly_divmod,
    _poly_gcd,
    _poly_is_squarefree,
    _poly_monic,
    _poly_mul,
    _poly_reconstruct_factorization,
    _poly_shift_alpha,
    _poly_squarefree_decomposition,
    _poly_sub,
    _poly_trim,
    _squarefree_decomposition_with_unit,
    _trager_shift_is_good,
    certified_factor_univariate,
)
