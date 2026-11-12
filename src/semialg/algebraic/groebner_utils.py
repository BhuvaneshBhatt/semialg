"""Small shared construction policy for exact Groebner bases."""

from __future__ import annotations

from collections.abc import Iterable, Sequence
from contextvars import ContextVar
from functools import wraps
from typing import Any

import sympy as sp

_construction_backend = ContextVar("semialg_groebner_backend", default="native")
_construction_limits = ContextVar("semialg_groebner_limits", default=None)


def _backend_policy(function):
    """Keep recursive adapter calls on the explicitly selected construction backend.

    A context variable isolates concurrent requests. The wrapper supplies an
    inherited choice only when the nested caller omits its keyword; explicit
    selections always win. Basis caches include the resolved backend.
    """

    @wraps(function)
    def scoped(*args, **kwargs):
        from algroots.basis import _validate_backend

        backend = kwargs.setdefault("groebner_backend", _construction_backend.get())
        inherited_limits = (
            _construction_limits.get() if backend == _construction_backend.get() else None
        )
        limits = kwargs.setdefault("groebner_limits", inherited_limits)
        _validate_backend(backend, limits)
        token = _construction_backend.set(backend)
        limits_token = _construction_limits.set(limits)
        try:
            return function(*args, **kwargs)
        finally:
            _construction_limits.reset(limits_token)
            _construction_backend.reset(token)

    return scoped


def _selected_backend(backend):
    return _construction_backend.get() if backend is None else backend


def _selected_limits(limits):
    return _construction_limits.get() if limits is None else limits


def polynomial_coordinates(generators, variables):
    """Shield bound root-polynomial symbols from coordinate-ring generators.

    SymPy's polynomial parser inspects ``has(x)``, including the bound x in
    CRootOf. Renaming only free coordinates avoids mistaking that bound symbol
    for a coordinate dependence, while retaining the selected embedding.
    """
    expressions, variables = tuple(map(sp.sympify, generators)), tuple(variables)
    roots = set().union(
        *(sp.sympify(expr).atoms(sp.CRootOf, sp.AlgebraicNumber) for expr in expressions)
    )
    if any(root.has(*variables) for root in roots):
        coordinates = tuple(sp.Dummy(variable.name) for variable in variables)
        replacements = dict(zip(variables, coordinates, strict=True))
        return tuple(
            expr.subs(replacements, simultaneous=True) for expr in expressions
        ), coordinates
    return expressions, variables


def algebraic_coefficient_field(generators, variables):
    """Infer one exact number field, preserving coefficient embeddings.

    Unlisted symbols are parameters, not algebraic constants. Rejecting them
    avoids treating a special parameter stratum as a generic function field.
    Rational input retains the fast QQ arithmetic path.
    """
    expressions, coordinates = polynomial_coordinates(generators, variables)
    _polynomials, options = sp.parallel_poly_from_expr(expressions, *coordinates, extension=True)
    domain = options.domain
    if not (domain.is_ZZ or domain.is_QQ or domain.is_AlgebraicField):
        raise sp.PolynomialError("coefficients must be exact algebraic numbers")
    return domain.get_field()


def _prefer_modular_fraction_field(
    generators: Sequence[sp.Expr], variables: Sequence[sp.Symbol], domain: Any
) -> bool:
    """Estimate when QQ(U) modular proposal overhead is likely worthwhile."""
    max_degree = 0
    coefficient_complexity = 0
    try:
        for generator in generators:
            poly = sp.Poly(generator, *variables, domain=domain)
            max_degree = max(max_degree, int(poly.total_degree()))
            for _monom, coefficient in poly.terms():
                expr = sp.cancel(domain.to_sympy(coefficient))
                numerator, denominator = sp.fraction(sp.together(expr))
                for part in (numerator, denominator):
                    if domain.symbols:
                        coefficient_complexity += max(
                            1, len(sp.Poly(part, *domain.symbols, domain=sp.QQ).terms())
                        )
    except (
        AttributeError,
        ValueError,
        TypeError,
        sp.PolynomialError,
        sp.polys.polyerrors.CoercionFailed,
    ):
        return False
    structural = len(generators) * max(1, len(variables)) * max(1, max_degree)
    return structural >= 12 or coefficient_complexity >= 24


def compute_groebner_basis(
    generators: Iterable[sp.Expr],
    variables: Sequence[sp.Symbol],
    *,
    order: str = "grevlex",
    domain: Any | None = None,
    extension: bool | None = None,
    modulus: int | None = None,
    modular: bool | None = None,
    backend=None,
    limits=None,
) -> sp.polys.polytools.GroebnerBasis:
    """Construct a Groebner basis with one consistent zero-variable guard.

    An omitted backend inherits the enclosing adapter's choice, defaulting to
    native F4. ``backend="sympy"`` selects portable construction;
    ``modular`` controls only its legacy modular proposals. Domain, extension
    and modulus options retain their coefficient-field meanings on both paths.
    """

    backend = _selected_backend(backend)
    if backend == _construction_backend.get():
        limits = _selected_limits(limits)
    vars_tuple = tuple(variables)
    if not vars_tuple:
        raise sp.PolynomialError("Groebner basis construction requires variables")
    generators_tuple = tuple(generators)
    if backend != "sympy" or limits is not None:
        if modulus is not None:
            if domain is not None or extension is not None:
                raise ValueError("modulus cannot be combined with domain or extension")
            domain = sp.GF(modulus)
        elif extension is not None:
            if domain is not None:
                raise ValueError("extension cannot be combined with domain")
            if extension is not True:
                _, options = sp.parallel_poly_from_expr(
                    generators_tuple, *vars_tuple, extension=extension
                )
                domain = options.domain.get_field()
        from algroots import groebner_basis

        return groebner_basis(
            generators_tuple, vars_tuple, domain=domain, order=order, backend=backend, limits=limits
        )

    if modular is not False and modulus is None and extension is None:
        qq_eligible = domain is None or domain == sp.QQ
        if qq_eligible:
            try:
                from .modular import modular_groebner_basis_qq

                result = modular_groebner_basis_qq(generators_tuple, vars_tuple, order=order)
            except (
                ArithmeticError,
                sp.PolynomialError,
                sp.polys.polyerrors.CoercionFailed,
                ValueError,
                TypeError,
            ):
                result = None
            if result is not None:
                return sp.groebner(result.basis, *vars_tuple, order=order, domain=sp.QQ)

    if modular is not False and modulus is None and extension is None and domain is not None:
        try:
            from sympy.polys.domains.fractionfield import FractionField

            is_qq_fraction = isinstance(domain, FractionField) and domain.domain == sp.QQ
        except (ImportError, AttributeError):
            is_qq_fraction = False
        if (
            is_qq_fraction
            and tuple(domain.symbols)
            and _prefer_modular_fraction_field(generators_tuple, vars_tuple, domain)
        ):
            try:
                from .modular import modular_groebner_basis_fraction_field

                result = modular_groebner_basis_fraction_field(
                    generators_tuple, vars_tuple, tuple(domain.symbols), order=order
                )
            except (
                ArithmeticError,
                sp.PolynomialError,
                sp.polys.polyerrors.CoercionFailed,
                ValueError,
                TypeError,
                ZeroDivisionError,
            ):
                result = None
            if result is not None:
                return sp.groebner(result.basis, *vars_tuple, order=order, domain=domain)

    options: dict[str, object] = {"order": order}
    if domain is not None:
        options["domain"] = domain
    if extension is not None:
        options["extension"] = extension
    if modulus is not None:
        options["modulus"] = modulus
    return sp.groebner(generators_tuple, *vars_tuple, **options)


def lex_basis_generators(
    generators: Iterable[sp.Expr],
    variables: Sequence[sp.Symbol],
    *,
    domain: Any | None = None,
    extension: bool | None = None,
) -> tuple[sp.Expr, ...]:
    """Return a lexicographic basis using semialg's triangular variable order.

    Triangular algorithms in semialg pass the ambient variables to the
    polynomial engine in reverse order, making the last listed ambient
    variable largest.  Keeping that convention here avoids duplicate basis
    construction policy in decomposition and real-algebraic solving.
    """
    basis = compute_groebner_basis(
        generators,
        tuple(reversed(tuple(variables))),
        order="lex",
        domain=domain,
        extension=extension,
    )
    return tuple(sp.expand(poly.as_expr()) for poly in basis.polys)


__all__ = ["compute_groebner_basis", "lex_basis_generators"]
