"""Lazy shared cache for reusable exact property proofs.

`PropertyProofContext` is the single internal entry point for proof primitives.
Every slot remains lazy: constructing a context does no mathematical work.
Recursive algorithms still call their uncached implementations where doing so
avoids cache-key overhead on transient subexpressions.
"""

from __future__ import annotations

import operator
from dataclasses import dataclass
from functools import lru_cache

import sympy as sp
from sympy import S
from sympy.core.function import PoleError

from .assumption_normalization import normalize_assumptions, normalize_domain


@dataclass(frozen=True)
class PropertyProofKey:
    expression: sp.Basic
    variable: sp.Basic
    domain: object = S.true
    assumptions: sp.Basic = S.true


def proof_key(expr, var, *, domain=S.true, assumptions=S.true):
    return PropertyProofKey(
        sp.sympify(expr),
        sp.sympify(var),
        normalize_domain(domain),
        normalize_assumptions(assumptions),
    )


@lru_cache(maxsize=2048)
def _derivative_cached(expr, var, order):
    return sp.diff(expr, var, order)


@lru_cache(maxsize=1024)
def _factored_derivative_cached(expr, var, order):
    return sp.factor(_derivative_cached(expr, var, order))


def _derivative_order(order):
    """Return a concrete nonnegative derivative order without truncation."""
    if isinstance(order, bool):
        raise ValueError("order must be a nonnegative integer")
    try:
        value = operator.index(order)
    except TypeError as exc:
        raise ValueError("order must be a nonnegative integer") from exc
    if value < 0:
        raise ValueError("order must be a nonnegative integer")
    return value


def _uncached_after_key_error(cache_fn, args):
    """Retry an LRU-backed call only when its arguments are unhashable."""
    try:
        hash(args)
    except TypeError:
        return cache_fn.__wrapped__(*args)
    raise


def cached_derivative(expr, var, *, order=1, factored=False):
    expr = sp.sympify(expr)
    var = sp.sympify(var)
    order = _derivative_order(order)
    if order == 1:
        from .registry_derivatives import registry_derivative

        cert = registry_derivative(expr, var)
        if cert.complete and cert.conditions is S.true:
            return sp.factor(cert.derivative) if factored else cert.derivative
    fn = _factored_derivative_cached if factored else _derivative_cached
    args = (expr, var, order)
    try:
        return fn(*args)
    except TypeError:
        return _uncached_after_key_error(fn, args)


@lru_cache(maxsize=1024)
def _domain_cached(expr, var, domain_kind, stratified):
    from .real_analysis.domain import _function_domain_impl

    dom = S.Reals if domain_kind == "real" else S.Complexes
    return _function_domain_impl(expr, var, domain=dom, stratified=stratified)


def cached_domain(expr, var, *, domain=S.Reals, stratified=False):
    if domain in (S.Reals, "real", "reals", "R", "RR"):
        kind = "real"
    elif domain in (S.Complexes, "complex", "complexes", "C", "CC"):
        kind = "complex"
    else:
        raise ValueError("domain must be Reals or Complexes")
    return _domain_cached(sp.sympify(expr), sp.sympify(var), kind, bool(stratified))


@lru_cache(maxsize=1024)
def _sign_cached(expr, var, domain, strict_inequalities):
    from .real_analysis.sign import _function_sign_public_impl

    return _function_sign_public_impl(
        expr, var, domain=domain, strict_inequalities=strict_inequalities
    )


def cached_sign(expr, var, *, domain=S.true, strict_inequalities=True):
    key = proof_key(expr, var, domain=domain)
    args = (key.expression, key.variable, key.domain, bool(strict_inequalities))
    try:
        return _sign_cached(*args)
    except TypeError:
        try:
            hash(args)
        except TypeError:
            from .real_analysis.sign import _function_sign_public_impl

            return _function_sign_public_impl(
                key.expression,
                key.variable,
                domain=key.domain,
                strict_inequalities=strict_inequalities,
            )
        raise


@lru_cache(maxsize=1024)
def _endpoint_limit_cached(expr, var, point, direction, attained):
    try:
        if point is S.NegativeInfinity:
            return sp.limit(expr, var, -sp.oo)
        if point is S.Infinity:
            return sp.limit(expr, var, sp.oo)
        if attained:
            return sp.simplify(expr.subs(var, point))
        return sp.limit(expr, var, point, dir=direction)
    except (PoleError, TypeError, ValueError, NotImplementedError):
        return None


def cached_endpoint_limit(expr, var, point, *, direction="+", attained=False):
    expression = sp.sympify(expr)
    variable = sp.sympify(var)
    point = sp.sympify(point)
    args = (expression, variable, point, str(direction), bool(attained))
    try:
        return _endpoint_limit_cached(*args)
    except TypeError:
        return _uncached_after_key_error(_endpoint_limit_cached, args)


@lru_cache(maxsize=1024)
def _zero_geometry_cached(expr, var, derivative_order, domain):
    from .transcendental.zero_geometry import _zero_geometry_impl

    return _zero_geometry_impl(
        expr, var, derivative_order=derivative_order, domain=domain
    )


def cached_zero_geometry(expr, var, *, derivative_order=0, domain=S.Reals):
    key = proof_key(expr, var, domain=domain)
    order = _derivative_order(derivative_order)
    args = (key.expression, key.variable, order, key.domain)
    try:
        return _zero_geometry_cached(*args)
    except TypeError:
        return _uncached_after_key_error(_zero_geometry_cached, args)


@lru_cache(maxsize=512)
def _singularity_cached(expr, var, domain, assumptions, ignore_unknown):
    from .singularities import _singularity_analysis_impl

    return _singularity_analysis_impl(
        expr, var, domain=domain, assumptions=assumptions, ignore_unknown=ignore_unknown
    )


def cached_singularity_analysis(
    expr, var, *, domain=S.Complexes, assumptions=S.true, ignore_unknown=False
):
    key = proof_key(expr, var, domain=domain, assumptions=assumptions)
    args = (
        key.expression,
        key.variable,
        key.domain,
        key.assumptions,
        bool(ignore_unknown),
    )
    try:
        return _singularity_cached(*args)
    except TypeError:
        return _uncached_after_key_error(_singularity_cached, args)


@lru_cache(maxsize=1024)
def _proof_context_cached(expr, var, domain, assumptions):
    return PropertyProofContext(expr, var, domain, assumptions)


@dataclass(frozen=True)
class PropertyProofContext:
    """Lazy normalized proof context shared by all global property APIs."""

    expression: sp.Basic
    variable: sp.Basic
    domain: object = S.true
    assumptions: sp.Basic = S.true

    def derivative(self, order=1, *, factored=False):
        return cached_derivative(
            self.expression, self.variable, order=order, factored=factored
        )

    def natural_domain(self, *, complex_domain=False, stratified=False):
        """Return the cached natural real or complex domain."""
        target = S.Complexes if complex_domain else S.Reals
        return cached_domain(
            self.expression,
            self.variable,
            domain=target,
            stratified=stratified,
        )

    def sign(self, *, domain=None, strict_inequalities=True):
        return cached_sign(
            self.expression,
            self.variable,
            domain=self.domain if domain is None else domain,
            strict_inequalities=strict_inequalities,
        )

    def endpoint_limit(self, point, *, direction="+", attained=False):
        return cached_endpoint_limit(
            self.expression,
            self.variable,
            point,
            direction=direction,
            attained=attained,
        )

    def zero_geometry(self, *, derivative_order=0, domain=None):
        return cached_zero_geometry(
            self.expression,
            self.variable,
            derivative_order=derivative_order,
            domain=self.domain if domain is None else domain,
        )

    def singularities(self, *, domain=S.Complexes, ignore_unknown=False):
        return cached_singularity_analysis(
            self.expression,
            self.variable,
            domain=domain,
            assumptions=self.assumptions,
            ignore_unknown=ignore_unknown,
        )

    def monotonic_partition(self, *, domain=None):
        from .real_analysis.monotonic_certificates import monotonic_partition

        return monotonic_partition(
            self.expression,
            self.variable,
            domain=self.domain if domain is None else domain,
            return_result=True,
        )


def property_proof_context(expr, var, *, domain=S.true, assumptions=S.true):
    key = proof_key(expr, var, domain=domain, assumptions=assumptions)
    args = (key.expression, key.variable, key.domain, key.assumptions)
    try:
        return _proof_context_cached(*args)
    except TypeError:
        try:
            hash(args)
        except TypeError:
            return PropertyProofContext(*args)
        raise


def clear_property_proof_cache():
    _derivative_cached.cache_clear()
    _factored_derivative_cached.cache_clear()
    _domain_cached.cache_clear()
    _sign_cached.cache_clear()
    _endpoint_limit_cached.cache_clear()
    _zero_geometry_cached.cache_clear()
    _singularity_cached.cache_clear()
    _proof_context_cached.cache_clear()
    from .real_analysis.monotonic_certificates import _monotonic_partition_on_set

    _monotonic_partition_on_set.cache_clear()


def property_proof_cache_info():
    return {
        "derivative": _derivative_cached.cache_info(),
        "factored_derivative": _factored_derivative_cached.cache_info(),
        "domain": _domain_cached.cache_info(),
        "sign": _sign_cached.cache_info(),
        "endpoint_limit": _endpoint_limit_cached.cache_info(),
        "zero_geometry": _zero_geometry_cached.cache_info(),
        "singularity": _singularity_cached.cache_info(),
        "context": _proof_context_cached.cache_info(),
    }
