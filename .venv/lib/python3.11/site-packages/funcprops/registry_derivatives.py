"""Native registry identities, gradients, and derivative certificates."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import sympy as sp
from sympy import S

from .exactness import certify_identity_equal
from .registry import FunctionSlot, _replace_slots


@dataclass(frozen=True)
class InstantiatedIdentityRule:
    rule: "IdentityRule"
    expression: sp.Basic
    replacement: Any
    conditions: Any
    substitutions: tuple[tuple[FunctionSlot, Any], ...]


@dataclass(frozen=True)
class IdentityRule:
    """A localized exact identity for a function head."""

    head: Any
    replacement: Any
    conditions: Any = S.true
    source: str | None = None
    arity: int | None = None

    def __post_init__(self):
        object.__setattr__(self, "conditions", sp.sympify(self.conditions))
        if self.arity is not None and (
            isinstance(self.arity, bool)
            or not isinstance(self.arity, int)
            or self.arity < 1
        ):
            raise ValueError("IdentityRule arity must be a positive integer or None")

    @property
    def slots(self):
        used = set()
        for value in (self.replacement, self.conditions):
            if isinstance(value, sp.Basic):
                used.update(value.atoms(FunctionSlot))
        return tuple(sorted(used, key=lambda s: s.position))

    def matches(self, expr) -> bool:
        expr = sp.sympify(expr)
        return (
            getattr(expr, "func", None) == self.head
            and (self.arity is None or len(expr.args) == self.arity)
            and all(s.position < len(expr.args) for s in self.slots)
        )

    def instantiate(self, expr) -> InstantiatedIdentityRule:
        expr = sp.sympify(expr)
        if not self.matches(expr):
            raise ValueError(f"Identity rule for {self.head} does not match {expr}")
        repl = {s: expr.args[s.position] for s in self.slots}
        return InstantiatedIdentityRule(
            self,
            expr,
            _replace_slots(self.replacement, repl),
            _replace_slots(self.conditions, repl),
            tuple((s, repl[s]) for s in self.slots),
        )


@dataclass(frozen=True)
class InstantiatedDerivativeRule:
    rule: "DerivativeRule"
    expression: sp.Basic
    argument: Any
    partial_derivative: Any
    conditions: Any
    substitutions: tuple[tuple[FunctionSlot, Any], ...]


@dataclass(frozen=True)
class DerivativeRule:
    """Exact partial derivative theorem for one positional function argument."""

    head: Any
    with_respect_to: FunctionSlot
    result: Any
    conditions: Any = S.true
    source: str | None = None
    arity: int | None = None

    def __post_init__(self):
        if not isinstance(self.with_respect_to, FunctionSlot):
            raise TypeError("with_respect_to must be a FunctionSlot")
        object.__setattr__(self, "conditions", sp.sympify(self.conditions))
        if self.arity is not None and (
            isinstance(self.arity, bool)
            or not isinstance(self.arity, int)
            or self.arity < 1
        ):
            raise ValueError("DerivativeRule arity must be a positive integer or None")
        positions = [slot.position for slot in self.slots]
        if len(positions) != len(set(positions)):
            raise ValueError("DerivativeRule cannot use distinct slots at one position")

    @property
    def slots(self):
        used = {self.with_respect_to}
        for value in (self.result, self.conditions):
            if isinstance(value, sp.Basic):
                used.update(value.atoms(FunctionSlot))
        return tuple(sorted(used, key=lambda s: s.position))

    def matches(self, expr) -> bool:
        expr = sp.sympify(expr)
        return (
            getattr(expr, "func", None) == self.head
            and (self.arity is None or len(expr.args) == self.arity)
            and all(s.position < len(expr.args) for s in self.slots)
        )

    def instantiate(self, expr) -> InstantiatedDerivativeRule:
        expr = sp.sympify(expr)
        if not self.matches(expr):
            raise ValueError(f"Derivative rule for {self.head} does not match {expr}")
        repl = {s: expr.args[s.position] for s in self.slots}
        return InstantiatedDerivativeRule(
            self,
            expr,
            repl[self.with_respect_to],
            _replace_slots(self.result, repl),
            _replace_slots(self.conditions, repl),
            tuple((s, repl[s]) for s in self.slots),
        )


@dataclass(frozen=True)
class GradientComponent:
    argument_position: int
    argument: Any
    partial_derivative: Any
    conditions: Any = S.true
    source: str | None = None


@dataclass(frozen=True)
class GradientCertificate:
    """Registered outer gradient with one component per requested argument."""

    expression: sp.Basic
    components: tuple[GradientComponent, ...]
    gradient: tuple[Any, ...] | None
    positions: tuple[int, ...]
    conditions: Any
    complete: bool
    missing_positions: tuple[int, ...] = ()
    source: str = "native-registry-gradient"

    def verify(self) -> bool:
        if not self.complete or self.gradient is None:
            return False
        return tuple(c.partial_derivative for c in self.components) == self.gradient


@dataclass(frozen=True)
class ChainRuleTerm:
    argument_position: int
    argument: Any
    partial_derivative: Any
    inner_derivative: Any
    product: Any
    conditions: Any = S.true
    source: str | None = None


@dataclass(frozen=True)
class ChainRuleCertificate:
    """Replayable first-derivative chain-rule construction."""

    expression: sp.Basic
    variable: sp.Symbol
    terms: tuple[ChainRuleTerm, ...]
    derivative: Any
    conditions: Any
    complete: bool
    source: str = "native-registry-chain-rule"
    gradient_certificate: GradientCertificate | None = None

    def verify(self) -> bool:
        if not self.complete:
            return False
        rebuilt = sp.Add(*(term.product for term in self.terms))
        return certify_identity_equal(rebuilt, self.derivative) is True


def registry_gradient(expr, *, registry=None, positions=None) -> GradientCertificate:
    """Build a certified gradient/partial vector from registered rules.

    Conditional strata remain explicit: each requested argument position must
    have exactly one matching derivative theorem. Ambiguous or missing partials
    make the certificate incomplete.
    """
    from .registry import registry as default_registry

    expr = sp.sympify(expr)
    reg = default_registry() if registry is None else registry
    if not getattr(expr, "args", None):
        return GradientCertificate(expr, (), (), (), S.true, True)

    if positions is None:
        positions = tuple(range(len(expr.args)))
    else:
        positions = tuple(positions)
    if any(
        isinstance(position, bool)
        or not isinstance(position, int)
        or position < 0
        or position >= len(expr.args)
        for position in positions
    ):
        raise ValueError("gradient positions must index concrete function arguments")
    if len(set(positions)) != len(positions):
        raise ValueError("gradient positions cannot contain duplicates")

    components = []
    conditions = []
    missing = []
    for position in positions:
        matches = reg.match_derivative_rules(expr, position=position)
        if len(matches) != 1:
            missing.append(position)
            continue
        inst = matches[0].instantiate(expr)
        components.append(
            GradientComponent(
                position,
                expr.args[position],
                inst.partial_derivative,
                inst.conditions,
                inst.rule.source,
            )
        )
        conditions.append(inst.conditions)

    complete = not missing and len(components) == len(positions)
    gradient = (
        tuple(component.partial_derivative for component in components)
        if complete
        else None
    )
    condition = sp.And(*conditions) if conditions else S.true
    return GradientCertificate(
        expr,
        tuple(components),
        gradient,
        positions,
        condition,
        complete,
        tuple(missing),
    )


def registry_derivative(expr, var, *, registry=None) -> ChainRuleCertificate:
    """Differentiate a function application via a registered outer gradient.

    For a path ``x -> F(g1(x), ..., gn(x))`` this certifies the scalar derivative
    as the dot product of the registered gradient of ``F`` with the path tangent.
    """
    expr = sp.sympify(expr)
    var = sp.sympify(var)
    if not isinstance(var, sp.Symbol):
        raise TypeError("var must be a Symbol")
    varying = tuple(i for i, arg in enumerate(expr.args) if arg.has(var))
    if not varying:
        if not expr.has(var):
            return ChainRuleCertificate(expr, var, (), S.Zero, S.true, True)
        # An atomic occurrence such as expr == var is not a function application
        # and therefore has no registered outer gradient.  Mark it incomplete so
        # cached_derivative falls back to SymPy's exact differentiation instead
        # of incorrectly certifying derivative zero.
        return ChainRuleCertificate(expr, var, (), None, S.true, False)

    gradient = registry_gradient(expr, registry=registry, positions=varying)
    if not gradient.complete or gradient.gradient is None:
        return ChainRuleCertificate(
            expr,
            var,
            (),
            None,
            gradient.conditions,
            False,
            gradient_certificate=gradient,
        )

    terms = []
    by_position = {c.argument_position: c for c in gradient.components}
    for position in varying:
        component = by_position[position]
        inner = sp.diff(expr.args[position], var)
        product = sp.simplify(component.partial_derivative * inner)
        terms.append(
            ChainRuleTerm(
                position,
                expr.args[position],
                component.partial_derivative,
                inner,
                product,
                component.conditions,
                component.source,
            )
        )
    derivative = sp.simplify(sp.Add(*(term.product for term in terms)))
    return ChainRuleCertificate(
        expr,
        var,
        tuple(terms),
        derivative,
        gradient.conditions,
        True,
        gradient_certificate=gradient,
    )


def rewrite_with_registry(expr, *, registry=None):
    """Instantiate an unambiguous native identity for ``expr`` if available."""
    from .registry import registry as default_registry

    expr = sp.sympify(expr)
    reg = default_registry() if registry is None else registry
    matches = reg.match_identity_rules(expr)
    if len(matches) != 1:
        return None
    return matches[0].instantiate(expr)


@dataclass(frozen=True)
class JacobianCertificate:
    """Exact Jacobian of a tuple-valued inner map."""

    expressions: tuple[Any, ...]
    variables: tuple[sp.Symbol, ...]
    jacobian: sp.ImmutableMatrix
    complete: bool = True
    source: str = "symbolic-inner-jacobian"

    def verify(self) -> bool:
        expected = sp.ImmutableMatrix(
            [
                [sp.diff(expr, var) for var in self.variables]
                for expr in self.expressions
            ]
        )
        return self.complete and expected == self.jacobian


@dataclass(frozen=True)
class HessianCertificate:
    """Registered outer Hessian derived from exact derivative rules."""

    expression: sp.Basic
    positions: tuple[int, ...]
    hessian: sp.ImmutableMatrix | None
    conditions: Any
    complete: bool
    source: str = "native-registry-hessian"

    def verify(self) -> bool:
        return (
            self.complete
            and self.hessian is not None
            and self.hessian == self.hessian.T
        )


@dataclass(frozen=True)
class MultivariateChainRuleCertificate:
    """Replayable first/second-order chain rule for F(g(x_1,...,x_m))."""

    expression: sp.Basic
    variables: tuple[sp.Symbol, ...]
    gradient: sp.ImmutableMatrix | None
    jacobian: sp.ImmutableMatrix | None
    hessian: sp.ImmutableMatrix | None
    conditions: Any
    complete: bool
    gradient_certificate: GradientCertificate | None = None
    hessian_certificate: HessianCertificate | None = None
    source: str = "native-registry-multivariate-chain-rule"

    def verify(self) -> bool:
        if not self.complete or self.gradient is None or self.jacobian is None:
            return False
        return self.gradient.rows == len(self.variables) and self.jacobian.cols == len(
            self.variables
        )


def registry_jacobian(expressions, variables) -> JacobianCertificate:
    """Build the exact Jacobian of an explicit inner map."""
    expressions = tuple(sp.sympify(expr) for expr in expressions)
    variables = tuple(sp.sympify(var) for var in variables)
    if not variables or not all(isinstance(var, sp.Symbol) for var in variables):
        raise TypeError("variables must be a nonempty tuple of Symbols")
    if len(set(variables)) != len(variables):
        raise ValueError("variables cannot contain duplicates")
    jac = sp.ImmutableMatrix(
        [[sp.diff(expr, var) for var in variables] for expr in expressions]
    )
    return JacobianCertificate(expressions, variables, jac)


def registry_hessian(expr, *, registry=None, positions=None) -> HessianCertificate:
    """Derive an outer Hessian from exact registered first-partial formulas.

    This is intentionally stricter than numerical differentiation: every requested
    first partial must come from one unambiguous DerivativeRule, and its formula is
    differentiated with respect to the rule-local FunctionSlots before concrete
    arguments are substituted.
    """
    from .registry import registry as default_registry

    expr = sp.sympify(expr)
    reg = default_registry() if registry is None else registry
    if positions is None:
        positions = tuple(range(len(expr.args)))
    else:
        positions = tuple(positions)
    if not positions:
        return HessianCertificate(expr, (), sp.ImmutableMatrix([]), S.true, True)

    rows = []
    conditions = []
    for i in positions:
        matches = reg.derivative_rules(head=expr.func, position=i)
        matches = tuple(rule for rule in matches if rule.matches(expr))
        if len(matches) != 1:
            return HessianCertificate(expr, positions, None, sp.And(*conditions), False)
        rule = matches[0]
        repl = {slot: expr.args[slot.position] for slot in rule.slots}
        conditions.append(_replace_slots(rule.conditions, repl))
        row = []
        for j in positions:
            slot = next((s for s in rule.slots if s.position == j), None)
            if slot is None:
                # The derivative formula is independent of an unmentioned slot.
                value = S.Zero
            else:
                value = _replace_slots(sp.diff(rule.result, slot), repl)
            row.append(sp.simplify(value))
        rows.append(row)
    matrix = sp.ImmutableMatrix(rows)
    # Mixed partial agreement is part of certification; do not symmetrize guesses.
    complete = matrix == matrix.T
    return HessianCertificate(
        expr, positions, matrix if complete else None, sp.And(*conditions), complete
    )


def registry_chain_rule(
    expr, variables, *, registry=None, second_order=False
) -> MultivariateChainRuleCertificate:
    """Transport registered derivatives through a multivariate inner map."""
    expr = sp.sympify(expr)
    variables = tuple(sp.sympify(v) for v in variables)
    if not variables or not all(isinstance(v, sp.Symbol) for v in variables):
        raise TypeError("variables must be a nonempty tuple of Symbols")
    varying = tuple(
        i for i, arg in enumerate(expr.args) if any(arg.has(v) for v in variables)
    )
    inner = tuple(expr.args[i] for i in varying)
    jac_cert = registry_jacobian(inner, variables)
    grad_cert = registry_gradient(expr, registry=registry, positions=varying)
    if not grad_cert.complete or grad_cert.gradient is None:
        return MultivariateChainRuleCertificate(
            expr,
            variables,
            None,
            jac_cert.jacobian,
            None,
            grad_cert.conditions,
            False,
            gradient_certificate=grad_cert,
        )
    outer_grad = sp.ImmutableMatrix([grad_cert.gradient])  # 1 x n
    transported_grad = sp.ImmutableMatrix(outer_grad * jac_cert.jacobian).T
    if not second_order:
        return MultivariateChainRuleCertificate(
            expr,
            variables,
            transported_grad,
            jac_cert.jacobian,
            None,
            grad_cert.conditions,
            True,
            gradient_certificate=grad_cert,
        )

    hess_cert = registry_hessian(expr, registry=registry, positions=varying)
    if not hess_cert.complete or hess_cert.hessian is None:
        return MultivariateChainRuleCertificate(
            expr,
            variables,
            transported_grad,
            jac_cert.jacobian,
            None,
            sp.And(grad_cert.conditions, hess_cert.conditions),
            False,
            gradient_certificate=grad_cert,
            hessian_certificate=hess_cert,
        )
    J = jac_cert.jacobian
    composite_hessian = J.T * hess_cert.hessian * J
    for k, position in enumerate(varying):
        inner_hessian = sp.hessian(expr.args[position], variables)
        composite_hessian += grad_cert.gradient[k] * inner_hessian
    composite_hessian = sp.ImmutableMatrix(composite_hessian.applyfunc(sp.simplify))
    return MultivariateChainRuleCertificate(
        expr,
        variables,
        transported_grad,
        J,
        composite_hessian,
        sp.And(grad_cert.conditions, hess_cert.conditions),
        True,
        gradient_certificate=grad_cert,
        hessian_certificate=hess_cert,
    )
