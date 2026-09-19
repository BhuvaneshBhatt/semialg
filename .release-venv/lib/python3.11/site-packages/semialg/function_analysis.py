"""Exact convexity and monotonicity classification for semialgebraic functions."""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field, replace

import sympy as sp

from .conditional import ParameterStratifiedResult
from .context import SemialgebraicContext
from .convexity import ConvexityCertificate, convexity_certificate
from .domain_solve import function_domain, normalize_domain_sensitive_constraints
from .formula import parse_formula
from .function_graph import (
    has_semialgebraic_graph_special,
    semialgebraic_formula_graph,
    semialgebraic_function_graph,
)
from .internal_symbols import fresh_real_dummy
from .matrix_analysis import MatrixDefinitenessResult
from .normalization import normalize_formula, normalize_variables
from .qe import qe_by_complete_cad
from .reasoning_signs import function_sign
from .structural_keys import symbol_identity_key


@dataclass(frozen=True)
class FunctionConvexityResult:
    """Exact convexity/concavity classification for one function on one domain."""

    classification: str
    expression: sp.Expr
    variables: tuple[sp.Symbol, ...]
    domain: sp.Expr
    natural_domain: sp.Expr
    domain_convex: bool | None
    convex: bool | None
    concave: bool | None
    method: str
    domain_certificate: ConvexityCertificate | None = None
    convexity_certificate: MatrixDefinitenessResult | None = None
    concavity_certificate: MatrixDefinitenessResult | None = None
    counterexample: Mapping[str, object] | None = None
    details: Mapping[str, object] = field(default_factory=dict)
    strictly_convex: bool | None = None
    strictly_concave: bool | None = None
    strongly_convex: bool | None = None
    strongly_concave: bool | None = None
    strong_convexity_modulus: sp.Expr | None = None
    strong_concavity_modulus: sp.Expr | None = None
    quasiconvex: bool | None = None
    quasiconcave: bool | None = None
    strictly_quasiconvex: bool | None = None
    strictly_quasiconcave: bool | None = None
    pseudoconvex: bool | None = None
    pseudoconcave: bool | None = None
    log_convex: bool | None = None
    log_concave: bool | None = None

    @property
    def certified(self) -> bool:
        """Return whether the primary classification is exact and decisive."""

        return self.classification != "unknown"


@dataclass(frozen=True)
class FunctionPropertyPartitionResult:
    """Exact univariate partition into regions with a certified function property."""

    expression: sp.Expr
    variable: sp.Symbol
    domain: sp.Expr
    natural_domain: sp.Expr
    property_name: str
    pieces: tuple[tuple[str, sp.Expr], ...]
    method: str
    diagnostics: Mapping[str, object] = field(default_factory=dict)

    @property
    def certified(self) -> bool:
        """Return whether every piece has a decisive classification."""

        return all(classification != "unknown" for classification, _ in self.pieces)


@dataclass
class _FunctionAnalysisContext:
    """Shared lazy exact products for one nonparametric function-property query."""

    expression: sp.Expr
    variables: tuple[sp.Symbol, ...]
    explicit_domain: sp.Expr
    _natural_domain: sp.Expr | None = None
    _effective_domain: sp.Expr | None = None
    _semialgebraic_context: SemialgebraicContext | None = None
    _derivatives: dict[tuple[sp.Symbol, int], sp.Expr] = field(default_factory=dict)
    _hessian: sp.ImmutableMatrix | None = None
    _gradient: sp.ImmutableMatrix | None = None
    _signs: dict[tuple[sp.Expr, sp.Expr], str] = field(default_factory=dict)
    _monotonic_partitions: dict[sp.Symbol, FunctionPropertyPartitionResult] = field(
        default_factory=dict
    )
    _smoothness: dict[int, object] = field(default_factory=dict)
    _domain_convexity: ConvexityCertificate | None = None

    @property
    def natural_domain(self) -> sp.Expr:
        if self._natural_domain is None:
            self._natural_domain = function_domain(self.expression, self.variables)
        return self._natural_domain

    @property
    def domain(self) -> sp.Expr:
        if self._effective_domain is None:
            self._effective_domain = _algebraic_domain_formula(
                normalize_formula(sp.And(self.explicit_domain, self.natural_domain)),
                self.variables,
            )
        return self._effective_domain

    @property
    def semialgebraic_context(self) -> SemialgebraicContext:
        if self._semialgebraic_context is None:
            self._semialgebraic_context = SemialgebraicContext(self.domain, self.variables)
        return self._semialgebraic_context

    def derivative(self, variable: sp.Symbol, order: int = 1) -> sp.Expr:
        key = (variable, order)
        if key not in self._derivatives:
            self._derivatives[key] = sp.diff(self.expression, variable, order)
        return self._derivatives[key]

    @property
    def gradient(self) -> sp.ImmutableMatrix:
        if self._gradient is None:
            self._gradient = sp.ImmutableMatrix(
                [sp.diff(self.expression, v) for v in self.variables]
            )
        return self._gradient

    @property
    def hessian(self) -> sp.ImmutableMatrix:
        if self._hessian is None:
            self._hessian = sp.ImmutableMatrix(sp.hessian(self.expression, self.variables))
        return self._hessian

    def sign(self, expression: sp.Expr, *, domain: sp.Expr | None = None) -> str:
        sign_domain = self.domain if domain is None else domain
        key = (sp.sympify(expression), sign_domain)
        if key not in self._signs:
            self._signs[key] = function_sign(expression, self.variables, assumptions=sign_domain)
        return self._signs[key]

    def domain_convexity(self) -> ConvexityCertificate:
        if self._domain_convexity is None:
            self._domain_convexity = convexity_certificate(self.domain, self.variables)
        return self._domain_convexity

    def monotonic_partition(self, variable: sp.Symbol) -> FunctionPropertyPartitionResult:
        if variable not in self._monotonic_partitions:
            from ._function_analysis_partitions import _partition_from_derivative_sign

            self._monotonic_partitions[variable] = _partition_from_derivative_sign(
                self.expression,
                variable,
                self.domain,
                derivative_order=1,
                property_name="monotonicity",
            )
        return self._monotonic_partitions[variable]

    def smoothness(self, max_order: int = 1):
        if max_order not in self._smoothness:
            from .function_properties import function_smoothness

            self._smoothness[max_order] = function_smoothness(
                self.expression, self.variables, domain=self.domain, max_order=max_order
            )
        return self._smoothness[max_order]

    def cache_diagnostics(self) -> dict[str, int]:
        """Return deterministic counts of lazily cached analysis products."""

        return {
            "derivatives": len(self._derivatives),
            "sign_queries": len(self._signs),
            "monotonic_partitions": len(self._monotonic_partitions),
            "smoothness_queries": len(self._smoothness),
        }


@dataclass(frozen=True)
class FunctionMonotonicityResult:
    """Exact monotonicity classification for one univariate function."""

    classification: str
    expression: sp.Expr
    variable: sp.Symbol
    domain: sp.Expr
    natural_domain: sp.Expr
    increasing: bool | None
    decreasing: bool | None
    strictly_increasing: bool | None
    strictly_decreasing: bool | None
    constant: bool | None
    method: str
    counterexample: Mapping[str, object] | None = None
    details: Mapping[str, object] = field(default_factory=dict)

    @property
    def certified(self) -> bool:
        """Return whether the primary classification is exact and decisive."""

        return self.classification != "unknown"


def _classify_monotonicity(
    *,
    increasing: bool | None,
    decreasing: bool | None,
    strictly_increasing: bool | None,
    strictly_decreasing: bool | None,
    constant: bool | None,
) -> str:
    if constant is True:
        return "constant"
    if strictly_increasing is True:
        return "strictly_increasing"
    if strictly_decreasing is True:
        return "strictly_decreasing"
    if increasing is True:
        return "increasing"
    if decreasing is True:
        return "decreasing"
    if increasing is False and decreasing is False:
        return "nonmonotonic"
    return "unknown"


def _classify(convex: bool | None, concave: bool | None) -> str:
    if convex is True and concave is True:
        return "affine"
    if convex is True:
        return "convex"
    if concave is True:
        return "concave"
    if convex is False and concave is False:
        return "neither"
    return "unknown"


def _normalize_function_variables(
    expression: sp.Expr,
    domain: sp.Expr,
    variables: Sequence[sp.Symbol | str] | None,
    parameters: Sequence[sp.Symbol | str] | None,
) -> tuple[tuple[sp.Symbol, ...], tuple[sp.Symbol, ...]]:
    combined = sp.Tuple(expression, domain)
    normalized_parameters = (
        normalize_variables(parameters, combined, append_context_symbols=False)
        if parameters is not None
        else ()
    )
    if variables is None:
        parameter_set = set(normalized_parameters)
        normalized_variables = tuple(
            sorted(combined.free_symbols - parameter_set, key=symbol_identity_key)
        )
    else:
        normalized_variables = normalize_variables(
            variables, combined, append_context_symbols=False
        )
    if set(normalized_variables) & set(normalized_parameters):
        raise ValueError("variables and parameters must be disjoint")
    return normalized_variables, normalized_parameters


def _truth_condition(result: ParameterStratifiedResult) -> sp.Expr:
    conditions = tuple(branch.condition for branch in result.branches if branch.value is True)
    return sp.Or(*conditions) if conditions else sp.false


def _algebraic_domain_formula(
    domain: sp.Expr,
    original_symbols: tuple[sp.Symbol, ...],
) -> sp.Expr:
    """Return an exact polynomial/Boolean formula in the original coordinates."""

    normalized = normalize_domain_sensitive_constraints(domain, original_symbols).formula
    if not has_semialgebraic_graph_special(normalized):
        return normalized
    graph = semialgebraic_formula_graph(normalized)
    if not graph.auxiliary_variables:
        return graph.formula
    free_symbols = tuple(
        sorted(
            set(original_symbols) | (graph.formula.free_symbols - set(graph.auxiliary_variables)),
            key=symbol_identity_key,
        )
    )
    result = qe_by_complete_cad(
        (*free_symbols, *graph.auxiliary_variables),
        tuple(("exists", auxiliary) for auxiliary in graph.auxiliary_variables),
        parse_formula(graph.formula),
        free_variables=free_symbols,
        return_result=True,
    )
    return normalize_formula(result.formula)


def _copy_variables(variables: tuple[sp.Symbol, ...], prefix: str) -> tuple[sp.Symbol, ...]:
    return tuple(fresh_real_dummy(f"{prefix}_{variable.name}") for variable in variables)


def _copy_expression(
    expression: sp.Expr,
    variables: tuple[sp.Symbol, ...],
    copied_variables: tuple[sp.Symbol, ...],
) -> sp.Expr:
    return expression.xreplace(dict(zip(variables, copied_variables, strict=True)))


def _domain_graph(formula: sp.Expr) -> tuple[sp.Expr, tuple[sp.Symbol, ...]]:
    graph = semialgebraic_formula_graph(formula)
    return graph.formula, graph.auxiliary_variables


def _function_graph(
    expression: sp.Expr, target_name: str
) -> tuple[sp.Expr, sp.Symbol, tuple[sp.Symbol, ...]]:
    target = fresh_real_dummy(target_name)
    graph = semialgebraic_function_graph(expression, target)
    return graph.formula, target, graph.auxiliary_variables


from ._function_analysis_partitions import (  # noqa: E402
    function_convex_partition,
    function_monotonic_partition,
    function_sign_partition,
)


def function_monotonicity(
    expression,
    variable: sp.Symbol | str | None = None,
    *,
    domain=sp.true,
    parameters: Sequence[sp.Symbol | str] | None = None,
    return_result: bool = False,
) -> str | FunctionMonotonicityResult | ParameterStratifiedResult:
    """Classify exact monotonicity of a supported semialgebraic function.

    The strongest canonical result is returned: ``'constant'``,
    ``'strictly_increasing'``, ``'increasing'``, ``'strictly_decreasing'``,
    ``'decreasing'``, ``'nonmonotonic'``, or ``'unknown'``.  The exact real
    function domain is intersected automatically with ``domain``.

    Monotonicity is with respect to one real variable.  Other free symbols are
    automatically treated as parameters when ``parameters`` is omitted, so a
    conditional :class:`ParameterStratifiedResult` is produced whenever the
    answer depends on them.  Derivative sign analysis is a fast sufficient
    path; exact pairwise graph/QE checks provide the defining fallback and
    correctly certify strict cases such as ``x**3`` despite a zero derivative
    at isolated points.
    """

    from ._function_analysis_monotonicity import (
        _function_monotonicity_with_parameters,
        _function_monotonicity_without_parameters,
    )

    normalized_expression = sp.sympify(expression)
    normalized_domain = normalize_formula(domain)
    combined = sp.Tuple(normalized_expression, normalized_domain)
    if variable is None:
        explicit_parameters = (
            normalize_variables(parameters, combined, append_context_symbols=False)
            if parameters is not None
            else ()
        )
        candidates = tuple(
            sorted(combined.free_symbols - set(explicit_parameters), key=symbol_identity_key)
        )
        if len(candidates) != 1:
            raise ValueError(
                "variable must be supplied unless exactly one non-parameter symbol is present"
            )
        normalized_variable = candidates[0]
    else:
        normalized_variable = normalize_variables(
            (variable,), combined, append_context_symbols=False
        )[0]
    if parameters is None:
        normalized_parameters = tuple(
            sorted(combined.free_symbols - {normalized_variable}, key=symbol_identity_key)
        )
    else:
        normalized_parameters = normalize_variables(
            parameters, combined, append_context_symbols=False
        )
    if normalized_variable in normalized_parameters:
        raise ValueError("variable and parameters must be disjoint")
    if normalized_parameters:
        return _function_monotonicity_with_parameters(
            normalized_expression,
            normalized_variable,
            normalized_parameters,
            normalized_domain,
        )
    result = _function_monotonicity_without_parameters(
        normalized_expression, normalized_variable, normalized_domain
    )
    return result if return_result else result.classification


from ._function_analysis_properties import _augment_convexity_result  # noqa: E402


def function_convexity(
    expression,
    variables: Sequence[sp.Symbol | str] | None = None,
    *,
    domain=sp.true,
    parameters: Sequence[sp.Symbol | str] | None = None,
    return_result: bool = False,
    properties: str = "primary",
) -> str | FunctionConvexityResult | ParameterStratifiedResult:
    """Classify exact convexity/concavity of a supported semialgebraic function.

    The primary result is the strongest certified canonical classification,
    including ``'strongly_convex'``, ``'strictly_convex'``, ``'convex'``,
    ``'affine'``, and the corresponding concave classifications, plus
    ``'neither'``, ``'nonconvex_domain'``, and ``'unknown'``.  The natural real
    domain recognized by :func:`function_domain` is automatically intersected
    with ``domain``.

    Polynomial functions use exact Hessian matrix definiteness first.  Safe
    affine equality presolve makes this relative to lower-dimensional affine
    domains where possible.  Any unresolved side is decided by the exact
    Jensen definition using semialgebraic function graphs and QE, which also
    supports nondifferentiable expressions such as ``Abs``, ``Min``, ``Max``,
    rational powers, real roots, and finite ``Piecewise`` expressions.

    With ``parameters`` supplied, the result is a certified
    :class:`ParameterStratifiedResult` whose branch values are the same primary
    classifications.  Polynomial full-dimensional problems use matrix
    definiteness stratification; other supported cases use quantified Jensen
    conditions.  Set ``properties="all"`` with ``return_result=True`` to
    additionally request quasi-/strict-quasi-, pseudo-, and algebraically
    provable log-convex/log-concave properties.
    """

    from ._function_analysis_convexity import (
        _function_convexity_with_parameters,
        _function_convexity_without_parameters,
    )

    normalized_expression = sp.sympify(expression)
    explicit_domain_norm = normalize_formula(domain)
    normalized_variables, normalized_parameters = _normalize_function_variables(
        normalized_expression,
        explicit_domain_norm,
        variables,
        parameters,
    )
    if parameters is None and variables is not None:
        normalized_parameters = tuple(
            sorted(
                (
                    sp.Tuple(normalized_expression, explicit_domain_norm).free_symbols
                    - set(normalized_variables)
                ),
                key=symbol_identity_key,
            )
        )
    if normalized_parameters:
        return _function_convexity_with_parameters(
            normalized_expression,
            normalized_variables,
            normalized_parameters,
            explicit_domain_norm,
        )
    if properties not in {"primary", "all"}:
        raise ValueError("properties must be 'primary' or 'all'")
    analysis = _FunctionAnalysisContext(
        normalized_expression, normalized_variables, explicit_domain_norm
    )
    result = _function_convexity_without_parameters(
        normalized_expression, normalized_variables, explicit_domain_norm, analysis=analysis
    )
    result = _augment_convexity_result(
        result,
        normalized_expression,
        normalized_variables,
        properties=properties,
        analysis=analysis,
    )
    result = replace(
        result,
        details={**result.details, "analysis_cache": analysis.cache_diagnostics()},
    )
    return result if return_result else result.classification


__all__ = [
    "FunctionConvexityResult",
    "FunctionMonotonicityResult",
    "FunctionPropertyPartitionResult",
    "function_convexity",
    "function_convex_partition",
    "function_monotonicity",
    "function_monotonic_partition",
    "function_sign_partition",
]
