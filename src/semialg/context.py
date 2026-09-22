from __future__ import annotations

from collections.abc import Hashable, Iterator
from contextlib import contextmanager
from contextvars import ContextVar
from dataclasses import dataclass, field
from functools import wraps
from threading import RLock


@dataclass
class ExactComputationContext:
    """Per-operation cache shared by exact CAD/algebraic subroutines.

    Context entries live only for one top-level solve/optimization/range call.
    Process-local bounded LRUs remain a second-level cache, but repeated work
    inside a solve is resolved here first and disappears when the operation
    finishes.  This avoids global cache pollution while allowing nested solver
    calls to reuse projection, root, sign, comparison, specialization and RUR
    results.
    """

    caches: dict[str, dict[Hashable, object]] = field(default_factory=dict)
    hits: dict[str, int] = field(default_factory=dict)
    misses: dict[str, int] = field(default_factory=dict)
    _lock: RLock = field(default_factory=RLock, init=False, repr=False)

    @staticmethod
    def _logical_namespace(namespace: Hashable) -> str:
        if isinstance(namespace, tuple) and len(namespace) == 2 and isinstance(namespace[0], str):
            return namespace[0]
        return str(namespace)

    def get(self, namespace: Hashable, key: Hashable) -> tuple[bool, object]:
        logical = self._logical_namespace(namespace)
        with self._lock:
            bucket = self.caches.get(namespace)
            if bucket is not None and key in bucket:
                self.hits[logical] = self.hits.get(logical, 0) + 1
                return True, bucket[key]
            self.misses[logical] = self.misses.get(logical, 0) + 1
            return False, None

    def put(self, namespace: Hashable, key: Hashable, value: object) -> None:
        with self._lock:
            self.caches.setdefault(namespace, {})[key] = value

    def cache_size(self, namespace: str | None = None) -> int:
        with self._lock:
            if namespace is not None:
                return sum(
                    len(bucket)
                    for key, bucket in self.caches.items()
                    if self._logical_namespace(key) == namespace
                )
            return sum(len(bucket) for bucket in self.caches.values())

    def prune_namespace(self, namespace: str, generation: int) -> None:
        """Discard mirrors for older generations of one process cache."""

        with self._lock:
            stale = [
                key
                for key in self.caches
                if (
                    isinstance(key, tuple)
                    and len(key) == 2
                    and key[0] == namespace
                    and key[1] != generation
                )
            ]
            for key in stale:
                self.caches.pop(key, None)

    def stats(self) -> dict[str, dict[str, int]]:
        with self._lock:
            sizes: dict[str, int] = {}
            for key, bucket in self.caches.items():
                logical = self._logical_namespace(key)
                sizes[logical] = sizes.get(logical, 0) + len(bucket)
            names = set(sizes) | set(self.hits) | set(self.misses)
            return {
                name: {
                    "hits": self.hits.get(name, 0),
                    "misses": self.misses.get(name, 0),
                    "size": sizes.get(name, 0),
                }
                for name in sorted(names)
            }


_CURRENT_CONTEXT: ContextVar[ExactComputationContext | None] = ContextVar(
    "semialg_exact_computation_context", default=None
)


def current_computation_context() -> ExactComputationContext | None:
    """Return the active exact-computation context, if any."""

    return _CURRENT_CONTEXT.get()


@contextmanager
def computation_context(
    context: ExactComputationContext | None = None,
) -> Iterator[ExactComputationContext]:
    """Create/reuse a context for one complete exact solve operation.

    Nested calls automatically reuse the active context.  Passing an explicit
    context installs it only when no context is already active.
    """

    active = _CURRENT_CONTEXT.get()
    if active is not None:
        yield active
        return
    owned = context if context is not None else ExactComputationContext()
    token = _CURRENT_CONTEXT.set(owned)
    try:
        yield owned
    finally:
        _CURRENT_CONTEXT.reset(token)


def with_computation_context(function):
    """Decorator that gives a top-level solver call one reusable cache scope."""

    @wraps(function)
    def wrapped(*args, **kwargs):
        with computation_context():
            return function(*args, **kwargs)

    return wrapped


def context_cache_get(namespace: Hashable, key: Hashable) -> tuple[bool, object]:
    context = _CURRENT_CONTEXT.get()
    if context is None:
        return False, None
    return context.get(namespace, key)


def context_cache_put(namespace: Hashable, key: Hashable, value: object) -> None:
    context = _CURRENT_CONTEXT.get()
    if context is not None:
        context.put(namespace, key, value)


__all__ = [
    "ExactComputationContext",
    "SemialgebraicContext",
    "computation_context",
    "current_computation_context",
    "with_computation_context",
]


@dataclass(frozen=True)
class SemialgebraicContext:
    """Reusable normalized semialgebraic problem plus exact computation cache.

    The context stores only immutable/certified structural data and
    lazily computed solver artifacts. Existing function APIs remain usable; the
    context provides an explicit lifecycle for callers performing several
    queries about the same region or quantified matrix.  ``variables=None``
    infers variables from the formula; an explicit empty sequence denotes a
    zero-dimensional variable list.
    """

    formula: object
    variables: tuple[object, ...] | None = None
    computation: ExactComputationContext = field(default_factory=ExactComputationContext)
    _memo: dict[Hashable, object] = field(default_factory=dict, init=False, repr=False)

    def __post_init__(self) -> None:
        import sympy as sp

        from .normalization import normalize_formula, normalize_variables

        expr = normalize_formula(self.formula)
        vars_ = normalize_variables(
            self.variables, expr, append_context_symbols=self.variables is None
        )
        object.__setattr__(self, "formula", expr)
        object.__setattr__(self, "variables", tuple(vars_))
        # Ensure callers cannot smuggle non-Symbol variable objects through the
        # dataclass constructor while retaining string normalization support.
        if not all(isinstance(v, sp.Symbol) for v in self.variables):
            raise TypeError("SemialgebraicContext variables must normalize to SymPy Symbols")

    @property
    def parsed_formula(self):
        from .formula import parse_formula

        if "parsed_formula" not in self._memo:
            self._memo["parsed_formula"] = parse_formula(self.formula)
        return self._memo["parsed_formula"]

    @property
    def polynomials(self):
        from .formula import formula_polynomials

        if "polynomials" not in self._memo:
            self._memo["polynomials"] = tuple(formula_polynomials(self.parsed_formula))
        return self._memo["polynomials"]

    @property
    def equational_constraints(self):
        from .formula import equational_constraints

        if "equational_constraints" not in self._memo:
            self._memo["equational_constraints"] = tuple(
                equational_constraints(self.parsed_formula)
            )
        return self._memo["equational_constraints"]

    @property
    def incidence(self):
        from .incidence import analyze_incidence

        if "incidence" not in self._memo:
            self._memo["incidence"] = analyze_incidence(self.polynomials, self.variables)
        return self._memo["incidence"]

    @property
    def presolved(self):
        from .presolve import presolve_semialgebraic

        if "presolved" not in self._memo:
            self._memo["presolved"] = presolve_semialgebraic(self.formula, self.variables)
        return self._memo["presolved"]

    @property
    def projection_tower(self):
        from .cad_algorithms.projection.collins import build_collins_proj_set

        if "projection_tower" not in self._memo:
            with computation_context(self.computation):
                self._memo["projection_tower"] = build_collins_proj_set(
                    self.polynomials or (1,), self.variable_order
                )
        return self._memo["projection_tower"]

    def exceptional_parameters(self, parameters):
        from .decomposition._parametric_support import analyze_parametric_boundaries

        key = ("exceptional_parameters", tuple(parameters))
        if key not in self._memo:
            self._memo[key] = analyze_parametric_boundaries(
                self.formula, self.variables, tuple(parameters)
            )
        return self._memo[key]

    @property
    def variable_order(self):
        from .planner.features import extract_problem_features
        from .planner.heuristics import choose_best_variable_order

        if "variable_order" not in self._memo:
            features = extract_problem_features(self.parsed_formula, variables=self.variables)
            self._memo["variable_order"] = choose_best_variable_order(
                features, self.polynomials, equational_constraints=self.equational_constraints
            )
        return self._memo["variable_order"]

    def complete_cad(self):
        from .cad_algorithms.decomposition import decomp_collins_complete

        if "complete_cad" not in self._memo:
            _ = self.projection_tower
            with computation_context(self.computation):
                self._memo["complete_cad"] = decomp_collins_complete(
                    self.polynomials or (1,), self.variable_order
                )
        return self._memo["complete_cad"]

    def satisfiable(self, *, return_result: bool = False):
        from .decision.api import is_satisfiable

        with computation_context(self.computation):
            return is_satisfiable(self.formula, self.variables, return_result=return_result)

    def with_formula(self, formula):
        """Return a derived context sharing this context's exact-computation cache."""
        return SemialgebraicContext(formula, self.variables, self.computation)

    def add_constraints(self, *constraints):
        """Return a derived context for this region intersected with extra constraints."""
        import sympy as sp

        return self.with_formula(sp.And(self.formula, *constraints))

    def implies(self, conclusion, *, return_result=False):
        from .decision import implies

        with computation_context(self.computation):
            return implies(self.formula, conclusion, self.variables, return_result=return_result)

    def function_sign(self, expression, *, return_result=False):
        from .reasoning_signs import function_sign

        with computation_context(self.computation):
            return function_sign(
                expression, self.variables, assumptions=self.formula, return_result=return_result
            )

    def matrix_definiteness(
        self, matrix, *, requested="positive_semidefinite", return_result=False
    ):
        from .matrix_analysis import matrix_definiteness

        with computation_context(self.computation):
            return matrix_definiteness(
                matrix,
                self.variables,
                domain=self.formula,
                requested=requested,
                return_result=return_result,
            )

    def strict_feasible(self, *, relative=True, return_result=False):
        from .strict_feasibility import strict_feasible

        with computation_context(self.computation):
            return strict_feasible(
                self.formula, self.variables, relative=relative, return_result=return_result
            )

    def qe(self, quantifiers, *, free_variables=None, variable_order_strategy="auto"):
        from .qe.complete import qe_by_complete_cad

        with computation_context(self.computation):
            return qe_by_complete_cad(
                self.variables,
                tuple(quantifiers),
                self.parsed_formula,
                free_variables=free_variables,
                variable_order_strategy=variable_order_strategy,
            )

    def stats(self) -> dict[str, object]:
        incidence_components = len(self.incidence.components)
        return {
            "structural_cache_entries": len(self._memo),
            "exact_computation": self.computation.stats(),
            "incidence_components": incidence_components,
        }
