from __future__ import annotations

from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass, field, replace
from typing import Any

import sympy as sp
from sympy.logic.boolalg import Boolean

FormulaLike = sp.Expr | Boolean | bool


def _normalize_condition(condition: FormulaLike) -> sp.Expr:
    if condition is True:
        return sp.true
    if condition is False:
        return sp.false
    expr = condition if isinstance(condition, (sp.Basic, Boolean)) else sp.sympify(condition)
    if expr is sp.true or expr == sp.true:
        return sp.true
    if expr is sp.false or expr == sp.false:
        return sp.false
    try:
        return sp.simplify_logic(sp.simplify(expr), form="dnf")
    except (TypeError, ValueError, NotImplementedError, AttributeError):
        return sp.simplify(expr)


def _normalize_parameters(
    parameters: Sequence[sp.Symbol | str],
    *,
    known_symbols: Iterable[sp.Symbol] = (),
) -> tuple[sp.Symbol, ...]:
    """Normalize parameters without inventing assumption-incompatible symbols.

    String names are resolved against symbols already present in the guarded
    expressions whenever possible.  This matters because ``Symbol("a")`` and
    ``Symbol("a", real=True)`` are distinct SymPy objects.
    """

    by_name: dict[str, list[sp.Symbol]] = {}
    for symbol in known_symbols:
        by_name.setdefault(symbol.name, []).append(symbol)

    out: list[sp.Symbol] = []
    seen: set[sp.Symbol] = set()
    for item in parameters:
        if isinstance(item, str):
            matches = tuple(dict.fromkeys(by_name.get(item, ())))
            if len(matches) > 1:
                raise ValueError(
                    f"parameter name {item!r} is ambiguous across symbols with different assumptions"
                )
            symbol = matches[0] if matches else sp.Symbol(item, real=True)
        else:
            symbol = item
        if not isinstance(symbol, sp.Symbol):
            raise TypeError("parameters must be SymPy symbols or symbol names")
        if symbol not in seen:
            out.append(symbol)
            seen.add(symbol)
    return tuple(out)


def _resolve_assignment_keys(
    assignments: Mapping[sp.Symbol | str, object],
    symbols: Iterable[sp.Symbol],
) -> dict[sp.Symbol, object]:
    """Resolve string assignment keys to the actual symbols in an API object."""

    symbols = tuple(dict.fromkeys(symbols))
    by_name: dict[str, list[sp.Symbol]] = {}
    for symbol in symbols:
        by_name.setdefault(symbol.name, []).append(symbol)

    resolved: dict[sp.Symbol, object] = {}
    for key, value in assignments.items():
        if isinstance(key, str):
            matches = tuple(dict.fromkeys(by_name.get(key, ())))
            if not matches:
                raise ValueError(f"unknown parameter name {key!r}")
            if len(matches) > 1:
                raise ValueError(
                    f"parameter name {key!r} is ambiguous across symbols with different assumptions"
                )
            symbol = matches[0]
        elif isinstance(key, sp.Symbol):
            symbol = key
        else:
            raise TypeError("assignment keys must be SymPy symbols or symbol names")
        resolved[symbol] = value
    return resolved


def _substitute_truth(
    condition: sp.Expr, assignments: Mapping[sp.Symbol | str, object]
) -> bool | None:
    substitutions = _resolve_assignment_keys(assignments, condition.free_symbols)
    reduced = sp.simplify(condition.subs(substitutions))
    if reduced is sp.true or reduced == sp.true:
        return True
    if reduced is sp.false or reduced == sp.false:
        return False
    try:
        if reduced == True:  # noqa: E712
            return True
        if reduced == False:  # noqa: E712
            return False
    except TypeError:
        pass
    return None


def _values_equal(left: object, right: object) -> bool:
    if left is right:
        return True
    try:
        eq = left == right
        if isinstance(eq, bool):
            return eq
    except Exception:
        pass
    if isinstance(left, sp.Basic) or isinstance(right, sp.Basic):
        try:
            return sp.simplify(sp.sympify(left) - sp.sympify(right)) == 0
        except (TypeError, ValueError, NotImplementedError, AttributeError):
            return False
    return False


@dataclass(frozen=True)
class ConditionalBranch:
    """One value guarded by a semialgebraic parameter condition.

    ``condition`` describes exactly where ``value`` is valid.  ``sample`` is
    optional supporting evidence and is never used as a substitute for the
    guard itself.  ``certified`` records whether the branch validity has been
    established by an exact semialgebraic argument.
    """

    condition: sp.Expr
    value: Any
    sample: Mapping[sp.Symbol, sp.Expr] = field(default_factory=dict)
    certified: bool = True
    metadata: Mapping[str, object] = field(default_factory=dict)

    def __post_init__(self) -> None:
        object.__setattr__(self, "condition", _normalize_condition(self.condition))

    @property
    def empty(self) -> bool:
        return self.condition is sp.false or self.condition == sp.false

    def applies(self, assignments: Mapping[sp.Symbol | str, object]) -> bool | None:
        """Return whether this branch applies after substituting parameters.

        ``None`` means the supplied assignments are insufficient to decide the
        guard exactly.
        """

        return _substitute_truth(self.condition, assignments)

    def specialize(self, assignments: Mapping[sp.Symbol | str, object]) -> ConditionalBranch:
        """Substitute parameter values into both the guard and symbolic value."""

        symbols = set(self.condition.free_symbols)
        if isinstance(self.value, sp.Basic):
            symbols.update(self.value.free_symbols)
        substitutions = _resolve_assignment_keys(assignments, symbols)
        value = self.value
        if hasattr(value, "subs"):
            try:
                value = value.subs(substitutions)
            except (TypeError, ValueError, NotImplementedError, AttributeError):
                pass
        return ConditionalBranch(
            _normalize_condition(self.condition.subs(substitutions)),
            value,
            self.sample,
            self.certified,
            self.metadata,
        )


@dataclass(frozen=True)
class ParameterStratifiedResult:
    """A first-class piecewise result over semialgebraic parameter strata.

    The branches are intended to be pairwise disjoint when ``disjoint`` is
    true. ``coverage_condition`` is the portion of parameter space covered by
    the result. ``complete`` means that this coverage is the full requested
    parameter domain; it is deliberately separate from per-branch
    certification.
    """

    parameters: tuple[sp.Symbol, ...]
    branches: tuple[ConditionalBranch, ...]
    coverage_condition: sp.Expr = sp.true
    complete: bool = True
    disjoint: bool = True
    certified: bool = True
    method: str = "parameter_stratification"
    diagnostics: Mapping[str, object] = field(default_factory=dict)

    def __post_init__(self) -> None:
        branches = tuple(branch for branch in self.branches if not branch.empty)
        known_symbols: set[sp.Symbol] = set(self.coverage_condition.free_symbols)
        for branch in branches:
            known_symbols.update(branch.condition.free_symbols)
            if isinstance(branch.value, sp.Basic):
                known_symbols.update(branch.value.free_symbols)
        object.__setattr__(
            self, "parameters", _normalize_parameters(self.parameters, known_symbols=known_symbols)
        )
        object.__setattr__(
            self, "coverage_condition", _normalize_condition(self.coverage_condition)
        )
        object.__setattr__(self, "branches", branches)
        if self.certified and any(not branch.certified for branch in self.branches):
            object.__setattr__(self, "certified", False)

    @property
    def stratum_count(self) -> int:
        return len(self.branches)

    @property
    def conditions(self) -> tuple[sp.Expr, ...]:
        return tuple(branch.condition for branch in self.branches)

    @property
    def values(self) -> tuple[Any, ...]:
        return tuple(branch.value for branch in self.branches)

    def normalize(self, *, merge_equal_values: bool = True) -> ParameterStratifiedResult:
        """Simplify guards and optionally merge branches with equal values."""

        branches = [
            ConditionalBranch(
                branch.condition, branch.value, branch.sample, branch.certified, branch.metadata
            )
            for branch in self.branches
            if branch.condition is not sp.false and branch.condition != sp.false
        ]
        if not merge_equal_values:
            return replace(self, branches=tuple(branches))

        merged: list[ConditionalBranch] = []
        for branch in branches:
            for index, current in enumerate(merged):
                if _values_equal(branch.value, current.value):
                    merged[index] = ConditionalBranch(
                        _normalize_condition(sp.Or(current.condition, branch.condition)),
                        current.value,
                        current.sample or branch.sample,
                        current.certified and branch.certified,
                        {**dict(current.metadata), **dict(branch.metadata)},
                    )
                    break
            else:
                merged.append(branch)
        return replace(self, branches=tuple(merged))

    def select(
        self,
        assignments: Mapping[sp.Symbol | str, object],
        *,
        require_unique: bool = True,
    ) -> Any:
        """Return the value selected by exact parameter substitution.

        Raises ``ValueError`` when no branch is known to apply, when the
        assignments are insufficient to decide the active stratum, or when
        multiple branches apply despite ``require_unique=True``.
        """

        matches: list[ConditionalBranch] = []
        undecided = False
        for branch in self.branches:
            applies = branch.applies(assignments)
            if applies is True:
                matches.append(branch)
            elif applies is None:
                undecided = True
        if not matches:
            if undecided:
                raise ValueError("parameter assignments are insufficient to select a stratum")
            raise ValueError("parameter assignments lie outside the covered strata")
        if require_unique and len(matches) != 1:
            raise ValueError("parameter assignments select multiple strata")
        selected = matches[0].specialize(assignments)
        return selected.value

    def condition_for_value(self, value: object) -> sp.Expr:
        """Return the union of guards on which the result equals ``value``."""

        pieces = [
            branch.condition for branch in self.branches if _values_equal(branch.value, value)
        ]
        if not pieces:
            return sp.false
        return _normalize_condition(sp.Or(*pieces))

    def as_piecewise(self, *, default: object = sp.nan) -> sp.Piecewise:
        """Convert symbolic-valued branches to a SymPy ``Piecewise`` object."""

        pieces = [(sp.sympify(branch.value), branch.condition) for branch in self.branches]
        pieces.append((sp.sympify(default), True))
        return sp.Piecewise(*pieces)

    def specialize(
        self, assignments: Mapping[sp.Symbol | str, object]
    ) -> ParameterStratifiedResult:
        """Substitute a partial parameter assignment into all branches."""

        substitutions = _resolve_assignment_keys(assignments, self.parameters)
        remaining = tuple(
            parameter for parameter in self.parameters if parameter not in substitutions
        )
        branches = tuple(branch.specialize(assignments) for branch in self.branches)
        coverage = _normalize_condition(self.coverage_condition.subs(substitutions))
        return ParameterStratifiedResult(
            remaining,
            branches,
            coverage,
            self.complete,
            self.disjoint,
            self.certified,
            self.method,
            self.diagnostics,
        ).normalize()


def conditional_result(
    parameters: Sequence[sp.Symbol | str],
    branches: Iterable[ConditionalBranch | tuple[FormulaLike, object]],
    *,
    coverage_condition: FormulaLike = sp.true,
    complete: bool = True,
    disjoint: bool = True,
    certified: bool = True,
    method: str = "conditional_result",
    diagnostics: Mapping[str, object] | None = None,
    normalize: bool = True,
) -> ParameterStratifiedResult:
    """Construct a normalized first-class parameter-stratified result."""

    converted: list[ConditionalBranch] = []
    for branch in branches:
        if isinstance(branch, ConditionalBranch):
            converted.append(branch)
        else:
            condition, value = branch
            converted.append(ConditionalBranch(_normalize_condition(condition), value))
    result = ParameterStratifiedResult(
        _normalize_parameters(
            parameters,
            known_symbols={
                symbol for branch in converted for symbol in branch.condition.free_symbols
            }
            | {
                symbol
                for branch in converted
                if isinstance(branch.value, sp.Basic)
                for symbol in branch.value.free_symbols
            },
        ),
        tuple(converted),
        _normalize_condition(coverage_condition),
        complete,
        disjoint,
        certified,
        method,
        diagnostics or {},
    )
    return result.normalize() if normalize else result


@dataclass(frozen=True)
class ParameterStratificationCertificate:
    """Exact verification summary for a parameter-stratified result."""

    pairwise_disjoint: bool
    coverage_verified: bool
    branches_certified: bool
    verified: bool
    uncovered_condition: sp.Expr = sp.false
    overlap_conditions: tuple[sp.Expr, ...] = ()
    diagnostics: Mapping[str, object] = field(default_factory=dict)

    def verify(self) -> bool:
        return self.verified


def verify_parameter_stratification(
    result: ParameterStratifiedResult,
    *,
    parameter_domain: FormulaLike | None = None,
) -> ParameterStratificationCertificate:
    """Verify disjointness and coverage of parameter branches exactly.

    The checks use semialgebraic satisfiability only when Boolean
    simplification cannot settle a question directly.  Coverage is checked
    against ``parameter_domain`` when supplied and otherwise against the
    result's own ``coverage_condition``.
    """

    from .decision import is_satisfiable

    domain = (
        result.coverage_condition
        if parameter_domain is None
        else _normalize_condition(parameter_domain)
    )
    overlaps: list[sp.Expr] = []
    pairwise = True
    for i, left in enumerate(result.branches):
        for right in result.branches[i + 1 :]:
            overlap = _normalize_condition(sp.And(domain, left.condition, right.condition))
            if overlap is sp.false or overlap == sp.false:
                continue
            try:
                satisfiable = bool(is_satisfiable(overlap, result.parameters))
            except (TypeError, ValueError, NotImplementedError, ArithmeticError):
                satisfiable = True
            if satisfiable:
                pairwise = False
                overlaps.append(overlap)

    union = (
        _normalize_condition(sp.Or(*(branch.condition for branch in result.branches)))
        if result.branches
        else sp.false
    )
    uncovered = _normalize_condition(sp.And(domain, sp.Not(union)))
    if uncovered is sp.false or uncovered == sp.false:
        coverage = True
    else:
        try:
            coverage = not bool(is_satisfiable(uncovered, result.parameters))
        except (TypeError, ValueError, NotImplementedError, ArithmeticError):
            coverage = False

    branches_certified = all(branch.certified for branch in result.branches)
    verified = pairwise and coverage and branches_certified
    return ParameterStratificationCertificate(
        pairwise,
        coverage,
        branches_certified,
        verified,
        sp.false if coverage else uncovered,
        tuple(overlaps),
        {"branch_count": len(result.branches), "parameter_count": len(result.parameters)},
    )


__all__ = [
    "ConditionalBranch",
    "ParameterStratifiedResult",
    "ParameterStratificationCertificate",
    "conditional_result",
    "verify_parameter_stratification",
]
