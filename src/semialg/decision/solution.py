from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass, field

import sympy as sp

from .sampling_helpers import _samples_from_cells, _samples_from_components


@dataclass(frozen=True)
class IntervalComponent:
    """Exact one-dimensional connected component of a semialgebraic set.

    Components are represented as closed/open intervals in one variable. A
    point component is encoded by ``lower == upper`` with both endpoints
    closed and has dimension 0.
    """

    variable: sp.Symbol
    lower: sp.Expr
    upper: sp.Expr
    lower_closed: bool
    upper_closed: bool

    @property
    def is_point(self) -> bool:
        return bool(self.lower == self.upper and self.lower_closed and self.upper_closed)

    @property
    def dimension(self) -> int:
        return 0 if self.is_point else 1

    @property
    def bounded(self) -> bool:
        return self.lower != -sp.oo and self.upper != sp.oo

    @property
    def closed(self) -> bool:
        return (self.lower == -sp.oo or self.lower_closed) and (
            self.upper == sp.oo or self.upper_closed
        )

    def as_formula(self) -> sp.Expr:
        v = self.variable
        if self.is_point:
            return sp.Eq(v, self.lower)
        pieces: list[sp.Expr] = []
        if self.lower != -sp.oo:
            pieces.append(v >= self.lower if self.lower_closed else v > self.lower)
        if self.upper != sp.oo:
            pieces.append(v <= self.upper if self.upper_closed else v < self.upper)
        return sp.And(*pieces) if pieces else sp.true

    def sample_point(self) -> sp.Expr:
        if self.is_point:
            return self.lower
        if self.lower == -sp.oo and self.upper == sp.oo:
            return sp.Integer(0)
        if self.lower == -sp.oo:
            return self.upper - 1
        if self.upper == sp.oo:
            return self.lower + 1
        return sp.simplify((self.lower + self.upper) / 2)


@dataclass(frozen=True)
class SatisfiabilityResult:
    """Structured result for a real satisfiability query."""

    satisfiable: bool
    formula: sp.Expr
    variables: tuple[sp.Symbol, ...]
    witness: Mapping[sp.Symbol, sp.Expr] | None = None
    method: str = "unknown"
    diagnostics: Mapping[str, object] = field(default_factory=dict)

    @property
    def status(self) -> str:
        return "sat" if self.satisfiable else "unsat"

    def __bool__(self) -> bool:
        return self.satisfiable


@dataclass(frozen=True)
class TautologyResult:
    """Structured result for a real validity/tautology query."""

    tautology: bool
    formula: sp.Expr
    variables: tuple[sp.Symbol, ...]
    counterexample: Mapping[sp.Symbol, sp.Expr] | None = None
    method: str = "unknown"
    diagnostics: Mapping[str, object] = field(default_factory=dict)

    @property
    def valid(self) -> bool:
        return self.tautology

    def __bool__(self) -> bool:
        return self.tautology


@dataclass(frozen=True)
class ImplicationResult:
    """Structured result for an implication query."""

    valid: bool
    premise: sp.Expr
    conclusion: sp.Expr
    variables: tuple[sp.Symbol, ...]
    counterexample: Mapping[sp.Symbol, sp.Expr] | None = None
    method: str = "unknown"
    diagnostics: Mapping[str, object] = field(default_factory=dict)

    def __bool__(self) -> bool:
        return self.valid


@dataclass(frozen=True)
class EquivalenceResult:
    """Structured result for semantic equivalence of two formulas."""

    equivalent: bool
    lhs: sp.Expr
    rhs: sp.Expr
    variables: tuple[sp.Symbol, ...]
    counterexample: Mapping[sp.Symbol, sp.Expr] | None = None
    failed_direction: str | None = None
    method: str = "unknown"
    diagnostics: Mapping[str, object] = field(default_factory=dict)

    def __bool__(self) -> bool:
        return self.equivalent


@dataclass(frozen=True)
class SemialgebraicSolution:
    """Structured solution summary for a semialgebraic constraint system.

    ``formula`` is the best currently available quantifier-free description of
    the solution set in ``variables``. Optional metadata fields are left as ``None`` or empty tuples when an
    analysis is unsupported, rather than guessed.
    """

    formula: sp.Expr
    variables: tuple[sp.Symbol, ...]
    samples: tuple[Mapping[sp.Symbol, sp.Expr], ...]
    satisfiable: bool
    method: str
    diagnostics: Mapping[str, object] = field(default_factory=dict)
    parameters: tuple[sp.Symbol, ...] = ()
    simplified_constraints: tuple[sp.Expr, ...] = ()
    parameter_conditions: sp.Expr | None = None
    parameter_decomposition: object | None = None
    dimension: int | sp.Expr | None = None
    bounded: bool | None = None
    closed: bool | None = None
    compact: bool | None = None
    components: tuple[object, ...] = ()
    cells: tuple[object, ...] = ()
    cylindrical_solution: object | None = None
    connectivity: object | None = None

    @property
    def sample(self) -> Mapping[sp.Symbol, sp.Expr] | None:
        return self.samples[0] if self.samples else None

    @property
    def empty(self) -> bool:
        return not self.satisfiable

    @property
    def nonempty(self) -> bool:
        return self.satisfiable

    def __bool__(self) -> bool:
        return self.satisfiable

    def as_formula(self, *, prefer: str = "auto", closed_cells: bool = False) -> sp.Expr:
        """Return the best available formula representation of the solution."""

        key = prefer.lower().replace("-", "_")
        if not self.satisfiable:
            return sp.false

        def from_components() -> sp.Expr | None:
            if not self.components:
                return None
            # Prefer interval components as reduced formulas. Higher-dimensional
            # connectivity components are graph/grouping objects; their formulas
            # may intentionally describe open CAD cells, so cells/cylindrical
            # source formulas are better for human-facing reduced output.
            if not all(hasattr(component, "variable") for component in self.components):
                return None
            forms: list[sp.Expr] = []
            for component in self.components:
                as_formula = getattr(component, "as_formula", None)
                if as_formula is None:
                    return None
                try:
                    forms.append(as_formula())
                except TypeError:
                    forms.append(as_formula)
                except Exception:
                    return None
            if not forms:
                return None
            return sp.Or(*forms) if len(forms) > 1 else forms[0]

        def from_cylindrical() -> sp.Expr | None:
            cyl = self.cylindrical_solution
            if cyl is None:
                return None
            as_formula = getattr(cyl, "as_formula", None)
            if as_formula is None:
                return None
            try:
                return as_formula(closed=closed_cells)
            except TypeError:
                try:
                    return as_formula()
                except Exception:
                    return None
            except Exception:
                return None

        def from_cells() -> sp.Expr | None:
            if not self.cells:
                return None
            forms: list[sp.Expr] = []
            for cell in self.cells:
                as_formula = getattr(cell, "as_formula", None)
                if as_formula is None:
                    return None
                try:
                    forms.append(as_formula(closed=closed_cells))
                except TypeError:
                    try:
                        forms.append(as_formula())
                    except Exception:
                        return None
                except Exception:
                    return None
            if not forms:
                return None
            return sp.Or(*forms) if len(forms) > 1 else forms[0]

        if key in {"component", "components"}:
            attempts = ("components", "formula")
        elif key in {"cell", "cells"}:
            attempts = ("cells", "cylindrical", "formula")
        elif key in {"cylindrical", "cad", "cad_cells"}:
            attempts = ("cylindrical", "cells", "formula")
        elif key in {"formula", "stored"}:
            attempts = ("formula",)
        elif key in {"auto", "reduced", "reduced_formula", "best"}:
            attempts = ("components", "cells", "cylindrical", "formula")
        else:
            raise ValueError(f"unsupported formula preference: {prefer!r}")

        for attempt in attempts:
            candidate = None
            if attempt == "components":
                candidate = from_components()
            elif attempt == "cylindrical":
                candidate = from_cylindrical()
            elif attempt == "cells":
                candidate = from_cells()
            elif attempt == "formula":
                candidate = self.formula
            if candidate is not None:
                return candidate
        return self.formula

    def as_piecewise(self) -> sp.Piecewise:
        """Return an indicator-style Piecewise representation."""

        return sp.Piecewise((sp.Integer(1), self.as_formula()), (sp.Integer(0), True))

    def as_components(self) -> tuple[object, ...]:
        return self.components

    def as_cells(self) -> tuple[object, ...]:
        return self.cells

    def as_cylindrical(self) -> object | None:
        return self.cylindrical_solution

    def sample_points(
        self, mode: str = "auto", count: int | None = None
    ) -> tuple[Mapping[sp.Symbol, sp.Expr], ...]:
        """Return samples from stored components or cells when available."""

        key = mode.lower().replace("-", "_")
        if key in {"component", "components"}:
            key = "per_component"
        if key in {"cell", "cells"}:
            key = "per_cell"
        if key == "per_component":
            samples = _samples_from_components(self.components)
            if not samples and self.connectivity is not None:
                samples = _samples_from_components(getattr(self.connectivity, "components", ()))
            if not samples:
                samples = self.samples
        elif key == "per_cell":
            samples = _samples_from_cells(self.cells)
            if not samples and self.cylindrical_solution is not None:
                samples = tuple(getattr(self.cylindrical_solution, "sample_points", ()))
            if not samples:
                samples = self.samples
        elif key == "auto":
            samples = (
                self.samples
                or _samples_from_components(self.components)
                or _samples_from_cells(self.cells)
            )
        else:
            raise ValueError(f"unsupported sample mode: {mode!r}")
        return samples if count is None else samples[:count]

    def to_dict(self) -> dict[str, object]:
        """Return a lightweight inspection summary."""

        return {
            "formula": sp.sstr(self.formula),
            "variables": tuple(sp.sstr(v) for v in self.variables),
            "parameters": tuple(sp.sstr(p) for p in self.parameters),
            "satisfiable": self.satisfiable,
            "sample_count": len(self.samples),
            "dimension": self.dimension,
            "bounded": self.bounded,
            "closed": self.closed,
            "compact": self.compact,
            "component_count": len(self.components),
            "cell_count": len(self.cells),
            "has_cylindrical_solution": self.cylindrical_solution is not None,
            "has_connectivity": self.connectivity is not None,
            "has_parameter_conditions": self.parameter_conditions is not None,
            "has_parameter_decomposition": self.parameter_decomposition is not None,
            "method": self.method,
            "diagnostics": dict(self.diagnostics),
        }

    def explain(self) -> dict[str, object]:
        """Return solver diagnostics in a stable, user-facing layout."""

        diagnostics = dict(self.diagnostics)
        return {
            "method": self.method,
            "satisfiable": self.satisfiable,
            "used_interval_decomposition": bool(
                diagnostics.get("used_interval_decomposition", False)
            ),
            "used_cad": bool(diagnostics.get("used_cad", False)),
            "used_qe": bool(diagnostics.get("used_qe", False)),
            "used_cylindrical_solution": bool(diagnostics.get("used_cylindrical_solution", False)),
            "used_parameter_decomposition": bool(
                diagnostics.get("used_parameter_decomposition", False)
            ),
            "used_domain_normalization": bool(diagnostics.get("domain_normalized", False)),
            "normalization_steps": tuple(diagnostics.get("normalization_steps", ())),
            "removed_redundant_constraints": tuple(
                diagnostics.get("removed_redundant_constraints", ())
            ),
            "unsupported_features": tuple(diagnostics.get("unsupported_features", ())),
            "raw": diagnostics,
        }

    def contains(self, point: Mapping[sp.Symbol, sp.Expr]) -> bool:
        """Return whether ``point`` satisfies this solution formula."""

        from ..instances.real_fallbacks import satisfies_formula

        try:
            return bool(sp.simplify(self.as_formula().subs(point)))
        except Exception:
            return satisfies_formula(self.as_formula(), point, strict=False)

    def is_subset_of(self, other: SemialgebraicSolution | sp.Expr) -> bool:
        """Return whether this solution set is contained in ``other``."""

        other_formula = (
            other.as_formula() if isinstance(other, SemialgebraicSolution) else sp.sympify(other)
        )
        from .api import _merge_variables, implies

        variables = _merge_variables(
            self.variables, getattr(other, "variables", ()), self.as_formula(), other_formula
        )
        return implies(self.as_formula(), other_formula, variables)

    def is_equal_to(self, other: SemialgebraicSolution | sp.Expr) -> bool:
        """Return whether this solution set equals ``other``."""

        other_formula = (
            other.as_formula() if isinstance(other, SemialgebraicSolution) else sp.sympify(other)
        )
        from .api import _merge_variables, equivalent

        variables = _merge_variables(
            self.variables, getattr(other, "variables", ()), self.as_formula(), other_formula
        )
        return equivalent(self.as_formula(), other_formula, variables)

    def is_disjoint_from(self, other: SemialgebraicSolution | sp.Expr) -> bool:
        """Return whether this solution set is disjoint from ``other``."""

        other_formula = (
            other.as_formula() if isinstance(other, SemialgebraicSolution) else sp.sympify(other)
        )
        from .api import _merge_variables, is_satisfiable

        variables = _merge_variables(
            self.variables, getattr(other, "variables", ()), self.as_formula(), other_formula
        )
        return not is_satisfiable(sp.And(self.as_formula(), other_formula), variables)

    def union(self, other: SemialgebraicSolution | sp.Expr) -> SemialgebraicSolution:
        """Return a solved object for the union with ``other``."""

        other_formula = (
            other.as_formula() if isinstance(other, SemialgebraicSolution) else sp.sympify(other)
        )
        from .api import _merge_variables, solve_semialgebraic

        variables = _merge_variables(
            self.variables, getattr(other, "variables", ()), self.as_formula(), other_formula
        )
        return solve_semialgebraic(sp.Or(self.as_formula(), other_formula), variables, count=0)

    def intersection(self, other: SemialgebraicSolution | sp.Expr) -> SemialgebraicSolution:
        """Return a solved object for the intersection with ``other``."""

        other_formula = (
            other.as_formula() if isinstance(other, SemialgebraicSolution) else sp.sympify(other)
        )
        from .api import _merge_variables, solve_semialgebraic

        variables = _merge_variables(
            self.variables, getattr(other, "variables", ()), self.as_formula(), other_formula
        )
        return solve_semialgebraic(sp.And(self.as_formula(), other_formula), variables, count=0)

    def difference(self, other: SemialgebraicSolution | sp.Expr) -> SemialgebraicSolution:
        """Return a solved object for this solution set minus ``other``."""

        other_formula = (
            other.as_formula() if isinstance(other, SemialgebraicSolution) else sp.sympify(other)
        )
        from .api import _merge_variables, solve_semialgebraic

        variables = _merge_variables(
            self.variables, getattr(other, "variables", ()), self.as_formula(), other_formula
        )
        return solve_semialgebraic(
            sp.And(self.as_formula(), sp.Not(other_formula)), variables, count=0
        )

    def complement(self) -> SemialgebraicSolution:
        """Return a solved object for the complement in the current ambient variables."""

        from .api import solve_semialgebraic

        return solve_semialgebraic(sp.Not(self.as_formula()), self.variables, count=0)

    def discretize(self, *, bounds=None, samples_per_curve: int = 33):
        """Return lightweight plotting data for supported 1D/2D solutions."""

        from ..solution_geometry import discretize_solution

        return discretize_solution(self, bounds=bounds, samples_per_curve=samples_per_curve)

    def plot(
        self,
        *,
        bounds=None,
        samples_per_curve: int = 33,
        ax=None,
        show: bool = False,
        **plot_kwargs,
    ):
        """Plot supported 1D/2D solutions using Matplotlib."""

        from ..solution_geometry import plot_solution

        return plot_solution(
            self,
            bounds=bounds,
            samples_per_curve=samples_per_curve,
            ax=ax,
            show=show,
            **plot_kwargs,
        )
