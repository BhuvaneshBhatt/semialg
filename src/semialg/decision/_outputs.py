from __future__ import annotations

from collections.abc import Mapping, Sequence

import sympy as sp

from ._inputs import normalize_symbols
from .solution import SemialgebraicSolution


class CellOutput(tuple):
    """Tuple-compatible cells view with a ``.cells`` alias."""

    @property
    def cells(self):
        return tuple(self)


def select_solution_output(result: SemialgebraicSolution, output: str | None) -> object:
    if output is None:
        return result
    key = output.lower().replace("-", "_")
    if key in {"result", "solution", "object"}:
        return result
    if key in {"formula", "constraints", "reduced_formula", "reduced", "best_formula", "reduce"}:
        return result.formula
    if key in {"piecewise", "indicator", "indicator_piecewise"}:
        return result.as_piecewise()
    if key in {"sample", "one_sample"}:
        return result.sample
    if key in {"samples", "points"}:
        return result.samples
    if key in {"components", "component"}:
        return result.components
    if key in {"cells", "cell"}:
        return CellOutput(result.cells)
    if key in {"cylindrical", "cylindrical_solution", "cylindrical_cells", "cad_cells"}:
        return result.cylindrical_solution
    if key in {"connectivity", "adjacency", "roadmap", "roadmap_graph", "components_graph"}:
        return result.connectivity
    if key in {"plot_data", "discretization", "discretized", "mesh_data"}:
        return result.discretize()
    if key in {"diagnostics", "explain", "explanation"}:
        return result.explain()
    if key in {
        "parameter_strata",
        "parameter_decomposition",
        "strata",
        "piecewise_solution",
        "parameterized_solution",
    }:
        return result.parameter_decomposition
    if key in {"conditions", "parameter_conditions", "solvability_conditions"}:
        if result.parameter_conditions is not None:
            return result.parameter_conditions
        if not result.parameters:
            return sp.true if result.satisfiable else sp.false
        return sp.false
    if key in {"satisfiable", "nonempty"}:
        return result.satisfiable
    if key in {"empty", "unsatisfiable"}:
        return result.empty
    raise ValueError(f"unsupported solve_semialgebraic output selector: {output!r}")


def add_standard_solver_diagnostics(
    diagnostics: dict[str, object],
    *,
    method: str,
    variables: Sequence[sp.Symbol],
    projection_order: Sequence[sp.Symbol | str] | None,
    domain_normalization: object | None,
    metadata: Mapping[str, object] | None = None,
    parameter_decomposition: object | None = None,
    solved: object | None = None,
) -> dict[str, object]:
    metadata = metadata or {}
    diagnostics.setdefault("backend", method)
    diagnostics.setdefault("normalization_steps", ())
    diagnostics.setdefault("removed_redundant_constraints", ())
    diagnostics.setdefault("unsupported_features", ())
    diagnostics["requested_method"] = method
    diagnostics["variable_order"] = tuple(sp.sstr(v) for v in variables)
    diagnostics["projection_order"] = tuple(sp.sstr(v) for v in normalize_symbols(projection_order))
    diagnostics["used_interval_decomposition"] = (
        bool(metadata.get("components")) and len(tuple(variables)) == 1
    )
    diagnostics["used_cad"] = bool(
        metadata.get("cells")
        or metadata.get("cylindrical_solution")
        or metadata.get("connectivity")
    ) or method in {"cad", "qe", "cylindrical"}
    diagnostics["used_qe"] = method in {"cad", "qe"} or solved is not None
    diagnostics["used_cylindrical_solution"] = metadata.get("cylindrical_solution") is not None
    diagnostics["used_connectivity"] = metadata.get("connectivity") is not None
    diagnostics["used_parameter_decomposition"] = parameter_decomposition is not None
    rewrites = (
        tuple(getattr(domain_normalization, "rewrites", ()))
        if domain_normalization is not None
        else ()
    )
    constraints = (
        tuple(sp.sstr(c) for c in getattr(domain_normalization, "domain_constraints", ()))
        if domain_normalization is not None
        else ()
    )
    active_domain_norm = bool(rewrites or constraints)
    diagnostics["domain_normalized"] = active_domain_norm
    if domain_normalization is not None:
        diagnostics["domain_rewrites"] = rewrites
        diagnostics["domain_constraints"] = constraints
        if active_domain_norm:
            diagnostics["normalization_steps"] = tuple(
                diagnostics.get("normalization_steps", ())
            ) + ("domain-sensitive-constraint-normalization",)
    return diagnostics


__all__ = ["CellOutput", "add_standard_solver_diagnostics", "select_solution_output"]
