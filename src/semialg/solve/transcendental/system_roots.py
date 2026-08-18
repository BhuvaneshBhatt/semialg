from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass, field
from itertools import product

import sympy as sp


@dataclass(frozen=True)
class CertifiedPoint:
    coordinates: tuple[sp.Expr, ...]
    residual_norm: sp.Expr
    certified: bool = False


@dataclass(frozen=True)
class CompletenessCertificate:
    complete: bool
    reason: str
    supporting_method: str
    details: dict = field(default_factory=dict)


@dataclass(frozen=True)
class SystemRootFallbackResult:
    variables: tuple[sp.Symbol, ...]
    points: tuple[tuple[sp.Expr, ...], ...] = ()
    certified_points: tuple[CertifiedPoint, ...] = ()
    completeness_certificate: CompletenessCertificate = CompletenessCertificate(
        False, "no_certificate", "none"
    )
    complete: bool = False
    method: str = "system_roots_fallback"
    metadata: dict = field(default_factory=dict)


@dataclass(frozen=True)
class SearchBox:
    lower_bounds: tuple[float, ...]
    upper_bounds: tuple[float, ...]


def _equations_from_formula(formula: sp.Expr) -> tuple[sp.Expr, ...]:
    atoms = list(formula.args) if isinstance(formula, sp.And) else [formula]
    eqs = []
    for atom in atoms:
        if isinstance(atom, sp.Equality):
            eqs.append(sp.simplify(atom.lhs - atom.rhs))
    return tuple(eqs)


def _certify_point(
    equations: Sequence[sp.Expr],
    variables: Sequence[sp.Symbol],
    point: Sequence[sp.Expr],
    tol: float = 1e-8,
) -> CertifiedPoint:
    subst = dict(zip(variables, point, strict=True))
    residuals = []
    for eq in equations:
        try:
            residuals.append(abs(complex(sp.N(eq.subs(subst), 30))))
        except Exception:
            residuals.append(float("inf"))
    residual_norm = max(residuals) if residuals else 0.0
    return CertifiedPoint(
        tuple(sp.nsimplify(v) for v in point),
        sp.nsimplify(residual_norm),
        certified=bool(residual_norm <= tol),
    )


def _unique_points(points):
    out = []
    seen = set()
    for pt in points:
        key = tuple(sp.nsimplify(v) for v in pt)
        if key not in seen:
            seen.add(key)
            out.append(key)
    return tuple(out)


def orchestrate_trans_search(
    formula: sp.Expr,
    variables: Sequence[sp.Symbol],
    *,
    search_box: SearchBox | None = None,
    grid_points_per_axis: int = 3,
) -> SystemRootFallbackResult:
    variables = tuple(variables)
    eqs = _equations_from_formula(formula)
    if not eqs or len(variables) < 2:
        cert = CompletenessCertificate(False, "not_applicable", "orchestrator")
        return SystemRootFallbackResult(
            variables=variables,
            points=(),
            completeness_certificate=cert,
            complete=False,
            method="system_roots_not_applicable",
        )

    try:
        direct = sp.nonlinsolve(eqs, variables)
        if isinstance(direct, sp.FiniteSet) and direct:
            pts = tuple(tuple(sp.simplify(v) for v in pt) for pt in direct)
            certs = tuple(_certify_point(eqs, variables, pt) for pt in pts)
            cert = CompletenessCertificate(
                True, "finite_symbolic_solution", "nonlinsolve", {"point_count": len(pts)}
            )
            return SystemRootFallbackResult(
                variables=variables,
                points=pts,
                certified_points=certs,
                completeness_certificate=cert,
                complete=True,
                method="nonlinsolve",
            )
    except Exception:
        pass

    if search_box is None:
        search_box = SearchBox(
            lower_bounds=tuple([-2.0] * len(variables)), upper_bounds=tuple([2.0] * len(variables))
        )

    seeds_by_axis = []
    for lo, hi in zip(search_box.lower_bounds, search_box.upper_bounds, strict=True):
        if grid_points_per_axis == 1:
            seeds_by_axis.append([(lo + hi) / 2.0])
        else:
            step = (hi - lo) / (grid_points_per_axis - 1)
            seeds_by_axis.append([lo + k * step for k in range(grid_points_per_axis)])

    found = []
    attempts = 0
    for seed in product(*seeds_by_axis):
        attempts += 1
        try:
            root = sp.nsolve(eqs, variables, seed, tol=1e-14, maxsteps=100, prec=50)
            if getattr(root, "shape", None) is not None:
                pt = tuple(sp.nsimplify(root[i]) for i in range(len(variables)))
            else:
                pt = (sp.nsimplify(root),)
            found.append(pt)
        except Exception:
            continue

    pts = _unique_points(found)
    certs = tuple(_certify_point(eqs, variables, pt) for pt in pts)
    complete = False
    cert = CompletenessCertificate(
        False,
        "numerical_grid_search_only",
        "nsolve_grid_search",
        {
            "attempt_count": attempts,
            "point_count": len(pts),
            "certified_count": sum(int(cp.certified) for cp in certs),
        },
    )
    return SystemRootFallbackResult(
        variables=variables,
        points=pts,
        certified_points=certs,
        completeness_certificate=cert,
        complete=complete,
        method="nsolve_grid_search",
        metadata={
            "search_box": search_box,
            "grid_points_per_axis": grid_points_per_axis,
        },
    )


def solve_bounded_trans_sys(
    formula: sp.Expr,
    variables: Sequence[sp.Symbol],
    *,
    search_box: SearchBox | None = None,
    grid_points_per_axis: int = 3,
) -> SystemRootFallbackResult:
    return orchestrate_trans_search(
        formula,
        variables,
        search_box=search_box,
        grid_points_per_axis=grid_points_per_axis,
    )


__all__ = [
    "CertifiedPoint",
    "CompletenessCertificate",
    "SystemRootFallbackResult",
    "SearchBox",
    "orchestrate_trans_search",
    "solve_bounded_trans_sys",
]
