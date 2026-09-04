"""Exact convexity certificates for semialgebraic sets and polynomials."""

from __future__ import annotations

from collections.abc import Callable, Mapping, Sequence
from dataclasses import dataclass, field
from functools import lru_cache

import sympy as sp

from .algebraic.groebner_utils import compute_groebner_basis
from .cad_algorithms.cells import CylindricalSolution, extract_cylindrical_solution
from .cad_algorithms.selected_samples import extract_selected_cad_samples
from .connectivity import (
    CADConnectivityGraph,
    build_cad_adjacency_graph,
    extract_cad_connectivity,
    factorized_equality_components,
)
from .context import with_computation_context
from .decomposition import cad
from .errors import CertificationFailure
from .formula import parse_formula
from .internal_symbols import fresh_real_dummy
from .matrix_analysis import constant_symmetric_inertia, matrix_psd_on, principal_minors
from .normalization import normalize_formula, normalize_problem_variables
from .qe import qe_by_complete_cad
from .sampling import sample_point, sample_points


@dataclass(frozen=True)
class ConvexityCertificate:
    """Exact certificate or exact inconclusive result for a set-convexity query."""

    outcome: bool | None
    method: str
    formula: sp.Expr
    variables: tuple[sp.Symbol, ...]
    witness: Mapping[sp.Symbol, sp.Expr] | None = None
    details: Mapping[str, object] = field(default_factory=dict)

    @property
    def certified(self) -> bool:
        """Return whether the result contains a decisive exact answer."""

        return self.outcome is not None


@dataclass(frozen=True)
class PolynomialConvexityCertificate:
    """Exact certificate for the requested polynomial Hessian sign on a domain.

    A positive certificate is a sufficient convexity/concavity proof whenever
    the domain used by the caller is already known to be convex.  A negative
    Hessian-sign result is not, by itself, a representation-independent
    nonconvexity theorem for a lower-dimensional or otherwise constrained set.
    """

    outcome: bool | None
    expression: sp.Expr
    variables: tuple[sp.Symbol, ...]
    domain: sp.Expr
    sense: str
    method: str
    principal_minors: tuple[sp.Expr, ...] = ()
    hessian_inertia: tuple[int, int, int] | None = None

    @property
    def certified(self) -> bool:
        """Return whether the requested Hessian sign was decided exactly."""

        return self.outcome is not None


@dataclass(frozen=True)
class QuadraticConvexityCertificate:
    """Exact structural certificate for a basic quadratic semialgebraic set."""

    outcome: bool | None
    formula: sp.Expr
    variables: tuple[sp.Symbol, ...]
    method: str
    atom_certificates: tuple[PolynomialConvexityCertificate, ...] = ()

    @property
    def certified(self) -> bool:
        """Return whether the quadratic structure gives a decisive answer."""

        return self.outcome is not None


@dataclass
class _ConvexityCADAnalysis:
    """Lazily shared certified CAD products for one convexity query."""

    formula: sp.Expr
    variables: tuple[sp.Symbol, ...]
    solution: CylindricalSolution | None = None
    connectivity: CADConnectivityGraph | None = None
    extractor: Callable[..., CylindricalSolution] = extract_cylindrical_solution

    def cylindrical_solution(self) -> CylindricalSolution:
        if self.solution is None:
            self.solution = self.extractor(self.formula, self.variables, selected_only=True)
        return self.solution

    def connectivity_graph(self) -> CADConnectivityGraph:
        if self.connectivity is None:
            self.connectivity = build_cad_adjacency_graph(
                self.cylindrical_solution(), formula=self.formula
            )
        return self.connectivity


def _exact_constant_sign(value: sp.Expr) -> int | None:
    """Return the exact sign of a constant real expression when decidable."""

    value = sp.simplify(value)
    if value == 0:
        return 0
    if value.is_positive is True:
        return 1
    if value.is_negative is True:
        return -1
    try:
        sign = sp.sign(value)
    except (TypeError, ValueError, NotImplementedError):
        return None
    if sign in (-1, 0, 1):
        return int(sign)
    return None


@dataclass(frozen=True)
class _AnalyzedAtom:
    atom: sp.Expr
    residual: sp.Expr
    poly: sp.Poly | None
    degree: int | None
    relation_kind: str


def _relation_kind(atom: sp.Expr) -> str:
    if isinstance(atom, sp.Equality):
        return "eq"
    if isinstance(atom, sp.Unequality):
        return "ne"
    if isinstance(atom, (sp.LessThan, sp.StrictLessThan)):
        return "le"
    if isinstance(atom, (sp.GreaterThan, sp.StrictGreaterThan)):
        return "ge"
    return "other"


def _analyze_atoms(
    formula: sp.Expr, variables: tuple[sp.Symbol, ...]
) -> tuple[_AnalyzedAtom, ...] | None:
    atoms = tuple(formula.args) if isinstance(formula, sp.And) else (formula,)
    if not all(getattr(atom, "is_Relational", False) for atom in atoms):
        return None
    result: list[_AnalyzedAtom] = []
    for atom in atoms:
        residual = sp.expand(atom.lhs - atom.rhs)
        try:
            poly = sp.Poly(residual, *variables)
            degree = int(poly.total_degree())
        except (sp.PolynomialError, ValueError, TypeError):
            poly = None
            degree = None
        result.append(_AnalyzedAtom(atom, residual, poly, degree, _relation_kind(atom)))
    return tuple(result)


def _atoms(formula: sp.Expr) -> tuple[sp.Expr, ...] | None:
    if isinstance(formula, sp.And):
        atoms = tuple(formula.args)
    else:
        atoms = (formula,)
    if all(getattr(atom, "is_Relational", False) for atom in atoms):
        return atoms
    return None


@lru_cache(maxsize=128)
def _hessian_data(
    expression: sp.Expr, variables: tuple[sp.Symbol, ...], sense: str
) -> tuple[sp.ImmutableMatrix, tuple[sp.Expr, ...], tuple[int, int, int] | None]:
    hessian = sp.ImmutableMatrix(sp.hessian(expression, variables))
    matrix = hessian if sense == "convex" else -hessian
    if all(not entry.free_symbols for entry in matrix):
        return sp.ImmutableMatrix(matrix), (), constant_symmetric_inertia(matrix)
    return sp.ImmutableMatrix(matrix), principal_minors(matrix), None


@lru_cache(maxsize=128)
def _cached_polynomial_convexity_certificate(
    expr: sp.Expr, vars_: tuple[sp.Symbol, ...], dom: sp.Expr, sense: str
) -> PolynomialConvexityCertificate:
    try:
        sp.Poly(expr, *vars_)
    except (sp.PolynomialError, ValueError, TypeError):
        return PolynomialConvexityCertificate(None, expr, vars_, dom, sense, "nonpolynomial")
    if not vars_:
        return PolynomialConvexityCertificate(
            True, expr, vars_, dom, sense, "constant-polynomial", ()
        )
    if dom is sp.false or dom == sp.false:
        return PolynomialConvexityCertificate(
            True, expr, vars_, dom, sense, "empty-domain", (), None
        )
    hessian_matrix, minors, inertia = _hessian_data(expr, vars_, sense)
    if inertia is not None:
        positive, negative, zero = inertia
        return PolynomialConvexityCertificate(
            negative == 0,
            expr,
            vars_,
            dom,
            sense,
            "constant-hessian-ldlt-inertia",
            (),
            (positive, negative, zero),
        )
    definiteness = matrix_psd_on(hessian_matrix, vars_, domain=dom, return_result=True)
    if definiteness.outcome is False:
        return PolynomialConvexityCertificate(
            False, expr, vars_, dom, sense, "principal-minor-counterexample", minors
        )
    if definiteness.outcome is True:
        return PolynomialConvexityCertificate(
            True, expr, vars_, dom, sense, "principal-minor-shared-cad", minors
        )
    return PolynomialConvexityCertificate(
        None, expr, vars_, dom, sense, "matrix-definiteness-inconclusive", minors
    )


def polynomial_convexity_certificate(
    expression,
    variables: Sequence[sp.Symbol | str] | None = None,
    domain=sp.true,
    *,
    sense: str = "convex",
) -> PolynomialConvexityCertificate:
    """Certify a polynomial Hessian sign globally or on a semialgebraic domain.

    ``sense='convex'`` proves positive-semidefiniteness of the Hessian;
    ``sense='concave'`` proves positive-semidefiniteness of its negative.  The
    proof uses the exact principal-minor criterion for symmetric matrices and
    complete CAD only for minor signs that SymPy cannot settle directly.
    """

    expr = sp.sympify(expression)
    dom = normalize_formula(domain)
    vars_ = normalize_problem_variables(variables, expr, dom)
    if sense not in {"convex", "concave"}:
        raise ValueError("sense must be 'convex' or 'concave'")
    return _cached_polynomial_convexity_certificate(expr, vars_, dom, sense)


def _is_affine_polyhedron(
    formula: sp.Expr,
    variables: tuple[sp.Symbol, ...],
    analyzed: tuple[_AnalyzedAtom, ...] | None = None,
) -> bool:
    """Return whether ``formula`` is an intersection of affine relations."""

    analyzed = analyzed if analyzed is not None else _analyze_atoms(formula, variables)
    if analyzed is None:
        return False
    for item in analyzed:
        if item.relation_kind == "ne" or item.degree is None or item.degree > 1:
            return False
    return True


def quadratic_convexity_certificate(
    region,
    variables: Sequence[sp.Symbol | str] | None = None,
    *,
    _analyzed: tuple[_AnalyzedAtom, ...] | None = None,
) -> QuadraticConvexityCertificate:
    """Certify convexity of an intersection of affine/quadratic relations.

    Equalities must be affine.  A quadratic sublevel must have a PSD constant
    Hessian and a quadratic superlevel an NSD constant Hessian.  This classifier
    is exact for the structures it accepts and otherwise returns inconclusive.
    """

    formula = normalize_formula(region)
    vars_ = normalize_problem_variables(variables, formula)
    analyzed = _analyzed if _analyzed is not None else _analyze_atoms(formula, vars_)
    if analyzed is None:
        return QuadraticConvexityCertificate(None, formula, vars_, "not-basic-quadratic")
    certificates: list[PolynomialConvexityCertificate] = []
    for item in analyzed:
        atom = item.atom
        residual = item.residual
        poly = item.poly
        if item.relation_kind == "ne":
            return QuadraticConvexityCertificate(None, formula, vars_, "unequality")
        if poly is None or item.degree is None:
            return QuadraticConvexityCertificate(None, formula, vars_, "nonpolynomial")
        if item.degree > 2:
            return QuadraticConvexityCertificate(None, formula, vars_, "degree-above-two")
        if isinstance(atom, sp.Equality):
            if item.degree > 1:
                return QuadraticConvexityCertificate(None, formula, vars_, "nonaffine-equality")
            continue
        sense = "convex" if isinstance(atom, (sp.LessThan, sp.StrictLessThan)) else "concave"
        cert = polynomial_convexity_certificate(residual, vars_, sense=sense)
        certificates.append(cert)
        if cert.outcome is not True:
            return QuadraticConvexityCertificate(
                None, formula, vars_, "indefinite-quadratic", tuple(certificates)
            )
    return QuadraticConvexityCertificate(
        True, formula, vars_, "quadratic-intersection", tuple(certificates)
    )


def _is_convex_polynomial_intersection(
    formula: sp.Expr,
    variables: tuple[sp.Symbol, ...],
    analyzed: tuple[_AnalyzedAtom, ...] | None = None,
) -> bool:
    """Certify a basic polynomial intersection using global Hessian signs."""

    analyzed = analyzed if analyzed is not None else _analyze_atoms(formula, variables)
    if analyzed is None:
        return False
    for item in analyzed:
        if item.relation_kind == "ne" or item.poly is None or item.degree is None:
            return False
        if item.relation_kind == "eq":
            if item.degree > 1:
                return False
            continue
        sense = "convex" if item.relation_kind == "le" else "concave"
        cert = polynomial_convexity_certificate(item.residual, variables, sense=sense)
        if cert.outcome is not True:
            return False
    return True


def _domain_relative_polynomial_certificate(
    formula: sp.Expr,
    variables: tuple[sp.Symbol, ...],
    analyzed: tuple[_AnalyzedAtom, ...] | None = None,
) -> ConvexityCertificate | None:
    """Certify a basic set by adding constraints over an already-convex domain.

    The growing domain starts with affine constraints.  A nonlinear sublevel is
    admitted only when its Hessian is PSD throughout the already-certified
    convex domain; superlevels use NSD Hessians.  This avoids circularly assuming
    convexity of the final feasible set while proving domain-relative curvature.
    """

    analyzed = analyzed if analyzed is not None else _analyze_atoms(formula, variables)
    if analyzed is None or any(item.relation_kind == "ne" for item in analyzed):
        return None
    affine: list[sp.Expr] = []
    pending: list[_AnalyzedAtom] = []
    for item in analyzed:
        if item.degree is None:
            return None
        if item.relation_kind == "eq":
            if item.degree > 1:
                return None
            affine.append(item.atom)
        elif item.degree <= 1:
            affine.append(item.atom)
        else:
            pending.append(item)
    if not pending:
        return None
    domain = sp.And(*affine) if affine else sp.true
    added: list[sp.Expr] = []
    remaining = list(pending)
    while remaining:
        progress = False
        remaining.sort(key=lambda item: (item.degree or 10**9, len(item.residual.free_symbols)))
        for item in tuple(remaining):
            atom = item.atom
            sense = "convex" if item.relation_kind == "le" else "concave"
            cert = polynomial_convexity_certificate(
                item.residual, variables, domain=domain, sense=sense
            )
            if cert.outcome is True:
                added.append(atom)
                domain = sp.And(domain, atom)
                remaining.remove(item)
                progress = True
        if not progress:
            return None
    return ConvexityCertificate(
        True,
        "domain-relative-hessian",
        formula,
        variables,
        details={"certified_constraints": tuple(added)},
    )


@lru_cache(maxsize=64)
def _is_convex_by_definition(formula: sp.Expr, variables: tuple[sp.Symbol, ...]) -> bool:
    """Decide convexity from the complete quantified segment definition."""

    xs = tuple(fresh_real_dummy(f"{v.name}_a") for v in variables)
    ys = tuple(fresh_real_dummy(f"{v.name}_b") for v in variables)
    t = fresh_real_dummy("convex_t")
    sx = formula.xreplace(dict(zip(variables, xs, strict=True)))
    sy = formula.xreplace(dict(zip(variables, ys, strict=True)))
    mixed = tuple(t * x + (1 - t) * y for x, y in zip(xs, ys, strict=True))
    smix = formula.xreplace(dict(zip(variables, mixed, strict=True)))
    counterexample = sp.to_nnf(sp.And(sx, sy, t >= 0, t <= 1, sp.Not(smix)), simplify=False)
    qvars = (*xs, *ys, t)
    qe = qe_by_complete_cad(
        qvars,
        tuple(("exists", v) for v in qvars),
        parse_formula(counterexample),
        return_result=True,
    )
    return qe.truth_value is False


def _truth_at(formula: sp.Expr, point: Mapping[sp.Symbol, sp.Expr]) -> bool | None:
    if formula is sp.true or formula == sp.true:
        return True
    if formula is sp.false or formula == sp.false:
        return False
    if isinstance(formula, sp.And):
        values = [_truth_at(arg, point) for arg in formula.args]
        if False in values:
            return False
        return True if all(value is True for value in values) else None
    if isinstance(formula, sp.Or):
        values = [_truth_at(arg, point) for arg in formula.args]
        if True in values:
            return True
        return False if all(value is False for value in values) else None
    if isinstance(formula, sp.Not):
        value = _truth_at(formula.args[0], point)
        return None if value is None else not value
    if getattr(formula, "is_Relational", False):
        try:
            value = sp.cancel((formula.lhs - formula.rhs).subs(point))
            sign = _exact_constant_sign(value)
        except (ArithmeticError, TypeError, ValueError, NotImplementedError):
            sign = None
        if sign is not None:
            if isinstance(formula, sp.Equality):
                return sign == 0
            if isinstance(formula, sp.Unequality):
                return sign != 0
            if isinstance(formula, (sp.LessThan, sp.StrictLessThan)):
                return sign <= 0 if isinstance(formula, sp.LessThan) else sign < 0
            if isinstance(formula, (sp.GreaterThan, sp.StrictGreaterThan)):
                return sign >= 0 if isinstance(formula, sp.GreaterThan) else sign > 0
    value = sp.simplify(formula.subs(point))
    if value is sp.true or value == sp.true:
        return True
    if value is sp.false or value == sp.false:
        return False
    return None


def _plausibly_zero_dimensional(
    formula: sp.Expr, variables: tuple[sp.Symbol, ...], analyzed: tuple[_AnalyzedAtom, ...] | None
) -> bool:
    if analyzed is None:
        return False
    equalities = [item for item in analyzed if item.relation_kind == "eq" and item.poly is not None]
    if not equalities:
        return False
    # Definite quadratic equal to its unique extremal value, e.g. x**2+y**2=0.
    if len(equalities) == 1 and equalities[0].degree == 2:
        _matrix, _minors, inertia = _hessian_data(equalities[0].residual, variables, "convex")
        if inertia is not None and (inertia[0] == len(variables) or inertia[1] == len(variables)):
            try:
                gradient = [sp.diff(equalities[0].residual, v) for v in variables]
                critical = sp.solve(gradient, variables, dict=True)
                if (
                    len(critical) == 1
                    and sp.simplify(equalities[0].residual.subs(critical[0])) == 0
                ):
                    return True
            except (TypeError, ValueError, NotImplementedError):
                pass
    try:
        gb = compute_groebner_basis((item.residual for item in equalities), variables)
        return bool(gb.is_zero_dimensional)
    except (sp.PolynomialError, TypeError, ValueError, NotImplementedError):
        return len(equalities) >= len(variables)


def _singleton_certificate(
    formula: sp.Expr,
    variables: tuple[sp.Symbol, ...],
    analyzed: tuple[_AnalyzedAtom, ...] | None = None,
) -> ConvexityCertificate | None:
    """Recognize singleton sets exactly when equality structure suggests one."""

    analyzed = analyzed if analyzed is not None else _analyze_atoms(formula, variables)
    if not _plausibly_zero_dimensional(formula, variables, analyzed):
        return None
    # Affine intersections are convex regardless of whether they collapse to a
    # singleton. Affine recognition is cheaper than the complete
    # sample-and-uniqueness query, so handle those sets separately.
    if _is_affine_polyhedron(formula, variables, analyzed):
        return None
    point = sample_point(formula, variables, strategy="complete", exact=True)
    if point is None:
        return ConvexityCertificate(True, "empty-set", formula, variables)
    distinct = sp.Or(*(sp.Ne(variable, point[variable]) for variable in variables))
    other = sp.And(formula, distinct)
    qe = qe_by_complete_cad(
        variables,
        tuple(("exists", variable) for variable in variables),
        parse_formula(other),
        return_result=True,
    )
    if qe.truth_value is False:
        return ConvexityCertificate(True, "singleton", formula, variables, witness=point)
    return None


def _one_dimensional_certificate(
    formula: sp.Expr, variables: tuple[sp.Symbol, ...], strategy: str | None
) -> ConvexityCertificate:
    """Decide one-dimensional convexity by exact ordered CAD-cell contiguity.

    In one real dimension a set is convex iff it is an interval (possibly
    empty, a singleton, open/closed/half-open, or unbounded).  General CAD
    connectivity can identify the two sides of a removed point through closure
    adjacency, so it is not the right criterion for interval convexity.  The
    ordered truth cells give the exact and cheaper test directly: no false CAD
    cell may lie between two selected cells.
    """

    del strategy
    decomposition = cad(formula, variables, output="cells", return_result=True)
    ordered_cells = tuple(sorted(decomposition.cad.cells_by_level[1], key=lambda cell: cell.index))
    selected = {cell.index for cell in decomposition.cells}
    selected_positions = [i for i, cell in enumerate(ordered_cells) if cell.index in selected]
    if not selected_positions:
        return ConvexityCertificate(True, "one-dimensional-empty", formula, variables)
    first, last = selected_positions[0], selected_positions[-1]
    gap_positions = tuple(
        i for i in range(first, last + 1) if ordered_cells[i].index not in selected
    )
    outcome = not gap_positions
    component_count = 1 + sum(
        1
        for left, right in zip(selected_positions, selected_positions[1:], strict=False)
        if right > left + 1
    )
    details = {
        "selected_cell_count": len(selected_positions),
        "gap_count": len(gap_positions),
        "component_count": component_count,
    }
    if gap_positions:
        details["gap_cells"] = tuple(ordered_cells[i].index for i in gap_positions)
    return ConvexityCertificate(
        outcome,
        "one-dimensional-ordered-cad-interval",
        formula,
        variables,
        details=details,
    )


def _topology_negative_certificate(
    formula: sp.Expr,
    variables: tuple[sp.Symbol, ...],
    analysis: _ConvexityCADAnalysis | None = None,
    *,
    factorized_only: bool = False,
) -> ConvexityCertificate | None:
    """Reject cases with a cheaply certifiable disconnected topology.

    Full connectivity reconstruction is intentionally reserved for algebraic
    equality sets here; for general Boolean regions it can cost as much as the
    final convexity query.  The later exact segment search and quantified
    fallback remain complete.
    """

    factors = factorized_equality_components(formula, variables)
    if factors is not None and len(factors) > 1:
        return ConvexityCertificate(
            False,
            "disconnected-factorized-variety",
            formula,
            variables,
            details={"component_count": len(factors)},
        )
    if factorized_only:
        return None
    atoms = _atoms(formula)
    if atoms is None or not atoms or not all(isinstance(atom, sp.Equality) for atom in atoms):
        return None
    graph = (
        analysis.connectivity_graph()
        if analysis is not None
        else extract_cad_connectivity(formula, variables)
    )
    if graph.component_count > 1:
        return ConvexityCertificate(
            False,
            "disconnected-cad",
            formula,
            variables,
            details={"component_count": graph.component_count},
        )
    return None


def _ordered_sample_pairs(cells, points):
    """Yield promising midpoint-witness pairs first, then all remaining pairs."""

    n = len(points)
    if n < 2:
        return
    emitted: set[tuple[int, int]] = set()

    def emit(i: int, j: int):
        pair = (min(i, j), max(i, j))
        if pair[0] != pair[1] and pair not in emitted:
            emitted.add(pair)
            return pair
        return None

    # Different coarse CAD branches/components are high-yield candidates.
    if cells and len(cells) == n:
        for i in range(n):
            for j in range(i + 1, n):
                left = cells[i].index
                right = cells[j].index
                common = 0
                for a, b in zip(left, right, strict=False):
                    if a != b:
                        break
                    common += 1
                if common <= max(0, len(left) - 2):
                    pair = emit(i, j)
                    if pair is not None:
                        yield pair

    # Extremal representatives often expose holes/nonconvex arcs immediately.
    extremal = tuple(dict.fromkeys((0, n - 1, n // 2, n // 4, (3 * n) // 4)))
    for i in extremal:
        for j in extremal:
            pair = emit(i, j)
            if pair is not None:
                yield pair

    # Complete fallback: ordering changes performance only, never semantics.
    for i in range(n):
        for j in range(i + 1, n):
            pair = emit(i, j)
            if pair is not None:
                yield pair


def _segment_counterexample_certificate(
    formula: sp.Expr,
    variables: tuple[sp.Symbol, ...],
    *,
    sample_limit: int = 32,
    analysis: _ConvexityCADAnalysis | None = None,
) -> ConvexityCertificate | None:
    """Search exact CAD samples for a verified midpoint convexity violation."""

    del analysis
    samples = extract_selected_cad_samples(formula, variables, limit=sample_limit)
    cells = list(samples)
    points = [sample.point for sample in samples]
    if len(points) < 2:
        points = sample_points(
            formula,
            variables,
            count=sample_limit,
            strategy="cad_cells",
            exact=True,
        )
    for left_index, right_index in _ordered_sample_pairs(cells, points):
        left, right = points[left_index], points[right_index]
        midpoint = {
            variable: sp.cancel((left[variable] + right[variable]) / 2) for variable in variables
        }
        if _truth_at(formula, midpoint) is False:
            witness = dict(midpoint)
            details = {"left": left, "right": right, "t": sp.Rational(1, 2)}
            return ConvexityCertificate(
                False,
                "exact-segment-counterexample",
                formula,
                variables,
                witness=witness,
                details=details,
            )
    return None


@with_computation_context
def convexity_certificate(
    region,
    variables: Sequence[sp.Symbol | str] | None = None,
    *,
    strategy: str | None = None,
) -> ConvexityCertificate:
    """Decide semialgebraic set convexity through an exact staged hierarchy.

    The stages are formula normalization, trivial/empty/singleton handling,
    complete one-dimensional classification, affine/polyhedral recognition,
    quadratic recognition, global polynomial Hessian certificates,
    domain-relative Hessian certificates, topology-based rejection, exact
    segment-counterexample search, and finally complete quantified QE.
    """

    formula = normalize_formula(region)
    vars_ = normalize_problem_variables(variables, formula)
    if formula is sp.false or formula == sp.false:
        return ConvexityCertificate(True, "empty-set", formula, vars_)
    if formula is sp.true or formula == sp.true or not vars_:
        return ConvexityCertificate(True, "trivial", formula, vars_)

    if len(vars_) == 1:
        return _one_dimensional_certificate(formula, vars_, strategy)

    analyzed = _analyze_atoms(formula, vars_)

    # Cheap reducible-variety/topological rejection comes before any expensive
    # nonlinear singleton uniqueness query.
    topology = _topology_negative_certificate(formula, vars_, None, factorized_only=True)
    if topology is not None:
        return topology

    singleton = _singleton_certificate(formula, vars_, analyzed)
    if singleton is not None:
        return singleton

    if _is_affine_polyhedron(formula, vars_, analyzed):
        return ConvexityCertificate(True, "affine-polyhedron", formula, vars_)

    quadratic = quadratic_convexity_certificate(formula, vars_, _analyzed=analyzed)
    if quadratic.outcome is True:
        return ConvexityCertificate(
            True,
            "quadratic-intersection",
            formula,
            vars_,
            details={"quadratic_certificate": quadratic},
        )

    if _is_convex_polynomial_intersection(formula, vars_, analyzed):
        return ConvexityCertificate(True, "global-polynomial-hessian", formula, vars_)

    relative = _domain_relative_polynomial_certificate(formula, vars_, analyzed)
    if relative is not None:
        return relative

    analysis = _ConvexityCADAnalysis(formula, vars_)
    topology = _topology_negative_certificate(formula, vars_, analysis)
    if topology is not None:
        return topology

    witness = _segment_counterexample_certificate(formula, vars_, analysis=analysis)
    if witness is not None:
        return witness

    outcome = _is_convex_by_definition(formula, vars_)
    return ConvexityCertificate(outcome, "quantified-segment-definition", formula, vars_)


def is_convex(
    region,
    variables: Sequence[sp.Symbol | str] | None = None,
    *,
    strategy: str | None = None,
) -> bool:
    """Return whether a semialgebraic region is convex, exactly."""

    certificate = convexity_certificate(region, variables, strategy=strategy)
    if certificate.outcome is None:
        raise CertificationFailure("convexity hierarchy could not certify a final result")
    return certificate.outcome


__all__ = [
    "ConvexityCertificate",
    "PolynomialConvexityCertificate",
    "QuadraticConvexityCertificate",
    "convexity_certificate",
    "polynomial_convexity_certificate",
    "quadratic_convexity_certificate",
    "is_convex",
]
