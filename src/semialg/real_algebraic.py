"""Certified real feasibility for rational polynomial equality systems.

The positive-dimensional backend follows the critical/polar reduction of
Aubry, Rouillier and Safey El Din.  Certified equidimensional pieces are
reduced recursively to zero-dimensional systems; the RUR solver then constructs exact real
witnesses.  Failure of a genericity search is reported
as incomplete rather than as infeasibility.
"""

from __future__ import annotations

from collections.abc import Iterable, Sequence
from dataclasses import dataclass
from itertools import combinations

import sympy as sp

from .algebraic.equality_ideal import EqualityIdealContext
from .algebraic.groebner_utils import lex_basis_generators
from .algebraic.rational_univariate import RationalUnivariateError
from .algebraic_decomposition import equidimensional_decomposition, verify_decomposition_certificate
from .solve.zero_dimensional import solve_zero_dimensional_system


@dataclass(frozen=True)
class RealAlgebraicFeasibilityResult:
    variables: tuple[sp.Symbol, ...]
    satisfiable: bool | None
    witness: tuple[tuple[sp.Symbol, sp.Expr], ...] | None
    complete: bool
    method: str
    zero_dimensional_systems: tuple[tuple[sp.Expr, ...], ...] = ()
    notes: tuple[str, ...] = ()

    @property
    def assignment(self) -> dict[sp.Symbol, sp.Expr] | None:
        return dict(self.witness) if self.witness is not None else None

    def require_decision(self) -> bool:
        if not self.complete or self.satisfiable is None:
            raise NotImplementedError("real-algebraic feasibility certification is incomplete")
        return self.satisfiable


def _normalize(
    equations: Iterable[sp.Expr | sp.Equality], variables: tuple[sp.Symbol, ...]
) -> tuple[sp.Expr, ...]:
    out: list[sp.Expr] = []
    for equation in equations:
        expr = sp.sympify(equation)
        if isinstance(expr, sp.Equality):
            expr = expr.lhs - expr.rhs
        expr = sp.expand(expr)
        if expr != 0:
            sp.Poly(expr, *variables, domain=sp.QQ)
            out.append(expr)
    return tuple(out)


def _extract_triangular(
    basis: Sequence[sp.Expr], variables: tuple[sp.Symbol, ...]
) -> tuple[tuple[sp.Expr, ...], tuple[sp.Symbol, ...], tuple[sp.Symbol, ...]]:
    remaining = list(basis)
    triangular: list[sp.Expr] = []
    high: list[sp.Symbol] = []
    low: list[sp.Symbol] = []
    for variable in reversed(variables):
        active = [poly for poly in remaining if poly.has(variable)]
        remaining = [poly for poly in remaining if not poly.has(variable)]
        if active:
            chosen = min(active, key=lambda poly: int(sp.degree(poly, variable)))
            high.append(variable)
            triangular.append(chosen)
        else:
            low.append(variable)
    return tuple(reversed(triangular)), tuple(reversed(high)), tuple(reversed(low))


def _gamma_a(
    triangular: Sequence[sp.Expr],
    point: Sequence[int],
    high: Sequence[sp.Symbol],
    low: Sequence[sp.Symbol],
) -> tuple[sp.Expr, ...]:
    codim = len(high)
    if codim == 0 or not triangular:
        return tuple()
    t = tuple(triangular)
    high_matrix = sp.Matrix([[sp.diff(poly, x) for x in high] for poly in t])
    augmented_high = high_matrix.col_join(
        sp.Matrix([[x - a for x, a in zip(high, point[-codim:], strict=True)]])
    )
    minors = []
    for rows in combinations(range(augmented_high.rows), codim):
        minors.append(sp.expand(augmented_high.extract(rows, range(codim)).det()))
    minors = list(reversed(minors))
    low_matrix = sp.Matrix([[sp.diff(poly, x) for x in low] for poly in t])
    augmented_low = low_matrix.col_join(
        sp.Matrix([[x - a for x, a in zip(low, point[:-codim], strict=True)]])
    )
    signs = tuple((-1) ** i for i in range(codim + 1))
    gamma: list[sp.Expr] = []
    for column in range(augmented_low.cols):
        values = [augmented_low[row, column] for row in range(augmented_low.rows)]
        if len(values) != len(minors):
            continue
        gamma.append(
            sp.expand(
                sum(
                    s * value * minor for s, value, minor in zip(signs, values, minors, strict=True)
                )
            )
        )
    return tuple(g for g in gamma if g != 0)


def _solve_real_zero_dimensional(
    system: Sequence[sp.Expr],
    variables: tuple[sp.Symbol, ...],
    *,
    solver=solve_zero_dimensional_system,
) -> tuple[tuple[tuple[sp.Expr, ...], ...] | None, str | None]:
    """Run the exact zero-dimensional backend with explicit incomplete semantics."""
    try:
        points = solver(system, vars=variables, real=True)
    except RationalUnivariateError:
        return None, "rur_failure"
    return tuple(points), None


def _generic_points(variable_count: int, attempts: int):
    point = [0] * variable_count
    yield tuple(point)
    for k in range(attempts - 1):
        point[k % variable_count] += 1
        yield tuple(point)


def _polar_zero_dimensional_systems(
    equations: tuple[sp.Expr, ...],
    variables: tuple[sp.Symbol, ...],
    *,
    max_generic_attempts: int,
    max_pieces: int,
    depth: int = 0,
) -> tuple[tuple[tuple[sp.Expr, ...], ...], bool, tuple[str, ...]]:
    if depth > len(variables) + 2:
        return tuple(), False, ("polar_recursion_limit",)
    context = EqualityIdealContext.build(equations, variables)
    if context.inconsistent:
        return tuple(), True, tuple()
    if context.zero_dimensional:
        return (tuple(context.generators),), True, tuple()

    decomposition = equidimensional_decomposition(equations, variables, max_pieces=max_pieces)
    if (
        not decomposition.complete
        or decomposition.certificate is None
        or not verify_decomposition_certificate(decomposition.certificate)
    ):
        return tuple(), False, ("uncertified_equidimensional_decomposition",)

    terminals: list[tuple[sp.Expr, ...]] = []
    notes: list[str] = []
    for piece in decomposition.pieces:
        piece_context = EqualityIdealContext.build(piece.equations, variables)
        if piece_context.zero_dimensional:
            terminals.append(tuple(piece_context.generators))
            continue
        basis = lex_basis_generators(piece_context.generators, variables, domain=sp.QQ)
        triangular, high, low = _extract_triangular(basis, variables)
        if not high:
            return tuple(terminals), False, (*notes, "triangular_extraction_failed")
        current_dim = piece_context.dimension
        reduced: tuple[sp.Expr, ...] | None = None
        for a in _generic_points(len(variables), max_generic_attempts):
            gamma = _gamma_a(triangular, a, high, low)
            if not gamma:
                continue
            candidate = EqualityIdealContext.build((*piece_context.generators, *gamma), variables)
            if candidate.inconsistent:
                continue
            if candidate.dimension < current_dim:
                reduced = tuple(candidate.generators)
                notes.append(f"polar_dimension_drop:{current_dim}->{candidate.dimension}")
                break
        if reduced is None:
            return tuple(terminals), False, (*notes, "generic_polar_search_exhausted")
        sub, complete, subnotes = _polar_zero_dimensional_systems(
            reduced,
            variables,
            max_generic_attempts=max_generic_attempts,
            max_pieces=max_pieces,
            depth=depth + 1,
        )
        terminals.extend(sub)
        notes.extend(subnotes)
        if not complete:
            return tuple(terminals), False, tuple(notes)
    return tuple(terminals), True, tuple(notes)


def real_algebraic_feasibility(
    equations: Iterable[sp.Expr | sp.Equality],
    variables: Sequence[sp.Symbol],
    *,
    max_generic_attempts: int = 64,
    max_pieces: int = 128,
) -> RealAlgebraicFeasibilityResult:
    """Decide real feasibility of rational polynomial equations when ARS certifies it."""
    vars_ = tuple(variables)
    if not vars_:
        raise ValueError("a nonempty variable list is required")
    if max_generic_attempts < 1:
        raise ValueError("max_generic_attempts must be positive")
    if max_pieces < 1:
        raise ValueError("max_pieces must be positive")
    try:
        eqs = _normalize(equations, vars_)
    except (sp.PolynomialError, ValueError, TypeError) as exc:
        raise RationalUnivariateError(
            "ARS feasibility requires rational polynomial equations"
        ) from exc
    if not eqs:
        witness = tuple((v, sp.Integer(0)) for v in vars_)
        return RealAlgebraicFeasibilityResult(
            vars_, True, witness, True, "ars", notes=("empty_equality_system",)
        )
    context = EqualityIdealContext.build(eqs, vars_)
    if context.inconsistent:
        return RealAlgebraicFeasibilityResult(
            vars_, False, None, True, "ars", notes=("unit_ideal",)
        )
    if context.zero_dimensional:
        systems = (tuple(context.generators),)
        points, failure = _solve_real_zero_dimensional(eqs, vars_)
        if failure is not None or points is None:
            return RealAlgebraicFeasibilityResult(
                vars_, None, None, False, "ars", systems, ("rur_zero_dimensional_failure",)
            )
        witness = tuple(zip(vars_, points[0], strict=True)) if points else None
        return RealAlgebraicFeasibilityResult(
            vars_, bool(points), witness, True, "ars+rur", systems
        )

    systems, complete, notes = _polar_zero_dimensional_systems(
        eqs,
        vars_,
        max_generic_attempts=max_generic_attempts,
        max_pieces=max_pieces,
    )
    for system in systems:
        points, failure = _solve_real_zero_dimensional(system, vars_)
        if failure is not None or points is None:
            return RealAlgebraicFeasibilityResult(
                vars_, None, None, False, "ars", systems, (*notes, "rur_terminal_failure")
            )
        for point in points:
            assignment = dict(zip(vars_, point, strict=True))
            if all(sp.simplify(eq.subs(assignment)) == 0 for eq in eqs):
                witness = tuple((v, sp.simplify(assignment[v])) for v in vars_)
                return RealAlgebraicFeasibilityResult(
                    vars_, True, witness, True, "ars+rur", systems, notes
                )
    if complete:
        return RealAlgebraicFeasibilityResult(vars_, False, None, True, "ars+rur", systems, notes)
    return RealAlgebraicFeasibilityResult(vars_, None, None, False, "ars", systems, notes)


def solve_real_algebraic_set(
    equations: Iterable[sp.Expr | sp.Equality],
    variables: Sequence[sp.Symbol],
    **kwargs,
) -> dict[sp.Symbol, sp.Expr] | None:
    """Return one exact real point, ``None`` for certified emptiness, or raise if incomplete."""
    result = real_algebraic_feasibility(equations, variables, **kwargs)
    decision = result.require_decision()
    return result.assignment if decision else None


__all__ = [
    "RealAlgebraicFeasibilityResult",
    "real_algebraic_feasibility",
    "solve_real_algebraic_set",
]
