"""Exact set-reduced decomposition of polynomial loci by dimension."""

from __future__ import annotations

from collections.abc import Iterable, Sequence
from itertools import combinations

import sympy as sp

from ._algebraic_decomposition_types import (
    AlgebraicLocusPiece,
    AssociatedPrimesResult,
    DecompositionCertificate,
    DecompositionPieceCertificate,
    EquidimensionalDecomposition,
    MinimalPrimeDecompositionCertificate,
    PrimaryComponent,
    PrimaryDecompositionCertificate,
    PrimaryDecompositionResult,
    RadicalComponent,
    RadicalIdealCertificate,
    RadicalIdealResult,
    RadicalMinimalPrimeDecomposition,
    RegularChainComponent,
    RegularChainDecomposition,
    RegularChainDecompositionCertificate,
    RegularChainSplitCertificate,
    TriangularPrimalityCertificate,
    TriangularPrimalityStage,
    _certificate_digest,
    _TriangularPrimalityAnalysis,
)
from .algebraic.equality_ideal import EqualityIdealContext
from .algebraic.groebner_utils import lex_basis_generators
from .algebraic.hilbert import ideal_degree
from .algebraic.rational_univariate.representation import RationalUnivariateError
from .algebraic.subresultants import subresultant_prs
from .algebraic_function_fields import (
    FunctionFieldError,
    MonogenicFunctionField,
    RationalFunctionField,
    certified_factor_univariate,
    field_element_to_expr,
    tower_evaluate,
)


def _normalize_equations(
    equations: Iterable[sp.Expr | sp.Equality], variables: Sequence[sp.Symbol]
) -> tuple[sp.Expr, ...]:
    raw: list[sp.Expr] = []
    for equation in equations:
        expr = sp.sympify(equation)
        if isinstance(expr, sp.Equality):
            expr = expr.lhs - expr.rhs
        expr = sp.expand(expr)
        if expr != 0:
            raw.append(expr)
    if not raw:
        return tuple()
    context = EqualityIdealContext.build(raw, variables)
    if context.inconsistent:
        return (sp.Integer(1),)
    if context.groebner_basis is None:
        return tuple()
    return tuple(sp.expand(poly.as_expr()) for poly in context.groebner_basis.polys)


def _canonical_factor(factor: sp.Expr, variables: Sequence[sp.Symbol]) -> sp.Expr:
    """Return a deterministic unit-normalized polynomial factor."""
    factor = sp.expand(factor)
    try:
        poly = sp.Poly(factor, *variables, extension=True)
        factor = sp.expand(poly.monic().as_expr())
    except (sp.PolynomialError, ValueError, TypeError, NotImplementedError):
        pass
    if factor.could_extract_minus_sign():
        factor = -factor
    return sp.expand(factor)


def _squarefree_factor_atoms(
    expression: sp.Expr, variables: Sequence[sp.Symbol]
) -> tuple[sp.Expr, ...] | None:
    """Return distinct exact nonconstant factors of a polynomial.

    The factors are used only through the identity ``V(product) = union V(f)``;
    no irreducibility or primality claim is required.  ``None`` means exact
    factorization was unavailable.
    """
    try:
        _unit, factors = sp.factor_list(expression, *variables, extension=True)
    except (sp.PolynomialError, ValueError, TypeError, NotImplementedError):
        return None
    atoms: list[sp.Expr] = []
    for factor, _multiplicity in factors:
        atom = _canonical_factor(factor, variables)
        if not atom.free_symbols.intersection(variables):
            continue
        if atom not in atoms:
            atoms.append(atom)
    return tuple(atoms)


def _distinct_factors(
    expression: sp.Expr, variables: Sequence[sp.Symbol]
) -> tuple[sp.Expr, ...] | None:
    """Return a nontrivial set-reducing factor split, if one is available."""
    atoms = _squarefree_factor_atoms(expression, variables)
    if atoms is None:
        return None
    # Detect multiplicity reduction separately because one repeated factor is
    # still a useful reduced-set normalization.
    try:
        _unit, factors = sp.factor_list(expression, *variables, extension=True)
    except (sp.PolynomialError, ValueError, TypeError, NotImplementedError):
        return None
    multiplicity_reduction = any(
        multiplicity > 1 and sp.sympify(factor).free_symbols.intersection(variables)
        for factor, multiplicity in factors
    )
    if len(atoms) > 1 or (len(atoms) == 1 and multiplicity_reduction):
        return atoms
    return None


def _minimal_vertex_covers(
    supports: Sequence[frozenset[int]], variable_count: int
) -> tuple[tuple[int, ...], ...]:
    """Return the minimal vertex covers of a finite support hypergraph.

    For a monomial ideal, these covers index the minimal coordinate primes of
    its radical.  Keeping only inclusion-minimal covers is exactly the
    reduced-set operation that discards embedded monomial primes.
    """
    edges = tuple(
        support
        for support in supports
        if support and not any(other < support for other in supports)
    )
    if not edges:
        return (tuple(),)

    covers: list[tuple[int, ...]] = []
    indices = tuple(range(variable_count))
    for size in range(variable_count + 1):
        for candidate in combinations(indices, size):
            candidate_set = frozenset(candidate)
            if any(set(existing).issubset(candidate_set) for existing in covers):
                continue
            if all(candidate_set.intersection(edge) for edge in edges):
                covers.append(candidate)
    return tuple(covers)


def _monomial_minimal_prime_generators(
    context: EqualityIdealContext,
) -> tuple[tuple[sp.Expr, ...], ...] | None:
    """Return the exact minimal-prime cover for a monomial ideal when available.

    If the reduced Groebner basis is monomial, the radical is a squarefree
    monomial ideal.  Its minimal primes are coordinate primes indexed by the
    minimal vertex covers of the generator supports.  This gives an exact
    reduced algebraic *and real* zero-set decomposition without requiring a
    general primary-decomposition backend.
    """
    basis = context.groebner_basis
    if basis is None or context.inconsistent:
        return None

    supports: list[frozenset[int]] = []
    for poly in basis.polys:
        terms = poly.terms()
        if len(terms) != 1:
            return None
        exponent, _coefficient = terms[0]
        support = frozenset(index for index, power in enumerate(exponent) if power)
        if not support:
            return tuple()
        supports.append(support)

    covers = _minimal_vertex_covers(tuple(supports), len(context.variables))
    return tuple(tuple(context.variables[index] for index in cover) for cover in covers)


def _factor_incidence_component_generators(
    context: EqualityIdealContext,
) -> tuple[tuple[sp.Expr, ...], ...] | None:
    """Return an exact minimal factor-incidence cover of ``V(context)``.

    For Groebner generators ``g_i = product_j f_ij`` we have, over both the
    real and complex numbers, ``V(g_i) = union_j V(f_ij)``. Distributing the
    intersection over these unions identifies component candidates with
    hitting sets of the factor supports. Inclusion-minimal hitting sets are
    sufficient under reduced-set semantics; supersets define loci already
    contained in a smaller candidate and therefore behave like embedded or
    redundant branches.

    This is an exact set-theoretic decomposition, not a claim that the
    candidate ideals are primary or prime. Candidates are recursively passed
    through the ordinary decomposition/certification pipeline.
    """
    basis = context.groebner_basis
    if basis is None or context.inconsistent or not basis.polys:
        return None

    atoms: list[sp.Expr] = []
    supports: list[frozenset[int]] = []
    has_nontrivial_product = False
    for poly in basis.polys:
        factors = _squarefree_factor_atoms(poly.as_expr(), context.variables)
        if factors is None or not factors:
            return None
        has_nontrivial_product |= len(factors) > 1
        support: set[int] = set()
        for factor in factors:
            try:
                index = atoms.index(factor)
            except ValueError:
                atoms.append(factor)
                index = len(atoms) - 1
            support.add(index)
        supports.append(frozenset(support))

    # A basis consisting only of one-factor generators would reproduce the
    # same ideal and must not recurse forever.
    if not has_nontrivial_product:
        return None

    covers = _minimal_vertex_covers(tuple(supports), len(atoms))
    return tuple(tuple(atoms[index] for index in cover) for cover in covers)


def _normalized_context_generators(context: EqualityIdealContext) -> tuple[sp.Expr, ...]:
    """Return deterministic Groebner generators for an equality context."""
    if context.groebner_basis is None:
        return tuple()
    return tuple(sp.expand(poly.as_expr()) for poly in context.groebner_basis.polys)


def _saturate_generators(
    generators: Sequence[sp.Expr],
    splitter: sp.Expr,
    variables: Sequence[sp.Symbol],
) -> tuple[sp.Expr, ...] | None:
    """Return generators of ``<generators> : splitter**infinity`` exactly.

    The Rabinowitsch identity
    ``I : h**infinity = (I + <1-u*h>) intersect k[x]`` is evaluated with a
    lexicographic Groebner basis and the fresh elimination variable first.
    ``None`` means the exact Groebner computation was unavailable.
    """
    u = sp.Dummy("semialg_saturation")
    try:
        basis = sp.groebner(
            (*generators, 1 - u * sp.expand(splitter)),
            u,
            *variables,
            order="lex",
            extension=True,
        )
    except (
        sp.PolynomialError,
        sp.polys.polyerrors.CoercionFailed,
        sp.polys.polyerrors.GeneratorsError,
        ValueError,
        TypeError,
        NotImplementedError,
    ):
        return None
    eliminated = tuple(
        sp.expand(poly.as_expr()) for poly in basis.polys if u not in poly.as_expr().free_symbols
    )
    return eliminated or tuple()


def _same_reduced_variety(left: EqualityIdealContext, right: EqualityIdealContext) -> bool:
    """Return whether two equality ideals have the same reduced zero set."""
    return _variety_subset(left, right) and _variety_subset(right, left)


def _lex_basis_generators(
    context: EqualityIdealContext,
) -> tuple[sp.Expr, ...] | None:
    """Return a deterministic lexicographic basis for triangular split tests."""
    try:
        return lex_basis_generators(context.generators, context.variables, extension=True)
    except (
        sp.PolynomialError,
        sp.polys.polyerrors.CoercionFailed,
        sp.polys.polyerrors.GeneratorsError,
        ValueError,
        TypeError,
        NotImplementedError,
    ):
        return None


def _triangular_splitter(
    context: EqualityIdealContext,
    *,
    kind: str,
    basis: tuple[sp.Expr, ...] | None = None,
) -> sp.Expr | None:
    """Return an initial or separant exposed by a lexicographic triangular slice."""
    if basis is None:
        basis = _lex_basis_generators(context)
    if not basis or len(context.variables) < 2:
        return None
    for count in range(2, len(context.variables) + 1):
        allowed = set(context.variables[:count])
        restricted = tuple(poly for poly in basis if poly.free_symbols.issubset(allowed))
        variable = context.variables[count - 1]
        candidates = tuple(poly for poly in restricted if poly.has(variable))
        if not candidates:
            continue
        polynomial = min(candidates, key=lambda poly: int(sp.degree(poly, variable)))
        if kind == "initial":
            degree = int(sp.degree(polynomial, variable))
            splitter = sp.expand(sp.Poly(polynomial, variable).coeff_monomial(variable**degree))
        elif kind == "separant":
            splitter = sp.expand(sp.diff(polynomial, variable))
        else:
            raise ValueError(f"unknown triangular splitter kind: {kind}")
        if splitter == 0 or not splitter.free_symbols.intersection(context.variables):
            continue
        return splitter
    return None


def _split_by_certified_splitter(
    context: EqualityIdealContext, splitter: sp.Expr
) -> tuple[tuple[sp.Expr, ...], tuple[sp.Expr, ...]] | None:
    """Apply the exact closed/saturated variety identity for one splitter."""
    if context.in_radical(splitter):
        return None
    closed = EqualityIdealContext.build((*context.generators, splitter), context.variables)
    if closed.inconsistent or _same_reduced_variety(closed, context):
        return None
    saturated_generators = _saturate_generators(context.generators, splitter, context.variables)
    if saturated_generators is None:
        return None
    saturated = EqualityIdealContext.build(saturated_generators, context.variables)
    if saturated.inconsistent or _same_reduced_variety(saturated, context):
        return None
    if _same_reduced_variety(closed, saturated):
        return None
    return _normalized_context_generators(closed), _normalized_context_generators(saturated)


def _initial_split_generators(
    context: EqualityIdealContext,
    *,
    basis: tuple[sp.Expr, ...] | None = None,
) -> tuple[tuple[sp.Expr, ...], tuple[sp.Expr, ...]] | None:
    """Split at a nonconstant triangular initial using exact saturation."""
    splitter = _triangular_splitter(context, kind="initial", basis=basis)
    return None if splitter is None else _split_by_certified_splitter(context, splitter)


def _separant_split_generators(
    context: EqualityIdealContext,
    *,
    basis: tuple[sp.Expr, ...] | None = None,
) -> tuple[tuple[sp.Expr, ...], tuple[sp.Expr, ...]] | None:
    """Split at a triangular separant using exact saturation."""
    splitter = _triangular_splitter(context, kind="separant", basis=basis)
    return None if splitter is None else _split_by_certified_splitter(context, splitter)


def _splitter_search_cost(splitter: sp.Expr, variables: Sequence[sp.Symbol]):
    """Cheap structural score for exact branch proposals.

    Lower-degree, sparser conditions involving fewer ambient variables are
    tried first.  The score never certifies a split; exact saturation still
    decides whether a proposal is admissible.
    """
    expression = sp.expand(splitter)
    active = tuple(v for v in variables if expression.has(v))
    try:
        poly = sp.Poly(expression, *variables, domain=sp.QQ)
        degree = int(poly.total_degree())
        terms = len(poly.terms())
    except (sp.PolynomialError, sp.polys.polyerrors.CoercionFailed, ValueError):
        degree = 10**9
        terms = 10**9
    return (degree, len(active), terms, sp.default_sort_key(expression))


def _saturation_split_generators(
    context: EqualityIdealContext,
) -> tuple[tuple[sp.Expr, ...], tuple[sp.Expr, ...]] | None:
    """Find a certified proper split from systematic triangular conditions.

    The fallback searches more than ambient coordinates: it includes every
    nonconstant initial, separant, and lower-variable coefficient exposed by a
    lexicographic triangular basis, plus exact squarefree factors of those
    conditions.  Each candidate remains only a proposal; the branch is accepted
    exclusively through the exact closed/saturation identity and strictness
    checks in :func:`_split_by_certified_splitter`.
    """
    if context.inconsistent or not context.generators:
        return None

    candidates: list[sp.Expr] = list(context.variables)
    basis = _lex_basis_generators(context) or _normalized_context_generators(context)
    for polynomial in basis:
        leader = _leader_of(polynomial, context.variables)
        if leader is None:
            continue
        try:
            poly = sp.Poly(polynomial, leader)
        except (sp.PolynomialError, ValueError, TypeError):
            continue
        degree = int(poly.degree())
        if degree > 0:
            initial = sp.expand(poly.LC())
            separant = sp.expand(sp.diff(polynomial, leader))
            if initial.free_symbols.intersection(context.variables):
                candidates.append(initial)
            if separant.free_symbols.intersection(context.variables):
                candidates.append(separant)
        for coefficient in poly.all_coeffs():
            coefficient = sp.expand(coefficient)
            if (
                coefficient != 0
                and coefficient.free_symbols.intersection(context.variables)
                and not coefficient.has(leader)
            ):
                candidates.append(coefficient)

        active = tuple(v for v in context.variables if polynomial.has(v))
        if len(active) == 1:
            factors = _squarefree_factor_atoms(polynomial, context.variables)
            if factors is not None:
                candidates.extend(factors)

    # Factoring coefficient conditions can expose a strict branch even when the
    # product itself does not.  This is still exact reduced-set branching.
    expanded: list[sp.Expr] = []
    for candidate in candidates:
        candidate = _canonical_factor(candidate, context.variables)
        expanded.append(candidate)
        factors = _squarefree_factor_atoms(candidate, context.variables)
        if factors is not None and len(factors) > 1:
            expanded.extend(factors)

    seen: set[sp.Expr] = set()
    for splitter in sorted(expanded, key=lambda h: _splitter_search_cost(h, context.variables)):
        splitter = _canonical_factor(splitter, context.variables)
        if splitter in seen or splitter == 0:
            continue
        seen.add(splitter)
        split = _split_by_certified_splitter(context, splitter)
        if split is not None:
            return split
    return None


def _strict_coefficient_sign(value: sp.Expr) -> int | None:
    """Return the certified strict sign of an exact scalar, if known."""
    value = sp.simplify(value)
    if value.is_positive is True:
        return 1
    if value.is_negative is True:
        return -1
    if value.is_zero is True:
        return 0
    return None


def _even_monomial_sum_zero_generators(
    expression: sp.Expr, variables: Sequence[sp.Symbol]
) -> tuple[sp.Expr, ...] | None:
    """Return an exact real-radical replacement for a definite even-monomial sum.

    If ``p = sum c_a x^(2a)`` and all nonzero ``c_a`` have the same strict
    sign, then over the reals ``p = 0`` iff every monomial ``x^a`` is zero.
    The returned equations therefore define exactly the same real zero set and
    may expose a lower-dimensional real radical or certify real emptiness.
    """
    try:
        poly = sp.Poly(sp.expand(expression), *variables, extension=True)
    except (sp.PolynomialError, ValueError, TypeError, NotImplementedError):
        return None
    if poly.is_zero:
        return None

    sign: int | None = None
    roots: list[sp.Expr] = []
    for exponent, coefficient in poly.terms():
        coefficient_sign = _strict_coefficient_sign(sp.sympify(coefficient))
        if coefficient_sign in (None, 0):
            return None
        if sign is None:
            sign = coefficient_sign
        elif coefficient_sign != sign:
            return None
        if any(int(power) % 2 for power in exponent):
            return None
        root = sp.Integer(1)
        for variable, power in zip(variables, exponent, strict=True):
            root *= variable ** (int(power) // 2)
        root = sp.expand(root)
        if root not in roots:
            roots.append(root)

    # A single nonconstant monomial is already handled by ordinary radical /
    # monomial decomposition and would otherwise cause pointless recursion.
    if len(roots) == 1 and roots[0] != 1:
        return None
    return tuple(roots)


def _quadratic_zero_generators(
    expression: sp.Expr, variables: Sequence[sp.Symbol]
) -> tuple[sp.Expr, ...] | None:
    """Return affine equations for the zero set of a certified PSD quadratic.

    For ``q(x)=x^T A x+b^T x+c`` with ``A`` positive semidefinite, ``q`` has
    minimum zero exactly when its real zero set is the affine stationary set
    ``grad(q)=0``.  All principal minors are checked exactly; no floating-point
    eigenvalue test is used.
    """
    try:
        poly = sp.Poly(sp.expand(expression), *variables, extension=True)
    except (sp.PolynomialError, ValueError, TypeError, NotImplementedError):
        return None
    if poly.total_degree() != 2:
        return None

    hessian = sp.Matrix([[sp.diff(poly.as_expr(), u, v) / 2 for v in variables] for u in variables])
    n = len(variables)
    # A symmetric matrix is PSD iff every principal minor is nonnegative.
    for size in range(1, n + 1):
        for indices in combinations(range(n), size):
            minor = sp.simplify(hessian.extract(indices, indices).det())
            if minor.is_nonnegative is not True:
                return None

    gradients = tuple(sp.expand(sp.diff(poly.as_expr(), variable)) for variable in variables)
    nonzero_gradients = tuple(gradient for gradient in gradients if gradient != 0)
    if not nonzero_gradients:
        return None
    try:
        stationary = sp.linsolve(nonzero_gradients, variables)
    except (ValueError, TypeError, NotImplementedError):
        return None
    if stationary is sp.EmptySet or stationary == sp.EmptySet:
        return None
    try:
        point = next(iter(stationary))
    except (StopIteration, TypeError):
        return None
    residual = sp.simplify(poly.as_expr().subs(dict(zip(variables, point, strict=True))))
    if residual != 0:
        return None
    return nonzero_gradients


def _real_radical_reduction_generators(
    context: EqualityIdealContext,
) -> tuple[sp.Expr, ...] | None:
    """Return a strictly stronger ideal defining the same *real* zero set.

    Supported certificates are definite sums of even monomials and PSD
    quadratics with zero minimum.  Replacements are accepted only when they
    change the reduced complex variety; this prevents recursive no-op rewrites.
    """
    basis = _normalized_context_generators(context)
    for index, polynomial in enumerate(basis):
        replacement = _even_monomial_sum_zero_generators(polynomial, context.variables)
        if replacement is None:
            replacement = _quadratic_zero_generators(polynomial, context.variables)
        if replacement is None:
            continue
        generators = (*basis[:index], *replacement, *basis[index + 1 :])
        reduced = EqualityIdealContext.build(generators, context.variables)
        if reduced.inconsistent:
            return (sp.Integer(1),)
        if _same_reduced_variety(reduced, context):
            continue
        return _normalized_context_generators(reduced)
    return None


def _is_linear_piece(context: EqualityIdealContext) -> bool:
    if context.groebner_basis is None:
        return True
    return all(poly.total_degree() <= 1 for poly in context.groebner_basis.polys)


def _ideal_generators_equal(
    left: Sequence[sp.Expr],
    right: Sequence[sp.Expr],
    variables: Sequence[sp.Symbol],
) -> bool:
    """Return exact ideal equality through canonical Groebner generators."""
    try:
        left_context = EqualityIdealContext.build(left, variables)
        right_context = EqualityIdealContext.build(right, variables)
    except RationalUnivariateError:
        return False
    return _normalized_context_generators(left_context) == _normalized_context_generators(
        right_context
    )


def _regular_chain_from_lex_basis(
    context: EqualityIdealContext,
) -> tuple[sp.Expr, ...] | None:
    """Extract a triangular chain with one polynomial per leader.

    The ambient variable order is used as a ranking.  The returned chain is
    ordered by increasing leader.  This is only a candidate; regularity,
    squarefreeness and equality with the branch are certified separately.
    """
    basis = _lex_basis_generators(context)
    if not basis:
        return tuple() if not context.generators else None
    index = {variable: i for i, variable in enumerate(context.variables)}
    by_leader: dict[sp.Symbol, list[sp.Expr]] = {}
    for polynomial in basis:
        active = [variable for variable in context.variables if polynomial.has(variable)]
        if not active:
            continue
        leader = max(active, key=index.__getitem__)
        by_leader.setdefault(leader, []).append(polynomial)
    if not by_leader:
        return None
    chain = []
    for leader in sorted(by_leader, key=index.__getitem__):
        candidates = by_leader[leader]
        chosen = min(
            candidates,
            key=lambda poly: (
                int(sp.degree(poly, leader)),
                sp.default_sort_key(poly),
            ),
        )
        chain.append(sp.expand(chosen))
    return tuple(chain)


def _chain_leader_data(
    chain: Sequence[sp.Expr],
    variables: Sequence[sp.Symbol],
) -> tuple[tuple[sp.Symbol, sp.Expr, sp.Expr], ...] | None:
    """Return ``(leader, initial, separant)`` for a triangular chain."""
    index = {variable: i for i, variable in enumerate(variables)}
    data: list[tuple[sp.Symbol, sp.Expr, sp.Expr]] = []
    leaders: set[sp.Symbol] = set()
    previous_index = -1
    for polynomial in chain:
        active = [variable for variable in variables if polynomial.has(variable)]
        if not active:
            return None
        leader = max(active, key=index.__getitem__)
        leader_index = index[leader]
        if leader in leaders or leader_index <= previous_index:
            return None
        degree = int(sp.degree(polynomial, leader))
        if degree <= 0:
            return None
        initial = sp.expand(sp.Poly(polynomial, leader).coeff_monomial(leader**degree))
        separant = sp.expand(sp.diff(polynomial, leader))
        if initial == 0 or separant == 0:
            return None
        data.append((leader, initial, separant))
        leaders.add(leader)
        previous_index = leader_index
    return tuple(data)


def _saturate_by_sequence(
    generators: Sequence[sp.Expr],
    splitters: Sequence[sp.Expr],
    variables: Sequence[sp.Symbol],
) -> tuple[sp.Expr, ...] | None:
    """Apply exact iterated saturation by a finite splitter sequence."""
    current = tuple(generators)
    for splitter in splitters:
        saturated = _saturate_generators(current, splitter, variables)
        if saturated is None:
            return None
        current = saturated
    try:
        return _normalized_context_generators(EqualityIdealContext.build(current, variables))
    except RationalUnivariateError:
        return None


def _leader_of(polynomial: sp.Expr, variables: Sequence[sp.Symbol]) -> sp.Symbol | None:
    index = {variable: i for i, variable in enumerate(variables)}
    active = [variable for variable in variables if polynomial.has(variable)]
    return max(active, key=index.__getitem__) if active else None


def _candidate_competing_leader_polynomials(
    context: EqualityIdealContext,
    *,
    basis: tuple[sp.Expr, ...] | None = None,
) -> tuple[sp.Expr, ...]:
    """Return exact equations worth comparing for regular-gcd branching.

    A reduced lex basis can hide the polynomial pair whose gcd degree changes
    on a lower-dimensional coefficient stratum.  Keep both the lex equations
    and the current branch generators, canonicalizing duplicates only after
    collecting them.  This is search data only; every resulting splitter is
    still accepted through the exact closed/saturation identity.
    """
    basis = _lex_basis_generators(context) if basis is None else basis
    candidates = list(basis or tuple()) + list(context.generators)
    unique: list[sp.Expr] = []
    for polynomial in candidates:
        polynomial = sp.expand(polynomial)
        if polynomial == 0 or polynomial in unique:
            continue
        unique.append(polynomial)
    return tuple(unique)


def _regular_gcd_split_generators(
    context: EqualityIdealContext,
    *,
    basis: tuple[sp.Expr, ...] | None = None,
) -> tuple[tuple[sp.Expr, ...], tuple[sp.Expr, ...]] | None:
    """Split on exact regular-gcd/subresultant obstructions.

    For two equations with the same leader, the degree of their gcd in the
    quotient field of the lower variables can change only where a principal
    subresultant coefficient vanishes.  Defective subresultants and
    pseudo-remainders expose additional coefficient conditions needed when
    leading coefficients themselves specialize.

    The coefficients are *not* trusted as a decomposition theorem by
    themselves.  They merely propose ``h``; :func:`_split_by_certified_splitter`
    proves ``V(I) = V(I + <h>) union V(I : h**inf)`` and that both children are
    strict before a branch is accepted.
    """
    polynomials = _candidate_competing_leader_polynomials(context, basis=basis)
    if not polynomials:
        return None
    groups: dict[sp.Symbol, list[sp.Expr]] = {}
    for polynomial in polynomials:
        leader = _leader_of(polynomial, context.variables)
        if leader is not None:
            groups.setdefault(leader, []).append(polynomial)

    candidates: list[sp.Expr] = []
    for leader, group in groups.items():
        if len(group) < 2:
            continue
        ordered = sorted(
            group,
            key=lambda expr: (int(sp.degree(expr, leader)), sp.default_sort_key(expr)),
        )
        for left, right in combinations(ordered, 2):
            if int(sp.degree(left, leader)) <= 0 or int(sp.degree(right, leader)) <= 0:
                continue
            try:
                prs = subresultant_prs(left, right, leader)
            except (
                ArithmeticError,
                TypeError,
                ValueError,
                NotImplementedError,
                sp.PolynomialError,
                sp.polys.polyerrors.CoercionFailed,
            ):
                prs = None
            if prs is not None:
                # PSCs stratify the possible regular-gcd degrees.  Coefficients
                # of defective PRS entries handle specializations where the
                # nominal leading coefficient disappears.
                for coefficient in prs.principal_coefficients:
                    coefficient = sp.expand(coefficient)
                    if (
                        coefficient != 0
                        and coefficient.free_symbols.intersection(context.variables)
                        and not coefficient.has(leader)
                    ):
                        candidates.append(coefficient)
                for entry in prs.polynomials[2:]:
                    try:
                        entry_poly = sp.Poly(entry.as_expr(), leader)
                    except (sp.PolynomialError, TypeError, ValueError):
                        continue
                    for coefficient in entry_poly.all_coeffs():
                        coefficient = sp.expand(coefficient)
                        if (
                            coefficient != 0
                            and coefficient.free_symbols.intersection(context.variables)
                            and not coefficient.has(leader)
                        ):
                            candidates.append(coefficient)
            try:
                left_poly = sp.Poly(left, leader)
                right_poly = sp.Poly(right, leader)
                for dividend, divisor in ((left_poly, right_poly), (right_poly, left_poly)):
                    if divisor.is_zero or divisor.degree() <= 0:
                        continue
                    remainder = dividend.prem(divisor)
                    for coefficient in remainder.all_coeffs():
                        coefficient = sp.expand(coefficient)
                        if (
                            coefficient != 0
                            and coefficient.free_symbols.intersection(context.variables)
                            and not coefficient.has(leader)
                        ):
                            candidates.append(coefficient)
            except (sp.PolynomialError, TypeError, ValueError):
                pass

    seen: set[sp.Expr] = set()
    for splitter in sorted(candidates, key=lambda h: _splitter_search_cost(h, context.variables)):
        splitter = _canonical_factor(splitter, context.variables)
        if splitter in seen:
            continue
        seen.add(splitter)
        split = _split_by_certified_splitter(context, splitter)
        if split is not None:
            return split
    return None


def _regular_chain_obstruction_split_generators(
    context: EqualityIdealContext,
    *,
    basis: tuple[sp.Expr, ...] | None = None,
) -> tuple[tuple[tuple[sp.Expr, ...], tuple[sp.Expr, ...]], str] | None:
    """Return a split at the first failed regular-chain proof obligation.

    This turns the squarefree-regular-chain *verifier* into a source of exact
    branch conditions.  Initials are tested against the prefix already
    saturated by earlier initials; separants are tested after adjoining the
    current chain polynomial.  If either is a zero divisor, the same polynomial
    is used as a global closed/saturation splitter.
    """
    chain = _regular_chain_from_lex_basis(context)
    if not chain:
        return None
    data = _chain_leader_data(chain, context.variables)
    if data is None:
        return None

    saturated_prefix: tuple[sp.Expr, ...] = tuple()
    initials: list[sp.Expr] = []
    for polynomial, (_leader, initial, separant) in zip(chain, data, strict=True):
        initial_sat = _saturate_generators(saturated_prefix, initial, context.variables)
        if initial_sat is None:
            return None
        if not _ideal_generators_equal(saturated_prefix, initial_sat, context.variables):
            split = _split_by_certified_splitter(context, initial)
            if split is not None:
                return split, "regular_zero_divisor_initial_split"

        prefix = (*saturated_prefix, polynomial)
        separant_sat = _saturate_generators(prefix, separant, context.variables)
        if separant_sat is None:
            return None
        if not _ideal_generators_equal(prefix, separant_sat, context.variables):
            split = _split_by_certified_splitter(context, separant)
            if split is not None:
                return split, "regular_zero_divisor_separant_split"

        initials.append(initial)
        updated = _saturate_by_sequence(prefix, initials, context.variables)
        if updated is None:
            return None
        saturated_prefix = updated

    # A triangular candidate may satisfy all local regularity checks but fail
    # to generate the whole current branch after saturation.  In that case,
    # let regular-gcd data from the omitted/competing equations propose the
    # next certified splitter.
    if not _ideal_generators_equal(saturated_prefix, context.generators, context.variables):
        split = _regular_gcd_split_generators(context, basis=basis)
        if split is not None:
            return split, "regular_competing_leader_gcd_split"
    return None


def _squarefree_regular_chain_certified(context: EqualityIdealContext) -> bool:
    """Certify that the branch is the saturation of a squarefree regular chain.

    See Kalkbrener (1993) and the regular-chain references [Kalkbrener1993],
    [Lazard1991], and [ChenEtAl2013] in ``docs/references.md``.

    Over a characteristic-zero field, the saturated ideal of a squarefree
    regular chain is radical and unmixed, hence its variety is equidimensional
    of codimension equal to the chain length.  We do not trust triangular shape
    alone: initials and separants are certified regular by exact saturation, and
    the resulting saturated chain is required to define exactly the same ideal
    as the current branch.
    """
    if context.inconsistent:
        return True
    chain = _regular_chain_from_lex_basis(context)
    if not chain:
        return False
    codimension = len(context.variables) - context.dimension
    if len(chain) != codimension:
        return False
    data = _chain_leader_data(chain, context.variables)
    if data is None:
        return False

    prefix: tuple[sp.Expr, ...] = tuple()
    initials: list[sp.Expr] = []
    saturated_prefix: tuple[sp.Expr, ...] = tuple()
    for polynomial, (_leader, initial, separant) in zip(chain, data, strict=True):
        # Regularity is relative to the prefix already saturated by all earlier
        # initials, not merely to the raw triangular equations.
        initial_sat = _saturate_generators(saturated_prefix, initial, context.variables)
        if initial_sat is None or not _ideal_generators_equal(
            saturated_prefix, initial_sat, context.variables
        ):
            return False
        prefix = (*prefix, polynomial)
        initials.append(initial)
        saturated = _saturate_by_sequence(prefix, initials, context.variables)
        if saturated is None:
            return False
        saturated_prefix = saturated
        # Squarefree regularity: the separant is a non-zero-divisor modulo the
        # saturated prefix. Exact ideal equality after saturation certifies it.
        separant_sat = _saturate_generators(saturated_prefix, separant, context.variables)
        if separant_sat is None or not _ideal_generators_equal(
            saturated_prefix, separant_sat, context.variables
        ):
            return False

    saturated_chain = saturated_prefix
    if saturated_chain is None:
        return False
    if not _ideal_generators_equal(saturated_chain, context.generators, context.variables):
        return False
    return True


def _equidimensionality_certified(context: EqualityIdealContext) -> bool:
    """Conservatively certify equidimensionality of the reduced zero set."""
    if context.inconsistent or context.dimension <= 0:
        return True
    if not context.generators:
        return True
    if _is_linear_piece(context):
        return True
    # A principal proper ideal defines a hypersurface.  Factor splitting has
    # already removed repeated and reducible factors before terminal pieces are
    # accepted, so its reduced components all have codimension one.
    if len(context.generators) == 1:
        return context.dimension == len(context.variables) - 1
    codimension = len(context.variables) - context.dimension
    # A proper ideal generated by height(I) elements is a complete intersection
    # in a polynomial ring over a field and is unmixed.  Use the reduced
    # Groebner basis here rather than the caller's original generator list:
    # decomposition branches retain inherited product equations that may be
    # redundant after a factor has selected one component.
    basis_size = (
        len(context.groebner_basis.polys)
        if context.groebner_basis is not None
        else len(context.generators)
    )
    if codimension > 0 and basis_size == codimension:
        return True
    # General triangular fallback: a squarefree regular chain has a radical,
    # unmixed saturation in characteristic zero.  This catches prime and
    # radical complete-intersection-like branches whose reduced Groebner basis
    # contains more generators than their height (for example monomial-curve
    # ideals), without guessing primality.
    if _squarefree_regular_chain_certified(context):
        return True
    # Any exactly certified prime ideal is equidimensional. This hook lets the
    # minimal-prime layer strengthen downstream ARS eligibility without making
    # primality a prerequisite for ordinary decomposition.
    if _prime_certificate_method(context) is not None:
        return True
    return False


def _piece_formula(context: EqualityIdealContext) -> sp.Expr:
    """Return the real zero-set formula represented by an equality ideal."""
    if context.inconsistent:
        return sp.false
    if not context.generators:
        return sp.true
    return sp.And(*(sp.Eq(generator, 0) for generator in context.generators))


def _real_dimension(context: EqualityIdealContext) -> int:
    """Return exact real dimension, avoiding CAD for affine-linear pieces."""
    if context.inconsistent:
        return -1
    if not context.generators:
        return len(context.variables)
    if _is_linear_piece(context):
        return context.dimension
    from .regions.operations import region_dimension

    return region_dimension(_piece_formula(context), context.variables)


def _real_point_contexts(context: EqualityIdealContext) -> tuple[EqualityIdealContext, ...]:
    """Replace a finite real zero set by exact maximal-ideal point pieces."""
    from .cad_algorithms.selected_samples import extract_selected_cad_samples

    samples = extract_selected_cad_samples(_piece_formula(context), context.variables)
    out: list[EqualityIdealContext] = []
    seen: set[tuple[sp.Expr, ...]] = set()
    for sample in samples:
        point = dict(sample.point)
        if not all(variable in point for variable in context.variables):
            continue
        equations = tuple(sp.expand(variable - point[variable]) for variable in context.variables)
        key = tuple(sp.simplify(point[variable]) for variable in context.variables)
        if key in seen:
            continue
        seen.add(key)
        out.append(EqualityIdealContext.build(equations, context.variables))
    return tuple(out)


def _variety_subset(smaller: EqualityIdealContext, larger: EqualityIdealContext) -> bool:
    """Return whether V(smaller) is contained in V(larger) exactly."""
    if smaller.inconsistent:
        return True
    if larger.inconsistent:
        return smaller.inconsistent
    return all(smaller.in_radical(generator) for generator in larger.generators)


def _is_squarefree_monomial_ideal(context: EqualityIdealContext) -> bool:
    basis = context.groebner_basis
    if basis is None or context.inconsistent:
        return False
    for poly in basis.polys:
        terms = poly.terms()
        if len(terms) != 1:
            return False
        exponent, _coefficient = terms[0]
        if any(power > 1 for power in exponent):
            return False
    return True


def _principal_squarefree(context: EqualityIdealContext) -> bool:
    if len(context.generators) != 1:
        return False
    try:
        _unit, factors = sp.factor_list(context.generators[0], *context.variables, extension=True)
    except (sp.PolynomialError, ValueError, TypeError, NotImplementedError):
        return False
    return bool(factors) and all(multiplicity == 1 for _factor, multiplicity in factors)


def _radicality_certified(context: EqualityIdealContext) -> bool:
    """Conservatively certify that the current ideal is radical."""
    if context.inconsistent or not context.generators or _is_linear_piece(context):
        return True
    if _principal_squarefree(context) or _is_squarefree_monomial_ideal(context):
        return True
    return _squarefree_regular_chain_certified(context)


def _irreducible_polynomial(expression: sp.Expr, variables: Sequence[sp.Symbol]) -> bool:
    try:
        _unit, factors = sp.factor_list(expression, *variables, extension=True)
    except (sp.PolynomialError, ValueError, TypeError, NotImplementedError):
        return False
    nonconstant = [
        factor
        for factor, _multiplicity in factors
        if sp.sympify(factor).free_symbols.intersection(variables)
    ]
    return len(nonconstant) == 1 and factors[0][1] == 1


def _constant_nonzero(expression: sp.Expr, variables: Sequence[sp.Symbol]) -> bool:
    expression = sp.cancel(expression)
    return not expression.free_symbols.intersection(variables) and expression != 0


def _clear_rational_function_denominator(expression: sp.Expr) -> sp.Expr:
    """Return a primitive polynomial numerator for a rational expression."""
    numerator, _denominator = sp.fraction(sp.cancel(expression))
    return sp.expand(numerator)


def certify_triangular_primality(
    equations: Iterable[sp.Expr | sp.Equality],
    variables: Sequence[sp.Symbol],
) -> TriangularPrimalityCertificate:
    """Certify a triangular ideal prime or expose an exact reducible stage.

    The proof is successive: after a squarefree-regular saturation certificate,
    each leader polynomial is factored over the exact fraction field of the
    preceding quotient.  Every factorization is independently reconstructed
    before irreducibility or a split is accepted.
    """
    vars_ = tuple(variables)
    context = EqualityIdealContext.build(_normalize_equations(equations, vars_), vars_)
    chain = _regular_chain_from_lex_basis(context) or tuple()
    data = _chain_leader_data(chain, vars_)
    if not chain or data is None or not _squarefree_regular_chain_certified(context):
        return TriangularPrimalityCertificate(
            vars_, chain, tuple(), False, False, method="not_squarefree_regular"
        )
    leaders = tuple(item[0] for item in data)
    free = tuple(variable for variable in vars_ if variable not in leaders)
    field = RationalFunctionField(free)
    stages: list[TriangularPrimalityStage] = []
    nonlinear = False
    try:
        for polynomial, (leader, _initial, _separant) in zip(chain, data, strict=True):
            poly = sp.Poly(polynomial, leader)
            degree = int(poly.degree())
            if degree <= 0:
                return TriangularPrimalityCertificate(
                    vars_, chain, tuple(stages), False, False, method="invalid_leader_degree"
                )
            coeffs = [
                tower_evaluate(poly.coeff_monomial(leader**i), field) for i in range(degree + 1)
            ]
            leading = coeffs[-1]
            if field.is_zero(leading):
                return TriangularPrimalityCertificate(
                    vars_, chain, tuple(stages), False, False, method="vanishing_initial"
                )
            leading_inv = field.inv(leading)
            monic = [coefficient * leading_inv for coefficient in coeffs]
            factorization = certified_factor_univariate(field, monic)
            factors = tuple(
                (factor, multiplicity)
                for factor, multiplicity in factorization.factors
                if len(factor) > 1
            )
            factor_degrees = tuple(len(factor) - 1 for factor, _m in factors)
            stages.append(
                TriangularPrimalityStage(
                    leader,
                    degree,
                    factor_degrees,
                    factorization.reconstruction_verified,
                    factorization.method,
                )
            )
            irreducible = (
                len(factors) == 1 and factors[0][1] == 1 and len(factors[0][0]) - 1 == degree
            )
            if not irreducible:
                splitters: list[sp.Expr] = []
                for factor, _multiplicity in factors:
                    expr = sp.expand(
                        sum(
                            field_element_to_expr(coefficient, field) * leader**i
                            for i, coefficient in enumerate(factor)
                        )
                    )
                    numerator, _denominator = sp.fraction(sp.cancel(expr))
                    if numerator != 0:
                        splitters.append(sp.expand(numerator))
                return TriangularPrimalityCertificate(
                    vars_,
                    chain,
                    tuple(stages),
                    len(splitters) >= 2,
                    False,
                    tuple(splitters),
                    "certified_reducible_stage",
                )
            factor = factors[0][0]
            nonlinear |= degree > 1
            field = MonogenicFunctionField(field, leader, tuple(factor))
    except (
        FunctionFieldError,
        ZeroDivisionError,
        TypeError,
        ValueError,
        sp.PolynomialError,
        NotImplementedError,
    ):
        return TriangularPrimalityCertificate(
            vars_, chain, tuple(stages), False, False, method="field_tower_incomplete"
        )
    return TriangularPrimalityCertificate(
        vars_,
        chain,
        tuple(stages),
        True,
        True,
        method=(
            "certified_algebraic_function_field_tower_prime"
            if nonlinear
            else "certified_fraction_field_graph_prime"
        ),
    )


def _tower_triangular_primality_analysis(
    context: EqualityIdealContext,
) -> _TriangularPrimalityAnalysis:
    """Analyze a squarefree regular chain in an exact algebraic function-field tower.

    Each leader polynomial is factored over the fraction field of the preceding
    quotient. Irreducible stages are adjoined as monogenic extensions; a
    reducible stage yields exact polynomial splitters in the ambient variables.
    """
    chain = _regular_chain_from_lex_basis(context)
    data = _chain_leader_data(chain or tuple(), context.variables)
    if not chain or data is None or not _squarefree_regular_chain_certified(context):
        return _TriangularPrimalityAnalysis()
    leaders = tuple(item[0] for item in data)
    free = tuple(variable for variable in context.variables if variable not in leaders)
    field = RationalFunctionField(free)
    nonlinear = False
    try:
        for polynomial, (leader, _initial, _separant) in zip(chain, data, strict=True):
            poly = sp.Poly(polynomial, leader)
            degree = poly.degree()
            if degree <= 0:
                return _TriangularPrimalityAnalysis()
            coeffs = [
                tower_evaluate(poly.coeff_monomial(leader**i), field) for i in range(degree + 1)
            ]
            leading = coeffs[-1]
            if field.is_zero(leading):
                return _TriangularPrimalityAnalysis()
            leading_inv = field.inv(leading)
            monic = [coefficient * leading_inv for coefficient in coeffs]
            factorization = certified_factor_univariate(field, monic)
            if not factorization.reconstruction_verified:
                return _TriangularPrimalityAnalysis()
            factors = factorization.factors
            nonconstant = [
                (factor, multiplicity) for factor, multiplicity in factors if len(factor) > 1
            ]
            if len(nonconstant) != 1 or nonconstant[0][1] != 1:
                splitters = []
                for factor, _multiplicity in nonconstant:
                    expr = sp.expand(
                        sum(
                            field_element_to_expr(coefficient, field) * leader**i
                            for i, coefficient in enumerate(factor)
                        )
                    )
                    numerator, _denominator = sp.fraction(sp.cancel(expr))
                    if numerator != 0:
                        splitters.append(sp.expand(numerator))
                if len(splitters) >= 2:
                    return _TriangularPrimalityAnalysis(splitters=tuple(splitters))
                return _TriangularPrimalityAnalysis()
            factor, multiplicity = nonconstant[0]
            if multiplicity != 1 or len(factor) - 1 != degree:
                return _TriangularPrimalityAnalysis()
            if degree > 1:
                nonlinear = True
            field = MonogenicFunctionField(field, leader, tuple(factor))
    except (
        FunctionFieldError,
        ZeroDivisionError,
        TypeError,
        ValueError,
        sp.PolynomialError,
        NotImplementedError,
    ):
        return _TriangularPrimalityAnalysis()
    return _TriangularPrimalityAnalysis(
        prime_method=(
            "algebraic_function_field_tower_prime"
            if nonlinear
            else "triangular_fraction_field_graph_prime"
        )
    )


def _triangular_primality_analysis(
    context: EqualityIdealContext,
) -> _TriangularPrimalityAnalysis:
    """Prove primality or expose a quotient-field factorization of a chain.

    The test proceeds in triangular leader order. Monic/unit-linear leaders are
    eliminated exactly, so the next nonlinear leader polynomial is tested over
    the fraction field of the remaining free variables, not in the ambient
    polynomial ring. This distinction is essential: ``(y-x**2, z**2-y)`` is
    reducible because the second equation becomes ``z**2-x**2`` after passing
    to the quotient by the first.

    SymPy currently cannot construct arbitrary towers of algebraic function
    fields. Therefore, after one genuinely nonlinear algebraic extension, later
    nonlinear stages are reported as undecided rather than guessed. Linear
    graph extensions after that stage remain certified because adjoining a
    linearly determined variable preserves a domain.
    """
    if context.inconsistent:
        return _TriangularPrimalityAnalysis()
    tower_analysis = _tower_triangular_primality_analysis(context)
    if tower_analysis.prime_method is not None or tower_analysis.splitters:
        return tower_analysis
    chain = _regular_chain_from_lex_basis(context)
    data = _chain_leader_data(chain or tuple(), context.variables)
    if not chain or data is None:
        return _TriangularPrimalityAnalysis()
    if not _squarefree_regular_chain_certified(context):
        return _TriangularPrimalityAnalysis()

    substitutions: dict[sp.Symbol, sp.Expr] = {}
    nonlinear_extension_seen = False

    for polynomial, (leader, _initial, _separant) in zip(chain, data, strict=True):
        reduced = sp.cancel(sp.expand(polynomial.subs(substitutions)))
        try:
            degree = int(sp.degree(reduced, leader))
        except (sp.PolynomialError, TypeError, ValueError):
            return _TriangularPrimalityAnalysis()
        if degree <= 0:
            return _TriangularPrimalityAnalysis()

        if degree == 1:
            poly = sp.Poly(reduced, leader)
            a = sp.cancel(poly.coeff_monomial(leader))
            b = sp.cancel(poly.coeff_monomial(1))
            # A linear polynomial over a domain's fraction field is irreducible.
            # We may solve it provided its leading coefficient is a nonzero
            # element of the already-certified domain. Before a nonlinear
            # extension, exact nonzeroness follows from rational-function-field
            # arithmetic. Afterwards require a scalar unit; otherwise deciding
            # whether the coefficient vanishes needs algebraic-extension
            # arithmetic that SymPy does not presently expose.
            if nonlinear_extension_seen:
                if not _constant_nonzero(a, context.variables):
                    return _TriangularPrimalityAnalysis()
            else:
                free_before = tuple(
                    variable
                    for variable in context.variables
                    if variable != leader
                    and variable not in substitutions
                    and variable in a.free_symbols
                )
                try:
                    domain = sp.QQ.frac_field(*free_before) if free_before else sp.QQ
                    if domain.convert(a) == domain.zero:
                        return _TriangularPrimalityAnalysis()
                except (TypeError, ValueError, sp.PolynomialError):
                    return _TriangularPrimalityAnalysis()
            substitutions[leader] = sp.cancel(-b / a)
            continue

        if nonlinear_extension_seen:
            return _TriangularPrimalityAnalysis()

        parameters = tuple(
            variable
            for variable in context.variables
            if variable != leader
            and variable not in substitutions
            and variable in reduced.free_symbols
        )
        try:
            domain = sp.QQ.frac_field(*parameters) if parameters else sp.QQ
            univariate = sp.Poly(reduced, leader, domain=domain)
            _unit, factors = univariate.factor_list()
        except (
            sp.PolynomialError,
            sp.polys.polyerrors.CoercionFailed,
            sp.polys.polyerrors.GeneratorsError,
            TypeError,
            ValueError,
            NotImplementedError,
        ):
            return _TriangularPrimalityAnalysis()

        nonconstant = [
            (factor, multiplicity) for factor, multiplicity in factors if factor.degree() > 0
        ]
        if len(nonconstant) != 1 or nonconstant[0][1] != 1:
            splitters: list[sp.Expr] = []
            for factor, _multiplicity in nonconstant:
                expression = _clear_rational_function_denominator(factor.as_expr())
                if expression != 0:
                    splitters.append(expression)
            if len(splitters) >= 2:
                return _TriangularPrimalityAnalysis(splitters=tuple(splitters))
            return _TriangularPrimalityAnalysis()

        factor, multiplicity = nonconstant[0]
        if multiplicity != 1 or factor.degree() != univariate.degree():
            return _TriangularPrimalityAnalysis()
        nonlinear_extension_seen = True

    method = (
        "successive_fraction_field_prime"
        if nonlinear_extension_seen
        else "triangular_fraction_field_graph_prime"
    )
    return _TriangularPrimalityAnalysis(prime_method=method)


def _prime_certificate_method(context: EqualityIdealContext) -> str | None:
    """Return the exact sufficient primality certificate used for ``context``.

    Triangular ideals are tested in the quotient/fraction field induced by the
    preceding chain, rather than by ambient factorization. See the regular-chain
    references in ``docs/references.md``.
    """
    if context.inconsistent:
        return None
    if not context.generators:
        return "zero_ideal"
    if _is_linear_piece(context):
        return "linear_prime"
    if len(context.generators) == 1 and _irreducible_polynomial(
        context.generators[0], context.variables
    ):
        return "irreducible_hypersurface"

    normalized = set(_normalized_context_generators(context))
    if normalized and all(expr in set(context.variables) for expr in normalized):
        return "coordinate_prime"

    return _triangular_primality_analysis(context).prime_method


def _ideal_intersection_generators(
    left: Sequence[sp.Expr], right: Sequence[sp.Expr], variables: Sequence[sp.Symbol]
) -> tuple[sp.Expr, ...] | None:
    """Compute ``<left> intersect <right>`` by exact elimination."""
    t = sp.Dummy("semialg_intersection")
    try:
        basis = sp.groebner(
            tuple(t * sp.expand(f) for f in left) + tuple((1 - t) * sp.expand(g) for g in right),
            t,
            *variables,
            order="lex",
            extension=True,
        )
    except (
        sp.PolynomialError,
        sp.polys.polyerrors.CoercionFailed,
        sp.polys.polyerrors.GeneratorsError,
        ValueError,
        TypeError,
        NotImplementedError,
    ):
        return None
    eliminated = tuple(
        sp.expand(poly.as_expr()) for poly in basis.polys if t not in poly.as_expr().free_symbols
    )
    try:
        return _normalized_context_generators(EqualityIdealContext.build(eliminated, variables))
    except RationalUnivariateError:
        return None


def _intersection_all(
    ideals: Sequence[Sequence[sp.Expr]], variables: Sequence[sp.Symbol]
) -> tuple[sp.Expr, ...] | None:
    if not ideals:
        return (sp.Integer(1),)
    current = tuple(ideals[0])
    for ideal in ideals[1:]:
        current = _ideal_intersection_generators(current, ideal, variables)
        if current is None:
            return None
    return current


def _recursive_regular_chain_prime_refinement(
    contexts: Sequence[EqualityIdealContext],
    variables: Sequence[sp.Symbol],
    *,
    max_pieces: int,
) -> tuple[EqualityIdealContext, ...]:
    """Recursively split certified reducible triangular stages.

    A splitter is used only when it comes from an exact factorization over the
    quotient/fraction field of the preceding regular-chain prefix. Each child is
    rebuilt through the ordinary certified equidimensional decomposition before
    further recursion. Final reconstruction in the public decomposition routine
    independently proves that the returned union still represents ``sqrt(I)``.
    """
    queue = list(contexts)
    terminal: list[EqualityIdealContext] = []
    seen: set[tuple[sp.Expr, ...]] = set()
    while queue:
        if len(queue) + len(terminal) > max_pieces:
            terminal.extend(queue)
            break
        context = queue.pop(0)
        key = _normalized_context_generators(context)
        if key in seen:
            continue
        seen.add(key)
        if _prime_certificate_method(context) is not None:
            terminal.append(context)
            continue
        analysis = _triangular_primality_analysis(context)
        if len(analysis.splitters) < 2:
            terminal.append(context)
            continue

        children: list[EqualityIdealContext] = []
        for splitter in analysis.splitters:
            try:
                branch_source = EqualityIdealContext.build(
                    (*context.generators, splitter), variables
                )
            except RationalUnivariateError:
                continue
            if branch_source.inconsistent or _ideal_generators_equal(
                branch_source.generators, context.generators, variables
            ):
                continue
            decomposition = _equidimensional_decomposition_impl(
                branch_source.generators,
                variables,
                max_pieces=max_pieces,
                _emit_certificate=True,
                _algebraic_semantics=True,
            )
            for piece in decomposition.pieces:
                try:
                    child = EqualityIdealContext.build(piece.equations, variables)
                except RationalUnivariateError:
                    continue
                if not child.inconsistent:
                    children.append(child)
        if len(children) < 2:
            terminal.append(context)
            continue
        queue.extend(children)

    # Remove duplicate/contained algebraic components exactly. If V(A) is
    # contained in V(B), A is redundant in a union whenever B is also present.
    reduced: list[EqualityIdealContext] = []
    for candidate in terminal:
        if any(
            _ideal_generators_equal(candidate.generators, other.generators, variables)
            for other in reduced
        ):
            continue
        reduced.append(candidate)
    keep: list[EqualityIdealContext] = []
    for i, candidate in enumerate(reduced):
        if any(i != j and _variety_subset(candidate, other) for j, other in enumerate(reduced)):
            continue
        keep.append(candidate)
    return tuple(keep)


def _canonical_ideal_tuple(
    generators: Sequence[sp.Expr], variables: Sequence[sp.Symbol]
) -> tuple[sp.Expr, ...]:
    return _normalized_context_generators(EqualityIdealContext.build(generators, variables))


def _recover_splitter_witness(
    parent: EqualityIdealContext,
    children: Sequence[Sequence[sp.Expr]],
) -> sp.Expr | None:
    """Recover an exact ``h`` when children are ``I+<h>`` and ``I:h^inf``."""
    if len(children) != 2:
        return None
    target = {_canonical_ideal_tuple(child, parent.variables) for child in children}
    candidates: list[sp.Expr] = []
    for child in children:
        for generator in child:
            candidate = _canonical_factor(generator, parent.variables)
            if candidate == 0 or parent.normal_form(candidate) == 0:
                continue
            if candidate not in candidates:
                candidates.append(candidate)
    for candidate in candidates:
        split = _split_by_certified_splitter(parent, candidate)
        if split is None:
            continue
        if {_canonical_ideal_tuple(piece, parent.variables) for piece in split} == target:
            return candidate
    return None


def recursive_regular_chain_decomposition(
    equations: Iterable[sp.Expr | sp.Equality],
    variables: Sequence[sp.Symbol],
    *,
    max_pieces: int | None = None,
) -> RegularChainDecomposition:
    """Recursively decompose a polynomial locus into certified regular chains.

    Branches are created only by exact closed/saturation identities or by
    exact quotient-field factorization.  Initials, separants and PRS principal
    coefficients drive the regular-gcd branching.  No terminal branch is called
    complete unless its ideal equals the saturation of a squarefree regular
    chain.
    """
    vars_ = tuple(variables)
    if not vars_:
        raise ValueError("regular-chain decomposition requires ambient variables")
    if max_pieces is not None and max_pieces < 1:
        raise ValueError("max_pieces must be positive when supplied")
    source = EqualityIdealContext.build(_normalize_equations(equations, vars_), vars_)
    source_degree = ideal_degree(source.generators, vars_)
    pending = [source]
    terminal: list[EqualityIdealContext] = []
    methods: set[str] = set()
    split_certificates: list[RegularChainSplitCertificate] = []
    seen: set[tuple[sp.Expr, ...]] = set()

    while pending:
        if max_pieces is not None and len(pending) + len(terminal) > max_pieces:
            terminal.extend(pending)
            pending.clear()
            break
        context = pending.pop()
        if context.inconsistent:
            continue
        key = _normalized_context_generators(context)
        if key in seen:
            continue
        seen.add(key)
        if _squarefree_regular_chain_certified(context):
            terminal.append(context)
            continue

        basis = _lex_basis_generators(context)
        obstruction = _regular_chain_obstruction_split_generators(context, basis=basis)
        if obstruction is not None:
            split, method = obstruction
        else:
            split = _initial_split_generators(context, basis=basis)
            method = "regular_initial_split"
        if split is None:
            split = _separant_split_generators(context, basis=basis)
            method = "regular_separant_split"
        if split is None:
            split = _regular_gcd_split_generators(context, basis=basis)
            method = "regular_gcd_subresultant_split"
        if split is None:
            analysis = _triangular_primality_analysis(context)
            if len(analysis.splitters) >= 2:
                children = []
                for splitter in analysis.splitters:
                    child = EqualityIdealContext.build((*context.generators, splitter), vars_)
                    if not child.inconsistent:
                        children.append(_normalized_context_generators(child))
                if len(children) >= 2:
                    split = tuple(children)  # type: ignore[assignment]
                    method = "regular_function_field_factor_split"
        if split is None:
            split = _saturation_split_generators(context)
            method = "regular_saturation_split"
        if split is None:
            for polynomial in _normalized_context_generators(context):
                factors = _distinct_factors(polynomial, vars_)
                if factors is None:
                    continue
                split = tuple(
                    _normalized_context_generators(
                        EqualityIdealContext.build((*context.generators, factor), vars_)
                    )
                    for factor in factors
                )
                method = "regular_squarefree_factor_split"
                break
        if split is None:
            terminal.append(context)
            continue
        methods.add(method)
        if max_pieces is not None and len(pending) + len(terminal) + len(split) > max_pieces:
            terminal.append(context)
            continue
        canonical_children = tuple(
            _canonical_ideal_tuple(generators, vars_) for generators in split
        )
        split_certificates.append(
            RegularChainSplitCertificate(
                parent=_normalized_context_generators(context),
                children=canonical_children,
                method=method,
                splitter=_recover_splitter_witness(context, canonical_children),
            )
        )
        for generators in canonical_children:
            child = EqualityIdealContext.build(generators, vars_)
            if not child.inconsistent:
                pending.append(child)

    # Exact union cleanup.
    reduced: list[EqualityIdealContext] = []
    for candidate in terminal:
        if any(_variety_subset(candidate, other) for other in reduced):
            continue
        reduced = [other for other in reduced if not _variety_subset(other, candidate)]
        reduced.append(candidate)

    components = tuple(
        RegularChainComponent(
            equations=_normalized_context_generators(context),
            chain=_regular_chain_from_lex_basis(context) or tuple(),
            dimension=context.dimension,
            degree=ideal_degree(context.generators, vars_),
            squarefree_regular=_squarefree_regular_chain_certified(context),
        )
        for context in sorted(
            reduced,
            key=lambda item: (-item.dimension, tuple(map(sp.default_sort_key, item.generators))),
        )
    )
    complete = all(component.squarefree_regular for component in components)

    # Reconstruct the reduced union independently.  This matters for inputs
    # with multiplicity: degree(I) need not equal degree(sqrt(I)).
    intersection = _intersection_all([component.equations for component in components], vars_)
    reconstruction_complete = False
    reduced_degree = None
    degree_complete = False
    if intersection is not None:
        reconstructed = EqualityIdealContext.build(intersection, vars_)
        source_in_reconstructed = all(
            reconstructed.normal_form(generator) == 0 for generator in source.generators
        )
        recon_in_radical = all(source.in_radical(generator) for generator in intersection)
        reconstruction_complete = source_in_reconstructed and recon_in_radical
        if reconstruction_complete:
            reduced_degree = ideal_degree(intersection, vars_)
            top_dimension = reconstructed.dimension
            top_degree = sum(
                component.degree for component in components if component.dimension == top_dimension
            )
            degree_complete = complete and top_degree == reduced_degree
    certificate = RegularChainDecompositionCertificate(
        variables=vars_,
        source_equations=_normalized_context_generators(source),
        splits=tuple(split_certificates),
        components=components,
        complete=complete,
        reconstruction_complete=reconstruction_complete,
        degree_complete=degree_complete,
        source_dimension=source.dimension,
        source_degree=source_degree,
        reduced_degree=reduced_degree,
    )
    return RegularChainDecomposition(
        variables=vars_,
        components=components,
        complete=complete,
        reconstruction_complete=reconstruction_complete,
        degree_complete=degree_complete,
        source_dimension=source.dimension,
        source_degree=source_degree,
        reduced_degree=reduced_degree,
        methods=tuple(sorted(methods)),
        certificate=certificate,
    )


def _same_canonical_ideal(
    left: Sequence[sp.Expr], right: Sequence[sp.Expr], variables: Sequence[sp.Symbol]
) -> bool:
    try:
        return _canonical_ideal_tuple(left, variables) == _canonical_ideal_tuple(right, variables)
    except RationalUnivariateError:
        return False


def _verify_radical_union(
    parent: EqualityIdealContext,
    children: Sequence[Sequence[sp.Expr]],
) -> bool:
    intersection = _intersection_all(children, parent.variables)
    if intersection is None:
        return False
    try:
        union_ideal = EqualityIdealContext.build(intersection, parent.variables)
    except RationalUnivariateError:
        return False
    parent_in_union = all(
        union_ideal.normal_form(generator) == 0 for generator in parent.generators
    )
    union_in_parent_radical = all(parent.in_radical(generator) for generator in intersection)
    return parent_in_union and union_in_parent_radical


def verify_regular_chain_decomposition_certificate(
    certificate: RegularChainDecompositionCertificate,
) -> bool:
    """Independently replay a recursive regular-chain proof tree.

    No decomposition search is rerun. Branch identities, terminal regular-chain
    claims, radical reconstruction and Hilbert-degree accounting are checked
    directly from the recorded algebraic data.
    """
    vars_ = certificate.variables
    try:
        source = EqualityIdealContext.build(certificate.source_equations, vars_)
    except (RationalUnivariateError, ValueError, TypeError):
        return False
    if _normalized_context_generators(source) != certificate.source_equations:
        return False
    if source.dimension != certificate.source_dimension:
        return False
    if ideal_degree(source.generators, vars_) != certificate.source_degree:
        return False

    active: list[tuple[sp.Expr, ...]] = [certificate.source_equations]
    for branch in certificate.splits:
        try:
            parent_key = _canonical_ideal_tuple(branch.parent, vars_)
        except RationalUnivariateError:
            return False
        if parent_key not in active or not branch.children:
            return False
        parent = EqualityIdealContext.build(parent_key, vars_)
        child_keys = tuple(_canonical_ideal_tuple(child, vars_) for child in branch.children)
        # Every child is a sublocus of the parent.
        for child_key in child_keys:
            child = EqualityIdealContext.build(child_key, vars_)
            if any(child.normal_form(generator) != 0 for generator in parent.generators):
                return False
        if branch.splitter is not None:
            split = _split_by_certified_splitter(parent, branch.splitter)
            if split is None:
                return False
            expected = {_canonical_ideal_tuple(item, vars_) for item in split}
            if expected != set(child_keys):
                return False
        elif not _verify_radical_union(parent, child_keys):
            return False
        active.remove(parent_key)
        active.extend(child_keys)

    final_keys: list[tuple[sp.Expr, ...]] = []
    for component in certificate.components:
        try:
            context = EqualityIdealContext.build(component.equations, vars_)
        except RationalUnivariateError:
            return False
        key = _normalized_context_generators(context)
        final_keys.append(key)
        chain = _regular_chain_from_lex_basis(context) or tuple()
        if chain != component.chain:
            return False
        squarefree = _squarefree_regular_chain_certified(context)
        if squarefree != component.squarefree_regular:
            return False
        if context.dimension != component.dimension:
            return False
        if ideal_degree(context.generators, vars_) != component.degree:
            return False

    # Cleanup may remove leaves contained in another leaf, but may not invent
    # a new terminal component or change the represented union.
    if any(key not in active for key in final_keys):
        return False
    for leaf in active:
        leaf_context = EqualityIdealContext.build(leaf, vars_)
        if not any(
            _variety_subset(leaf_context, EqualityIdealContext.build(key, vars_))
            for key in final_keys
        ):
            return False

    complete = all(component.squarefree_regular for component in certificate.components)
    if complete != certificate.complete:
        return False
    intersection = _intersection_all(final_keys, vars_)
    reconstruction_complete = False
    reduced_degree = None
    degree_complete = False
    if intersection is not None:
        reconstructed = EqualityIdealContext.build(intersection, vars_)
        reconstruction_complete = all(
            reconstructed.normal_form(generator) == 0 for generator in source.generators
        ) and all(source.in_radical(generator) for generator in intersection)
        if reconstruction_complete:
            reduced_degree = ideal_degree(intersection, vars_)
            top_dimension = reconstructed.dimension
            top_degree = sum(
                component.degree
                for component in certificate.components
                if component.dimension == top_dimension
            )
            degree_complete = complete and top_degree == reduced_degree
    return (
        reconstruction_complete == certificate.reconstruction_complete
        and reduced_degree == certificate.reduced_degree
        and degree_complete == certificate.degree_complete
    )


def verify_radical_ideal_certificate(certificate: RadicalIdealCertificate) -> bool:
    """Verify a radical computation without rerunning its search algorithm."""
    if not verify_regular_chain_decomposition_certificate(certificate.decomposition):
        return False
    if (
        not certificate.decomposition.complete
        or not certificate.decomposition.reconstruction_complete
    ):
        return False
    vars_ = certificate.variables
    if vars_ != certificate.decomposition.variables:
        return False
    if certificate.source_equations != certificate.decomposition.source_equations:
        return False
    source = EqualityIdealContext.build(certificate.source_equations, vars_)
    candidate = EqualityIdealContext.build(certificate.generators, vars_)
    if _normalized_context_generators(candidate) != tuple(certificate.generators):
        return False
    intersection = _intersection_all(
        [component.equations for component in certificate.decomposition.components], vars_
    )
    if intersection is None or not _same_canonical_ideal(
        intersection, certificate.generators, vars_
    ):
        return False
    return all(candidate.normal_form(generator) == 0 for generator in source.generators) and all(
        source.in_radical(generator) for generator in certificate.generators
    )


def verify_triangular_primality_certificate(
    certificate: TriangularPrimalityCertificate,
) -> bool:
    """Replay successive function-field primality/factorization evidence."""
    try:
        fresh = certify_triangular_primality(certificate.chain, certificate.variables)
    except (FunctionFieldError, RationalUnivariateError, ValueError, TypeError):
        return False
    return fresh == certificate


def verify_minimal_prime_decomposition_certificate(
    certificate: MinimalPrimeDecompositionCertificate,
) -> bool:
    """Verify radical reconstruction, primality, irredundancy and degree accounting."""
    vars_ = certificate.variables
    source = EqualityIdealContext.build(certificate.source_equations, vars_)
    contexts = [
        EqualityIdealContext.build(component.equations, vars_)
        for component in certificate.components
    ]
    all_radical = all(_radicality_certified(context) for context in contexts)
    intersection = _intersection_all(
        [component.equations for component in certificate.components], vars_
    )
    radical_complete = False
    if all_radical and intersection is not None:
        inter = EqualityIdealContext.build(intersection, vars_)
        radical_complete = all(
            inter.normal_form(generator) == 0 for generator in source.generators
        ) and all(source.in_radical(generator) for generator in intersection)
    if radical_complete != certificate.radical_complete:
        return False

    if len(certificate.primality_certificates) != len(certificate.components):
        return False
    if any(
        item is not None and not isinstance(item, TriangularPrimalityCertificate)
        for item in certificate.primality_certificates
    ):
        return False

    all_prime = True
    for index, (component, context) in enumerate(
        zip(certificate.components, contexts, strict=True)
    ):
        if (
            context.dimension != component.dimension
            or ideal_degree(context.generators, vars_) != component.degree
        ):
            return False
        method = _prime_certificate_method(context)
        if (method is not None) != component.prime or method != component.prime_method:
            return False
        if component.radical != _radicality_certified(context):
            return False
        tri = certificate.primality_certificates[index]
        if tri is not None and not verify_triangular_primality_certificate(tri):
            return False
        all_prime &= method is not None
    irredundant = not any(
        i != j and _variety_subset(left, right)
        for i, left in enumerate(contexts)
        for j, right in enumerate(contexts)
    )
    minimal_complete = radical_complete and all_prime and irredundant
    if minimal_complete != certificate.minimal_primes_complete:
        return False
    degree_complete = False
    if radical_complete and intersection is not None:
        radical_context = EqualityIdealContext.build(intersection, vars_)
        radical_degree = ideal_degree(intersection, vars_)
        top_dimension = radical_context.dimension
        degree_complete = (
            sum(
                component.degree
                for component in certificate.components
                if component.dimension == top_dimension
            )
            == radical_degree
        )
    return degree_complete == certificate.degree_complete


def radical_ideal(
    equations: Iterable[sp.Expr | sp.Equality],
    variables: Sequence[sp.Symbol],
    *,
    max_pieces: int | None = None,
) -> RadicalIdealResult:
    """Compute ``sqrt(<equations>)`` by certified recursive regular chains.

    With the default ``max_pieces=None`` there is no artificial branch budget.
    Every terminal component must be a certified squarefree regular-chain
    saturation.  Their exact ideal intersection is then checked in both
    directions against the source radical before generators are returned.
    """
    vars_ = tuple(variables)
    source = EqualityIdealContext.build(_normalize_equations(equations, vars_), vars_)
    decomposition = recursive_regular_chain_decomposition(
        source.generators, vars_, max_pieces=max_pieces
    )
    if not decomposition.complete or not decomposition.reconstruction_complete:
        return RadicalIdealResult(
            vars_, tuple(), False, "recursive_regular_chain_incomplete", decomposition
        )
    intersection = _intersection_all(
        [component.equations for component in decomposition.components], vars_
    )
    if intersection is None:
        return RadicalIdealResult(
            vars_, tuple(), False, "intersection_reconstruction_failed", decomposition
        )
    candidate = EqualityIdealContext.build(intersection, vars_)
    generators = _normalized_context_generators(candidate)
    source_in_candidate = all(
        candidate.normal_form(generator) == 0 for generator in source.generators
    )
    cand_in_radical = all(source.in_radical(generator) for generator in generators)
    complete = source_in_candidate and cand_in_radical and _radicality_certified(candidate)
    radical_certificate = None
    if complete and decomposition.certificate is not None:
        radical_certificate = RadicalIdealCertificate(
            variables=vars_,
            source_equations=_normalized_context_generators(source),
            generators=generators,
            decomposition=decomposition.certificate,
        )
    return RadicalIdealResult(
        vars_,
        generators if complete else tuple(),
        complete,
        "certified_squarefree_regular_chain_intersection" if complete else "verification_failed",
        decomposition,
        radical_certificate,
    )


def certified_radical_minimal_prime_decomposition(
    equations: Iterable[sp.Expr | sp.Equality],
    variables: Sequence[sp.Symbol],
    *,
    max_pieces: int | None = None,
) -> RadicalMinimalPrimeDecomposition:
    """Certify ``sqrt(I)`` and, when possible, all minimal primes of ``I``.

    The geometric splitting engine supplies an exact finite cover.  Each branch
    must independently pass a radicality certificate.  The intersection of the
    branch ideals is then reconstructed by elimination and checked in both
    directions against ``sqrt(I)`` using exact ideal/radical membership.

    Minimal-prime completeness is stronger: every surviving radical branch must
    have an exact primality certificate and the prime cover must be pairwise
    irredundant.  Unsupported primality cases remain useful radical components
    but set ``minimal_primes_complete=False``.
    """
    vars_ = tuple(variables)
    source = EqualityIdealContext.build(_normalize_equations(equations, vars_), vars_)
    decomposition = _equidimensional_decomposition_impl(
        source.generators,
        vars_,
        max_pieces=128 if max_pieces is None else max_pieces,
        _emit_certificate=True,
        _algebraic_semantics=True,
    )
    contexts = [
        EqualityIdealContext.build(piece.equations, vars_) for piece in decomposition.pieces
    ]
    contexts = list(
        _recursive_regular_chain_prime_refinement(
            contexts, vars_, max_pieces=128 if max_pieces is None else max_pieces
        )
    )
    components: list[RadicalComponent] = []
    all_radical = True
    all_prime = True
    for context in contexts:
        radical = _radicality_certified(context)
        prime_method = _prime_certificate_method(context) if radical else None
        prime = prime_method is not None
        all_radical &= radical
        all_prime &= prime
        components.append(
            RadicalComponent(
                _normalized_context_generators(context),
                context.dimension,
                radical,
                prime,
                ideal_degree(context.generators, vars_),
                prime_method,
            )
        )

    intersection = _intersection_all([component.equations for component in components], vars_)
    radical_complete = False
    if all_radical and intersection is not None:
        intersection_context = EqualityIdealContext.build(intersection, vars_)
        source_in_intersection = all(
            intersection_context.normal_form(generator) == 0 for generator in source.generators
        )
        inter_in_radical = all(source.in_radical(generator) for generator in intersection)
        radical_complete = source_in_intersection and inter_in_radical

    irredundant = True
    for i, left in enumerate(contexts):
        for j, right in enumerate(contexts):
            if i != j and _variety_subset(left, right):
                irredundant = False
                break
        if not irredundant:
            break
    minimal_complete = radical_complete and all_prime and irredundant
    degree_complete = False
    if radical_complete and intersection is not None:
        radical_context = EqualityIdealContext.build(intersection, vars_)
        radical_degree = ideal_degree(intersection, vars_)
        top_dimension = radical_context.dimension
        component_degree = sum(
            component.degree for component in components if component.dimension == top_dimension
        )
        degree_complete = component_degree == radical_degree
    component_tuple = tuple(components)
    primality_certificates: list[TriangularPrimalityCertificate | None] = []
    for component, context in zip(component_tuple, contexts, strict=True):
        if component.prime_method and (
            "triangular" in component.prime_method or "function_field" in component.prime_method
        ):
            tri = certify_triangular_primality(context.generators, vars_)
            primality_certificates.append(tri if tri.complete else None)
        else:
            primality_certificates.append(None)
    certificate = MinimalPrimeDecompositionCertificate(
        variables=vars_,
        source_equations=_normalized_context_generators(source),
        components=component_tuple,
        primality_certificates=tuple(primality_certificates),
        radical_complete=radical_complete,
        minimal_primes_complete=minimal_complete,
        degree_complete=degree_complete,
    )
    return RadicalMinimalPrimeDecomposition(
        variables=vars_,
        components=component_tuple,
        radical_complete=radical_complete,
        minimal_primes_complete=minimal_complete,
        degree_complete=degree_complete,
        certificate=certificate,
    )


def _is_zero_dimensional_maximal_prime(
    generators: Sequence[sp.Expr], variables: Sequence[sp.Symbol]
) -> bool:
    context = EqualityIdealContext.build(generators, variables)
    if context.dimension != 0:
        return False
    method = _prime_certificate_method(context)
    return method is not None


def _radical_relation_certified(ideal: EqualityIdealContext, radical: EqualityIdealContext) -> bool:
    """Certify ``radical.generators`` generate exactly ``sqrt(ideal)``.

    If R is radical, I subset R and every generator of R lies in sqrt(I), then
    sqrt(I)=R.  This is the ideal-theoretic relation needed by every primary
    component certificate, independent of the decomposition algorithm used.
    """

    if not _radicality_certified(radical):
        return False
    ideal_in_radical = all(radical.normal_form(g) == 0 for g in ideal.generators)
    radical_in_sqrt = all(ideal.in_radical(g) for g in radical.generators)
    return ideal_in_radical and radical_in_sqrt


def _primary_component_certificate(
    generators: Sequence[sp.Expr],
    radical_generators: Sequence[sp.Expr],
    variables: Sequence[sp.Symbol],
) -> tuple[bool, str]:
    """Return a conservative exact primary certificate for ``Q``."""
    vars_ = tuple(variables)
    qctx = EqualityIdealContext.build(generators, vars_)
    rctx = EqualityIdealContext.build(radical_generators, vars_)

    # Every prime ideal is primary, including the zero ideal in a polynomial
    # ring over a field.
    if _ideal_generators_equal(qctx.generators, rctx.generators, vars_):
        if _prime_certificate_method(rctx) is not None:
            return True, "prime_is_primary"

    # UFD theorem: if p is irreducible, (p) is prime and (p**e) is
    # (p)-primary.  Check this before invoking the general radical engine.
    if len(rctx.generators) == 1 and len(qctx.generators) == 1:
        try:
            p = sp.Poly(rctx.generators[0], *vars_, domain=sp.QQ)
            q = sp.Poly(qctx.generators[0], *vars_, domain=sp.QQ)
            _unit, factors = sp.factor_list(p.as_expr())
        except (sp.PolynomialError, ValueError):
            factors = []
        if len(factors) == 1 and factors[0][1] == 1:
            base = sp.Poly(factors[0][0], *vars_, domain=sp.QQ).monic()
            for exponent in range(1, max(2, q.total_degree() + 1)):
                candidate = sp.Poly(base.as_expr() ** exponent, *vars_, domain=sp.QQ)
                if q.monic() == candidate.monic():
                    return True, "principal_prime_power"

    # Irreducible monomial ideals <x_i**a_i> are primary, with radical
    # generated by exactly those variables.
    q_basis = qctx.groebner_basis
    if q_basis is not None:
        powers: dict[sp.Symbol, int] = {}
        monomial_primary = True
        for poly in q_basis.polys:
            terms = poly.terms()
            if len(terms) != 1:
                monomial_primary = False
                break
            exponent, _coeff = terms[0]
            support = [i for i, e in enumerate(exponent) if e]
            if len(support) != 1:
                monomial_primary = False
                break
            i = support[0]
            powers[vars_[i]] = exponent[i]
        if monomial_primary:
            expected_radical = tuple(powers)
            if _ideal_generators_equal(expected_radical, rctx.generators, vars_):
                return True, "irreducible_monomial_primary"

    # For zero-dimensional ideals use the independently certified reduced
    # decomposition: a single certified minimal prime P proves sqrt(Q)=P.
    if qctx.dimension == 0 and _prime_certificate_method(rctx) is not None:
        reduced = certified_radical_minimal_prime_decomposition(qctx.generators, vars_)
        if (
            reduced.radical_complete
            and reduced.minimal_primes_complete
            and len(reduced.components) == 1
            and _ideal_generators_equal(reduced.components[0].equations, rctx.generators, vars_)
        ):
            return True, "zero_dimensional_maximal_radical"
        return False, "radical_mismatch"
    return False, "unsupported_primary_certificate"


def _principal_primary_decomposition(
    source: EqualityIdealContext, variables: Sequence[sp.Symbol]
) -> tuple[PrimaryComponent, ...] | None:
    vars_ = tuple(variables)
    generators = _normalized_context_generators(source)
    if len(generators) != 1 or generators[0] in (0, 1):
        return None
    try:
        _unit, factors = sp.factor_list(generators[0], *vars_)
    except (sp.PolynomialError, ValueError, TypeError, NotImplementedError):
        return None
    if not factors:
        return None
    out = []
    for factor, exponent in factors:
        pctx = EqualityIdealContext.build((sp.expand(factor),), vars_)
        if _prime_certificate_method(pctx) is None:
            return None
        qctx = EqualityIdealContext.build((sp.expand(factor**exponent),), vars_)
        ok, method = _primary_component_certificate(qctx.generators, pctx.generators, vars_)
        if not ok:
            return None
        out.append(
            PrimaryComponent(
                _normalized_context_generators(qctx),
                _normalized_context_generators(pctx),
                qctx.dimension,
                ideal_degree(qctx.generators, vars_),
                True,
                method,
            )
        )
    return tuple(out)


def _monomial_primary_decomposition(
    source: EqualityIdealContext, variables: Sequence[sp.Symbol]
) -> tuple[PrimaryComponent, ...] | None:
    """Return the canonical irreducible decomposition of a monomial ideal.

    The ideal is Artinianized by adding artificial powers x_i**c_i with c_i
    larger than every exponent already occurring. Maximal standard monomials
    of that finite box encode the irreducible monomial components; artificial
    powers are omitted again on de-Artinianization.
    """
    vars_ = tuple(variables)
    basis = source.groebner_basis
    if basis is None or source.inconsistent:
        return None
    exponents: list[tuple[int, ...]] = []
    for poly in basis.polys:
        terms = poly.terms()
        if len(terms) != 1:
            return None
        exponent, _coeff = terms[0]
        exponents.append(tuple(int(e) for e in exponent))
    if not exponents:
        return None
    bounds = tuple(max(exp[i] for exp in exponents) + 1 for i in range(len(vars_)))

    from itertools import product

    def in_ideal(monomial: tuple[int, ...]) -> bool:
        return any(all(monomial[i] >= gen[i] for i in range(len(vars_))) for gen in exponents)

    maximal_standard: list[tuple[int, ...]] = []
    for u in product(*(range(bound) for bound in bounds)):
        if in_ideal(u):
            continue
        maximal = True
        for i, bound in enumerate(bounds):
            if u[i] + 1 >= bound:
                continue
            v = list(u)
            v[i] += 1
            if not in_ideal(tuple(v)):
                maximal = False
                break
        # Reaching an artificial boundary also counts as blocked: in the
        # Artinianized ideal x_i**bound is present.
        if maximal:
            maximal_standard.append(tuple(u))

    components: list[PrimaryComponent] = []
    seen: set[tuple[sp.Expr, ...]] = set()
    for u in maximal_standard:
        qgens = tuple(
            sp.expand(var ** (u[i] + 1)) for i, var in enumerate(vars_) if u[i] + 1 < bounds[i]
        )
        radical = tuple(var for i, var in enumerate(vars_) if u[i] + 1 < bounds[i])
        qctx = EqualityIdealContext.build(qgens, vars_)
        rctx = EqualityIdealContext.build(radical, vars_)
        key = _normalized_context_generators(qctx)
        if key in seen:
            continue
        seen.add(key)
        ok, method = _primary_component_certificate(qctx.generators, rctx.generators, vars_)
        if not ok:
            return None
        components.append(
            PrimaryComponent(
                key,
                _normalized_context_generators(rctx),
                qctx.dimension,
                ideal_degree(qctx.generators, vars_),
                True,
                method,
            )
        )
    if not components:
        return None
    intersection = _intersection_all([c.equations for c in components], vars_)
    if intersection is None or not _ideal_generators_equal(intersection, source.generators, vars_):
        return None
    return tuple(components)


def _separator_for_prime(index: int, primes: Sequence[EqualityIdealContext]) -> sp.Expr | None:
    """Choose an exact separator vanishing on every other maximal prime."""
    target = primes[index]
    factors = []
    for j, other in enumerate(primes):
        if j == index:
            continue
        chosen = None
        for generator in other.generators:
            if target.normal_form(generator) != 0:
                chosen = generator
                break
        if chosen is None:
            return None
        factors.append(chosen)
    return sp.expand(sp.prod(factors)) if factors else sp.Integer(1)


def _zero_dimensional_primary_decomposition(
    source: EqualityIdealContext, variables: Sequence[sp.Symbol]
) -> tuple[PrimaryComponent, ...] | None:
    vars_ = tuple(variables)
    if source.dimension != 0:
        return None
    minimal = certified_radical_minimal_prime_decomposition(source.generators, vars_)
    if not minimal.minimal_primes_complete:
        return None
    primes = [EqualityIdealContext.build(c.equations, vars_) for c in minimal.components]
    if not all(p.dimension == 0 and _prime_certificate_method(p) is not None for p in primes):
        return None
    components = []
    for i, prime in enumerate(primes):
        separator = _separator_for_prime(i, primes)
        if separator is None:
            return None
        saturated = _saturate_generators(source.generators, separator, vars_)
        if saturated is None:
            return None
        qctx = EqualityIdealContext.build(saturated, vars_)
        ok, method = _primary_component_certificate(qctx.generators, prime.generators, vars_)
        if not ok:
            return None
        components.append(
            PrimaryComponent(
                _normalized_context_generators(qctx),
                _normalized_context_generators(prime),
                qctx.dimension,
                ideal_degree(qctx.generators, vars_),
                True,
                method,
            )
        )
    return tuple(components)


def _primary_components_irredundant(
    components: Sequence[PrimaryComponent], variables: Sequence[sp.Symbol]
) -> bool:
    vars_ = tuple(variables)
    if len({tuple(c.radical) for c in components}) != len(components):
        return False
    for i, component in enumerate(components):
        others = [c.equations for j, c in enumerate(components) if j != i]
        if not others:
            continue
        intersection = _intersection_all(others, vars_)
        if intersection is None:
            return False
        # Component i is redundant iff intersection(other components) is
        # already contained in Q_i.
        other_ctx = EqualityIdealContext.build(intersection, vars_)
        qctx = EqualityIdealContext.build(component.equations, vars_)
        if all(qctx.normal_form(g) == 0 for g in other_ctx.generators):
            return False
    return True


def verify_primary_decomposition_certificate(
    certificate: PrimaryDecompositionCertificate,
) -> bool:
    """Replay a primary decomposition using exact ideal arithmetic only."""
    try:
        from .algebraic.gtz_primary import (
            GTZPrimaryDecompositionCertificate,
            verify_gtz_primary_decomposition_certificate,
        )

        if isinstance(certificate, GTZPrimaryDecompositionCertificate):
            return verify_gtz_primary_decomposition_certificate(certificate)
    except ImportError:
        pass
    vars_ = certificate.variables
    source = EqualityIdealContext.build(certificate.source_equations, vars_)
    if not certificate.components:
        return source.inconsistent and certificate.complete
    for component in certificate.components:
        if not component.primary:
            return False
        qctx = EqualityIdealContext.build(component.equations, vars_)
        rctx = EqualityIdealContext.build(component.radical, vars_)
        if not _radical_relation_certified(qctx, rctx):
            return False
        ok, method = _primary_component_certificate(qctx.generators, rctx.generators, vars_)
        if not ok or method != component.method:
            return False
        if qctx.dimension != component.dimension:
            return False
        if ideal_degree(qctx.generators, vars_) != component.degree:
            return False
    intersection = _intersection_all([c.equations for c in certificate.components], vars_)
    if intersection is None or not _ideal_generators_equal(intersection, source.generators, vars_):
        return False
    irredundant = _primary_components_irredundant(certificate.components, vars_)
    return certificate.complete and irredundant == certificate.irredundant


def primary_decomposition(
    equations: Iterable[sp.Expr | sp.Equality],
    variables: Sequence[sp.Symbol],
) -> PrimaryDecompositionResult:
    """Compute a certified primary decomposition when the proof is supported.

    The implementation is complete for certified prime ideals, principal
    ideals over ``QQ``, arbitrary monomial ideals, and zero-dimensional ideals
    whose minimal primes can be certified.  These regimes cover arbitrary
    multiplicities and monomial embedded primes.  Positive-dimensional
    nonmonomial positive-dimensional cases fall through to the full certified
    recursive GTZ localization/contraction driver.
    """
    vars_ = tuple(variables)
    source = EqualityIdealContext.build(_normalize_equations(equations, vars_), vars_)
    if source.inconsistent:
        components: tuple[PrimaryComponent, ...] = tuple()
        method = "unit_ideal"
    else:
        prime_method = _prime_certificate_method(source)
        if prime_method is not None:
            generators = _normalized_context_generators(source)
            components = (
                PrimaryComponent(
                    generators,
                    generators,
                    source.dimension,
                    ideal_degree(source.generators, vars_),
                    True,
                    "prime_is_primary",
                ),
            )
            method = "prime_ideal"
        else:
            components = _principal_primary_decomposition(source, vars_) or tuple()
            method = "principal_factorization" if components else ""
        if not components:
            monomial = _monomial_primary_decomposition(source, vars_)
            if monomial is not None:
                components = monomial
                method = "monomial_irreducible_decomposition"
        if not components:
            zero_dim = _zero_dimensional_primary_decomposition(source, vars_)
            if zero_dim is not None:
                components = zero_dim
                method = "zero_dimensional_localization"
    if not components and not source.inconsistent:
        from .algebraic.gtz_primary import gtz_primary_decomposition

        gtz = gtz_primary_decomposition(source.generators, vars_)
        if gtz.complete and gtz.certificate is not None:
            return PrimaryDecompositionResult(
                vars_,
                tuple(gtz.components),
                True,
                gtz.irredundant,
                gtz.method,
                gtz.certificate,
            )
        return PrimaryDecompositionResult(vars_, tuple(), False, False, "gtz_incomplete")
    intersection = _intersection_all([c.equations for c in components], vars_)
    complete = source.inconsistent or (
        intersection is not None and _ideal_generators_equal(intersection, source.generators, vars_)
    )
    irredundant = _primary_components_irredundant(components, vars_) if components else True
    cert = None
    if complete:
        cert = PrimaryDecompositionCertificate(
            vars_, _normalized_context_generators(source), components, True, irredundant
        )
    return PrimaryDecompositionResult(vars_, components, complete, irredundant, method, cert)


def associated_primes(
    equations: Iterable[sp.Expr | sp.Equality],
    variables: Sequence[sp.Symbol],
) -> AssociatedPrimesResult:
    """Return all associated primes when a certified primary decomposition exists."""
    decomposition = primary_decomposition(equations, variables)
    primes = tuple(component.radical for component in decomposition.components)
    complete = decomposition.complete and decomposition.irredundant
    return AssociatedPrimesResult(
        tuple(variables),
        primes if complete else tuple(),
        complete,
        "certified_primary_decomposition" if complete else "primary_decomposition_incomplete",
        decomposition,
    )


def _equidimensional_decomposition_impl(
    equations: Iterable[sp.Expr | sp.Equality],
    variables: Sequence[sp.Symbol],
    *,
    max_pieces: int = 128,
    _emit_certificate: bool = True,
    _algebraic_semantics: bool = False,
) -> EquidimensionalDecomposition:
    """Decompose a polynomial zero set into exact dimension-labelled pieces.

    The regular-chain refinement follows the triangular/unmixed decomposition
    tradition of [Kalkbrener1993] and [Lazard1991]; see ``docs/references.md``.

    The algorithm has reduced-set semantics. Monomial Groebner bases are
    decomposed exactly through the minimal coordinate primes of their radical.
    More generally, exactly factorizable Groebner generators are decomposed by
    minimal hitting sets of their square-free factor supports. This extends the
    certified reduced-set decomposition to nonmonomial factor-incidence ideals
    and removes redundant/embedded factor branches without branch guessing.
    Radical-containment tests remove any remaining duplicate and contained
    branches. The returned union is exact.

    General primary decomposition is not guessed. When terminal
    ideals cannot be proved equidimensional by the supported exact criteria,
    their pieces are retained with ``equidimensional=False`` and the overall
    result has ``complete=False``.  This conservative contract lets downstream
    singularity analysis use only certified decompositions.
    """
    vars_ = tuple(variables)
    if not vars_:
        raise ValueError("equidimensional decomposition requires ambient variables")
    if max_pieces < 1:
        raise ValueError("max_pieces must be positive")

    initial = _normalize_equations(equations, vars_)
    methods_used: set[str] = set()
    # The optional integer carried with a pending branch records the Krull
    # dimension before a real-radical rewrite. Algebraic splitting preserves
    # it so each final real piece can report both dimensions faithfully.
    pending: list[tuple[tuple[sp.Expr, ...], int | None]] = [(initial, None)]
    terminal: list[tuple[EqualityIdealContext, int | None]] = []
    seen_basis: set[tuple[sp.Expr, ...]] = set()

    while pending:
        generators, algebraic_dim_override = pending.pop()
        try:
            context = EqualityIdealContext.build(generators, vars_)
        except RationalUnivariateError as exc:
            raise NotImplementedError(
                "equidimensional decomposition requires exact polynomial coefficients"
            ) from exc
        if context.inconsistent:
            continue
        basis = (
            tuple(sp.expand(poly.as_expr()) for poly in context.groebner_basis.polys)
            if context.groebner_basis is not None
            else tuple()
        )
        if basis in seen_basis:
            continue
        seen_basis.add(basis)

        real_reduction = (
            None if _algebraic_semantics else _real_radical_reduction_generators(context)
        )
        if real_reduction is not None:
            methods_used.add("real_radical_reduction")
            if len(pending) + len(terminal) + 1 > max_pieces:
                terminal.append((context, algebraic_dim_override))
                continue
            source_dimension = (
                context.dimension if algebraic_dim_override is None else algebraic_dim_override
            )
            pending.append((real_reduction, source_dimension))
            continue

        monomial_primes = _monomial_minimal_prime_generators(context)
        if monomial_primes is not None:
            methods_used.add("monomial_minimal_primes")
            current = frozenset(basis)
            prime_sets = tuple(
                frozenset(sp.expand(eq) for eq in prime) for prime in monomial_primes
            )
            already_minimal_prime = len(prime_sets) == 1 and prime_sets[0] == current
            if not already_minimal_prime:
                if len(pending) + len(terminal) + len(monomial_primes) > max_pieces:
                    terminal.append((context, algebraic_dim_override))
                    continue
                pending.extend((prime, algebraic_dim_override) for prime in monomial_primes)
                continue

        factor_components = _factor_incidence_component_generators(context)
        if factor_components is not None:
            methods_used.add("factor_incidence")
            if len(pending) + len(terminal) + len(factor_components) > max_pieces:
                terminal.append((context, algebraic_dim_override))
                continue
            pending.extend((component, algebraic_dim_override) for component in factor_components)
            continue

        triangular_basis = _lex_basis_generators(context)
        initial_split = _initial_split_generators(context, basis=triangular_basis)
        if initial_split is not None:
            methods_used.add("initial_split")
            if len(pending) + len(terminal) + len(initial_split) > max_pieces:
                terminal.append((context, algebraic_dim_override))
                continue
            pending.extend((branch, algebraic_dim_override) for branch in initial_split)
            continue

        separant_split = _separant_split_generators(context, basis=triangular_basis)
        if separant_split is not None:
            methods_used.add("separant_split")
            if len(pending) + len(terminal) + len(separant_split) > max_pieces:
                terminal.append((context, algebraic_dim_override))
                continue
            pending.extend((branch, algebraic_dim_override) for branch in separant_split)
            continue

        subresultant_split = _regular_gcd_split_generators(context, basis=triangular_basis)
        if subresultant_split is not None:
            methods_used.add("regular_subresultant_split")
            if len(pending) + len(terminal) + len(subresultant_split) > max_pieces:
                terminal.append((context, algebraic_dim_override))
                continue
            pending.extend((branch, algebraic_dim_override) for branch in subresultant_split)
            continue

        saturation_split = _saturation_split_generators(context)
        if saturation_split is not None:
            methods_used.add("saturation_split")
            if len(pending) + len(terminal) + len(saturation_split) > max_pieces:
                terminal.append((context, algebraic_dim_override))
                continue
            pending.extend((branch, algebraic_dim_override) for branch in saturation_split)
            continue

        split = None
        for polynomial in basis:
            split = _distinct_factors(polynomial, vars_)
            if split is not None:
                break
        if split is None:
            terminal.append((context, algebraic_dim_override))
            continue
        methods_used.add("factor_split")
        if len(pending) + len(terminal) + len(split) > max_pieces:
            terminal.append((context, algebraic_dim_override))
            continue
        for factor in split:
            pending.append(((*context.generators, factor), algebraic_dim_override))

    # Remove branches whose real/complex algebraic zero set is contained in a
    # different surviving branch.  Equality is resolved deterministically by
    # keeping the structurally smaller Groebner representation.
    ordered = sorted(
        terminal,
        key=lambda item: (
            -item[0].dimension,
            len(item[0].generators),
            tuple(sp.default_sort_key(expr) for expr in item[0].generators),
        ),
    )
    kept: list[tuple[EqualityIdealContext, int | None]] = []
    for candidate, candidate_override in ordered:
        redundant = False
        for existing, _existing_override in kept:
            if _variety_subset(candidate, existing):
                redundant = True
                break
        if redundant:
            continue
        kept = [item for item in kept if not _variety_subset(item[0], candidate)]
        kept.append((candidate, candidate_override))

    real_contexts: list[tuple[EqualityIdealContext, int, int, bool]] = []
    for context, algebraic_dim_override in kept:
        # Algebraic ideal algorithms must retain components over the algebraic
        # closure, including components with no real points.  Public geometric
        # decomposition instead uses exact real dimension and may discard such
        # branches.  Keeping this distinction here also avoids an unnecessary
        # CAD construction in algebraic primary/minimal-prime decomposition.
        real_dimension = context.dimension if _algebraic_semantics else _real_dimension(context)
        if real_dimension < 0:
            continue
        if not _algebraic_semantics and real_dimension == 0 and context.dimension > 0:
            point_contexts = _real_point_contexts(context)
            if point_contexts:
                source_dimension = (
                    context.dimension if algebraic_dim_override is None else algebraic_dim_override
                )
                real_contexts.extend((point, 0, source_dimension, True) for point in point_contexts)
                continue
        certified = _equidimensionality_certified(context) and real_dimension == context.dimension
        if certified and real_dimension > 0:
            codimension = len(context.variables) - context.dimension
            basis_size = (
                len(context.groebner_basis.polys)
                if context.groebner_basis is not None
                else len(context.generators)
            )
            if basis_size != codimension and _squarefree_regular_chain_certified(context):
                methods_used.add("squarefree_regular_chain")
        if real_dimension == 0:
            certified = True
        source_dimension = (
            context.dimension if algebraic_dim_override is None else algebraic_dim_override
        )
        real_contexts.append((context, real_dimension, source_dimension, certified))

    pieces = tuple(
        AlgebraicLocusPiece(
            variables=vars_,
            equations=(
                tuple(sp.expand(poly.as_expr()) for poly in context.groebner_basis.polys)
                if context.groebner_basis is not None
                else tuple()
            ),
            dimension=real_dimension,
            algebraic_dimension=algebraic_dimension,
            equidimensional=certified,
        )
        for context, real_dimension, algebraic_dimension, certified in sorted(
            real_contexts,
            key=lambda item: (
                -item[1],
                tuple(sp.default_sort_key(expr) for expr in item[0].generators),
            ),
        )
    )
    complete = all(piece.equidimensional for piece in pieces)
    certificate = None
    if _emit_certificate:
        certified_pieces = tuple(
            DecompositionPieceCertificate(
                equations=piece.equations,
                dimension=piece.dimension,
                algebraic_dimension=piece.algebraic_dimension,
                equidimensional=piece.equidimensional,
            )
            for piece in pieces
        )
        methods = tuple(sorted(methods_used))
        certificate = DecompositionCertificate(
            variables=vars_,
            normalized_equations=initial,
            pieces=certified_pieces,
            complete=complete,
            methods=methods,
            max_pieces=max_pieces,
            digest=_certificate_digest(
                vars_, initial, certified_pieces, complete, methods, max_pieces
            ),
        )
    return EquidimensionalDecomposition(
        variables=vars_,
        pieces=pieces,
        complete=complete,
        certificate=certificate,
    )


def equidimensional_decomposition(
    equations: Iterable[sp.Expr | sp.Equality],
    variables: Sequence[sp.Symbol],
    *,
    max_pieces: int = 128,
    _emit_certificate: bool = True,
) -> EquidimensionalDecomposition:
    """Public wrapper around the certified reduced-locus decomposition core."""
    return _equidimensional_decomposition_impl(
        equations, variables, max_pieces=max_pieces, _emit_certificate=_emit_certificate
    )


def verify_decomposition_certificate(certificate: DecompositionCertificate) -> bool:
    """Replay and verify a decomposition certificate.

    Verification is strict: the normalized input is decomposed
    again with certificate emission disabled, and every canonical piece, real
    dimension, algebraic provenance dimension, equidimensionality flag, and
    completeness claim must agree. The digest additionally detects mutation of
    any recorded certificate field.
    """
    expected_digest = _certificate_digest(
        certificate.variables,
        certificate.normalized_equations,
        certificate.pieces,
        certificate.complete,
        certificate.methods,
        certificate.max_pieces,
    )
    if expected_digest != certificate.digest:
        return False
    try:
        replay = equidimensional_decomposition(
            certificate.normalized_equations,
            certificate.variables,
            max_pieces=certificate.max_pieces,
            _emit_certificate=True,
        )
    except (NotImplementedError, ValueError, RationalUnivariateError):
        return False
    replay_pieces = tuple(
        DecompositionPieceCertificate(
            equations=piece.equations,
            dimension=piece.dimension,
            algebraic_dimension=piece.algebraic_dimension,
            equidimensional=piece.equidimensional,
        )
        for piece in replay.pieces
    )
    replay_methods = replay.certificate.methods if replay.certificate is not None else tuple()
    return (
        replay.complete == certificate.complete
        and replay_pieces == certificate.pieces
        and replay_methods == certificate.methods
    )


__all__ = [
    "AlgebraicLocusPiece",
    "RegularChainSplitCertificate",
    "RegularChainDecompositionCertificate",
    "RadicalIdealCertificate",
    "MinimalPrimeDecompositionCertificate",
    "RadicalComponent",
    "RadicalMinimalPrimeDecomposition",
    "RadicalIdealResult",
    "DecompositionPieceCertificate",
    "DecompositionCertificate",
    "EquidimensionalDecomposition",
    "RegularChainComponent",
    "RegularChainDecomposition",
    "TriangularPrimalityStage",
    "TriangularPrimalityCertificate",
    "recursive_regular_chain_decomposition",
    "certify_triangular_primality",
    "radical_ideal",
    "certified_radical_minimal_prime_decomposition",
    "equidimensional_decomposition",
    "verify_regular_chain_decomposition_certificate",
    "verify_radical_ideal_certificate",
    "verify_triangular_primality_certificate",
    "verify_minimal_prime_decomposition_certificate",
    "PrimaryComponent",
    "PrimaryDecompositionCertificate",
    "PrimaryDecompositionResult",
    "AssociatedPrimesResult",
    "primary_decomposition",
    "associated_primes",
    "verify_primary_decomposition_certificate",
    "verify_decomposition_certificate",
]
