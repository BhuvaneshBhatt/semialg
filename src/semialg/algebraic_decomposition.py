"""Exact set-reduced decomposition of polynomial loci by dimension."""

from __future__ import annotations

from collections.abc import Iterable, Sequence
from dataclasses import dataclass
from hashlib import sha256
from itertools import combinations

import sympy as sp

from .algebraic.equality_ideal import EqualityIdealContext
from .algebraic.groebner_utils import lex_basis_generators
from .algebraic.rational_univariate.representation import RationalUnivariateError


@dataclass(frozen=True)
class AlgebraicLocusPiece:
    """One exact polynomial-locus piece with a certified global dimension.

    ``dimension`` is the exact real semialgebraic dimension;
    ``algebraic_dimension`` records the Krull dimension of the polynomial
    ideal before real-set reduction. ``equidimensional`` is true only when the
    current implementation can prove that the represented reduced real set has
    no component of a different dimension. Multiplicity in the input ideal has
    no geometric significance.
    """

    variables: tuple[sp.Symbol, ...]
    equations: tuple[sp.Expr, ...]
    dimension: int
    algebraic_dimension: int
    equidimensional: bool

    @property
    def codimension(self) -> int:
        """Return the codimension in the recorded ambient polynomial ring."""
        return len(self.variables) - self.dimension


@dataclass(frozen=True)
class DecompositionPieceCertificate:
    """Replay data for one reduced-real decomposition piece."""

    equations: tuple[sp.Expr, ...]
    dimension: int
    algebraic_dimension: int
    equidimensional: bool


@dataclass(frozen=True)
class DecompositionCertificate:
    """Replayable certificate for an equidimensional decomposition result.

    The certificate records the normalized input, canonical final pieces, the
    exact reduction/splitting mechanisms exercised, and a deterministic digest.
    :func:`verify_decomposition_certificate` replays the certified decomposition
    from the normalized input and checks that all recorded semantic claims match.
    """

    variables: tuple[sp.Symbol, ...]
    normalized_equations: tuple[sp.Expr, ...]
    pieces: tuple[DecompositionPieceCertificate, ...]
    complete: bool
    methods: tuple[str, ...]
    max_pieces: int
    digest: str


def _certificate_digest(
    variables: Sequence[sp.Symbol],
    normalized_equations: Sequence[sp.Expr],
    pieces: Sequence[DecompositionPieceCertificate],
    complete: bool,
    methods: Sequence[str],
    max_pieces: int,
) -> str:
    payload = repr(
        (
            tuple(map(str, variables)),
            tuple(sp.srepr(eq) for eq in normalized_equations),
            tuple(
                (
                    tuple(sp.srepr(eq) for eq in piece.equations),
                    piece.dimension,
                    piece.algebraic_dimension,
                    piece.equidimensional,
                )
                for piece in pieces
            ),
            bool(complete),
            tuple(methods),
            int(max_pieces),
        )
    ).encode("utf-8")
    return sha256(payload).hexdigest()


@dataclass(frozen=True)
class EquidimensionalDecomposition:
    """Exact finite cover of a reduced polynomial locus grouped by dimension.

    ``complete`` means every returned piece is certified equidimensional, so
    the result is suitable as the substrate for component-relative Jacobian
    regularity.  If it is false, the union is still exactly the same geometric
    zero set, but at least one piece may contain components of several local
    dimensions and callers must use a more complete method before reasoning
    componentwise.
    """

    variables: tuple[sp.Symbol, ...]
    pieces: tuple[AlgebraicLocusPiece, ...]
    complete: bool
    certificate: DecompositionCertificate | None = None

    @property
    def dimensions(self) -> tuple[int, ...]:
        return tuple(sorted({piece.dimension for piece in self.pieces}, reverse=True))

    def pieces_of_dimension(self, dimension: int) -> tuple[AlgebraicLocusPiece, ...]:
        """Return all pieces whose certified global dimension is ``dimension``."""
        return tuple(piece for piece in self.pieces if piece.dimension == dimension)


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


def _saturation_split_generators(
    context: EqualityIdealContext,
) -> tuple[tuple[sp.Expr, ...], tuple[sp.Expr, ...]] | None:
    """Find one certified proper split using ``V(I)=V(I+h) union V(I:h^inf)``.

    Candidate splitters are deliberately modest: ambient coordinates and exact
    affine-linear factors occurring in univariate Groebner-basis elements.  A
    candidate is accepted only after radical-containment checks prove that both
    resulting branches are strict subloci of the current reduced variety.
    This makes the procedure a proof-producing fallback rather than a
    primary-decomposition heuristic.
    """
    if context.inconsistent or not context.generators:
        return None

    candidates: list[sp.Expr] = list(context.variables)
    basis = _normalized_context_generators(context)
    for polynomial in basis:
        active = tuple(v for v in context.variables if polynomial.has(v))
        if len(active) != 1:
            continue
        factors = _squarefree_factor_atoms(polynomial, context.variables)
        if factors is None:
            continue
        for factor in factors:
            try:
                poly = sp.Poly(factor, *context.variables, extension=True)
            except (sp.PolynomialError, ValueError, TypeError, NotImplementedError):
                continue
            if poly.total_degree() == 1 and factor not in candidates:
                candidates.append(factor)

    for splitter in candidates:
        # If h already vanishes on all of V(I), I+<h> is not a strict branch.
        if context.in_radical(splitter):
            continue
        closed = EqualityIdealContext.build((*context.generators, splitter), context.variables)
        if closed.inconsistent or _same_reduced_variety(closed, context):
            continue
        saturated_generators = _saturate_generators(context.generators, splitter, context.variables)
        if saturated_generators is None:
            continue
        saturated = EqualityIdealContext.build(saturated_generators, context.variables)
        if saturated.inconsistent or _same_reduced_variety(saturated, context):
            continue
        if _same_reduced_variety(closed, saturated):
            continue
        return (
            _normalized_context_generators(closed),
            _normalized_context_generators(saturated),
        )
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


def equidimensional_decomposition(
    equations: Iterable[sp.Expr | sp.Equality],
    variables: Sequence[sp.Symbol],
    *,
    max_pieces: int = 128,
    _emit_certificate: bool = True,
) -> EquidimensionalDecomposition:
    """Decompose a polynomial zero set into exact dimension-labelled pieces.

    The algorithm has reduced-set semantics. Monomial Groebner bases are
    decomposed exactly through the minimal coordinate primes of their radical.
    More generally, exactly factorizable Groebner generators are decomposed by
    minimal hitting sets of their square-free factor supports. This extends the
    certified reduced-set decomposition to nonmonomial factor-incidence ideals
    and removes redundant/embedded factor branches without branch guessing.
    Radical-containment tests remove any remaining duplicate and contained
    branches. The returned union is exact.

    General primary decomposition is deliberately not guessed. When terminal
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

        real_reduction = _real_radical_reduction_generators(context)
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
        real_dimension = _real_dimension(context)
        if real_dimension < 0:
            continue
        if real_dimension == 0 and context.dimension > 0:
            point_contexts = _real_point_contexts(context)
            if point_contexts:
                source_dimension = (
                    context.dimension if algebraic_dim_override is None else algebraic_dim_override
                )
                real_contexts.extend((point, 0, source_dimension, True) for point in point_contexts)
                continue
        certified = _equidimensionality_certified(context) and real_dimension == context.dimension
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


def verify_decomposition_certificate(certificate: DecompositionCertificate) -> bool:
    """Replay and verify a decomposition certificate.

    Verification is intentionally strict: the normalized input is decomposed
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
    "DecompositionPieceCertificate",
    "DecompositionCertificate",
    "EquidimensionalDecomposition",
    "equidimensional_decomposition",
    "verify_decomposition_certificate",
]
