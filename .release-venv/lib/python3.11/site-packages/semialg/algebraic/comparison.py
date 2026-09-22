from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass

import sympy as sp

from .._zero_testing import certified_zero
from .cache import (
    CACHE,
    expr_key,
    mark_refinement_decisive,
    poly_key,
    sample_certificate_key,
    sample_identity_key,
)
from .exact_sign import exact_algebraic_sign
from .roots import _ordinary_root_count, refine_root_once
from .samples import AlgebraicRoot, RationalSample, Sample, sample_to_expr


def _sample_cache_identity(sample: Sample) -> object:
    # Comparison depends on the current isolating certificate, not only the
    # abstract algebraic root identity.
    return sample_certificate_key(sample)


def _structural_cache_pair(left: Sample, right: Sample) -> tuple[tuple[object, ...], int]:
    left_key = _sample_cache_identity(left)
    right_key = _sample_cache_identity(right)
    left_first = repr(left_key) <= repr(right_key)
    return ((left_key, right_key) if left_first else (right_key, left_key)), (
        1 if left_first else -1
    )


def _fiber_contexts_share_parent(left: AlgebraicRoot, right: AlgebraicRoot) -> bool:
    left_context = left.fiber_context
    right_context = right.fiber_context
    return (
        left_context is not None
        and right_context is not None
        and left_context.fiber_variable == right_context.fiber_variable
        and left_context.variables == right_context.variables
        and tuple(sample_identity_key(s) for s in left_context.parent_samples)
        == tuple(sample_identity_key(s) for s in right_context.parent_samples)
    )


@dataclass(frozen=True)
class FiberCoprimalityCandidate:
    """Compiled degree-zero subresultant for one symbolic fiber-family pair."""

    expression: sp.Expr
    parent_poly: sp.Poly


@dataclass(frozen=True)
class FiberSubresultantCandidate:
    """Compiled coefficient scan for one symbolic subresultant candidate."""

    expression: sp.Expr
    degree: int
    coefficient_expressions: tuple[sp.Expr, ...]


@dataclass
class FiberSubresultantSectionScan:
    """Persistent coefficient signs for one family pair at one parent section."""

    coefficient_signs: dict[tuple[int, int], int]
    candidate_index: int | None = None
    coefficient_index: int | None = None
    coefficient_sign: int = 0
    complete: bool = False


@dataclass
class FiberCommonFactorCertificate:
    """Persistent common-factor certificate for two families over one parent section."""

    gcd_poly: sp.Poly | None
    computed: bool = True
    root_memberships: dict[object, frozenset[int]] | None = None
    root_results: dict[object, dict[int, bool]] | None = None
    boundary_signs: dict[object, dict[sp.Rational, int]] | None = None

    def __post_init__(self) -> None:
        if self.root_memberships is None:
            self.root_memberships = {}
        if self.root_results is None:
            self.root_results = {}
        if self.boundary_signs is None:
            self.boundary_signs = {}


@dataclass
class SharedFiberPartition:
    """One ranked rational partition shared by every fiber family at a parent cell."""

    contexts: dict[object, object]
    roots: dict[object, tuple[AlgebraicRoot, ...]]
    boundary_ranks: dict[sp.Rational, dict[object, int]]


def _shared_fiber_partition_key(root: AlgebraicRoot) -> object | None:
    context = root.fiber_context
    if context is None:
        return None
    return (
        "shared-fiber-partition",
        tuple(sample_identity_key(sample) for sample in context.parent_samples),
        context.variables,
        context.fiber_variable,
    )


def _shared_fiber_partition(root: AlgebraicRoot) -> SharedFiberPartition | None:
    key = _shared_fiber_partition_key(root)
    context = root.fiber_context
    if key is None or context is None:
        return None
    partition = CACHE.shared_fiber_partitions.get(key)
    if not isinstance(partition, SharedFiberPartition):
        partition = SharedFiberPartition({}, {}, {})
        CACHE.shared_fiber_partitions.put(key, partition)
    family = context._sturm_identity
    if family not in partition.contexts:
        partition.contexts[family] = context

    return partition


def _rank_from_known_root_intervals(
    roots: tuple[AlgebraicRoot, ...], value: sp.Rational
) -> int | None:
    """Derive a family rank from isolated root nodes without a Sturm endpoint."""
    below = 0
    for root in roots:
        if root.interval.right <= value:
            below += 1
        elif value < root.interval.left:
            break
        else:
            return None
    return below


def _classify_shared_boundary(
    partition: SharedFiberPartition,
    value: sp.Rational,
    *,
    required: tuple[object, ...] = (),
) -> dict[object, int]:
    """Insert one global boundary and classify all families that are free to classify.

    Root-node topology classifies every participating family whenever possible.
    A new Sturm endpoint is paid only for families whose rank is required by the
    current comparison; unrelated overlapping root nodes remain intentionally
    unresolved until a later query actually needs them.
    """
    value = sp.Rational(value)
    ranks = partition.boundary_ranks.get(value)
    if ranks is None:
        ranks = {}
        partition.boundary_ranks[value] = ranks
    required_set = set(required)
    for family, context in tuple(partition.contexts.items()):
        if family in ranks:
            continue
        rank = _rank_from_known_root_intervals(partition.roots.get(family, ()), value)
        if rank is None and family in required_set:
            rank = context.roots_left_of(value)
        if rank is not None:
            context._remember_boundary(value, rank)
            ranks[family] = rank
    return ranks


def _root_below_shared_boundary(
    partition: SharedFiberPartition, root: AlgebraicRoot, value: sp.Rational
) -> bool:
    """Place only the selected ranked root relative to a shared boundary."""
    context = root.fiber_context
    if context is None:
        raise ArithmeticError("fiber root has no context")
    family = context._sturm_identity
    ranks = partition.boundary_ranks.setdefault(sp.Rational(value), {})
    rank = ranks.get(family)
    if rank is not None:
        return root.root_index < rank
    cert = context.root_rank_cert(root.root_index)
    if cert is not None:
        if value < cert.interval.left:
            return False
        if value >= cert.interval.right:
            return True
        if (
            cert.node_root_count == 1
            and root.multiplicity % 2 == 1
            and context.polynomial.degree(context.fiber_variable) >= 3
        ):
            try:
                value_sign = context.sign_at_rational(value)
                if value_sign == 0:
                    return True
                left_sign = context.sign_at_rational(cert.interval.left)
                if left_sign == 0:
                    return True
                return value_sign != left_sign
            except (ArithmeticError, ValueError, TypeError):
                pass
    below = context.root_on_or_below(root.root_index, value)
    boundary = context._boundary_cert(value)
    if boundary is not None:
        ranks[family] = boundary.root_rank
    return below


def prepare_shared_fiber_partition(roots: Sequence[AlgebraicRoot]) -> None:
    """Register all root families in one K-way merge before comparisons begin."""
    by_family: dict[object, list[AlgebraicRoot]] = {}
    partition = None
    for root in roots:
        current = _shared_fiber_partition(root)
        if current is None:
            continue
        partition = current
        context = root.fiber_context
        if context is not None:
            by_family.setdefault(context._sturm_identity, []).append(root)
    if partition is None:
        return
    for family, family_roots in by_family.items():
        partition.roots[family] = tuple(sorted(family_roots, key=lambda root: root.root_index))


def _order_from_shared_boundaries(
    partition: SharedFiberPartition, left: AlgebraicRoot, right: AlgebraicRoot
) -> tuple[int, sp.Rational] | None:
    left_context = left.fiber_context
    right_context = right.fiber_context
    if left_context is None or right_context is None:
        return None
    left_family = left_context._sturm_identity
    right_family = right_context._sturm_identity
    for boundary in sorted(partition.boundary_ranks):
        ranks = partition.boundary_ranks[boundary]
        if left_family not in ranks or right_family not in ranks:
            continue
        left_below = left.root_index < ranks[left_family]
        right_below = right.root_index < ranks[right_family]
        if left_below != right_below:
            return (-1 if left_below else 1), boundary
    return None


def _one_sided_shared_boundary(
    partition: SharedFiberPartition, left: AlgebraicRoot, right: AlgebraicRoot
) -> tuple[int, sp.Rational] | None:
    """Reuse a boundary already ranked for one family with at most one new decision.

    A successful rank ``r`` classifies every sibling root of that family at the
    boundary.  Later family pairs can therefore probe only the still-unknown
    selected root instead of paying two endpoint decisions at a fresh midpoint.
    """
    left_context = left.fiber_context
    right_context = right.fiber_context
    if left_context is None or right_context is None:
        return None
    left_family = left_context._sturm_identity
    right_family = right_context._sturm_identity
    overlap_left = max(left.interval.left, right.interval.left)
    overlap_right = min(left.interval.right, right.interval.right)
    if overlap_left >= overlap_right:
        return None
    candidates = []
    center = sp.Rational(overlap_left + overlap_right, 2)
    for boundary, ranks in partition.boundary_ranks.items():
        if not (overlap_left < boundary < overlap_right):
            continue
        left_known = left_family in ranks
        right_known = right_family in ranks
        if left_known == right_known:
            continue
        candidates.append((abs(boundary - center), boundary, left_known))
    for _, boundary, left_known in sorted(candidates):
        ranks = partition.boundary_ranks[boundary]
        if left_known:
            left_below = left.root_index < ranks[left_family]
            right_below = _root_below_shared_boundary(partition, right, boundary)
        else:
            left_below = _root_below_shared_boundary(partition, left, boundary)
            right_below = right.root_index < ranks[right_family]
        if left_below != right_below:
            return (-1 if left_below else 1), boundary
    return None


@dataclass
class FiberFamilyRelationship:
    """Cached algebraic relationship between two fiber-polynomial families."""

    resultant_zero: bool | None
    common_root_pairs: set[tuple[int, int]]
    root_order_pairs: dict[tuple[int, int], int]
    common_gcd: sp.Poly | None = None
    gcd_computed: bool = False


def _fiber_family_pair_key(left: AlgebraicRoot, right: AlgebraicRoot) -> tuple[object, int] | None:
    if not _fiber_contexts_share_parent(left, right):
        return None
    left_context = left.fiber_context
    right_context = right.fiber_context
    if left_context is None or right_context is None:
        return None
    parent_key = tuple(sample_identity_key(s) for s in left_context.parent_samples)
    left_key = poly_key(left_context.polynomial)
    right_key = poly_key(right_context.polynomial)
    left_first = repr(left_key) <= repr(right_key)
    pair = (left_key, right_key) if left_first else (right_key, left_key)
    return ("fiber-family", parent_key, left_context.fiber_variable, pair), (
        1 if left_first else -1
    )


def _fiber_family_relationship(
    left: AlgebraicRoot, right: AlgebraicRoot
) -> tuple[FiberFamilyRelationship, int] | None:
    """Return one cached resultant certificate for a stable pair of fiber families."""
    keyed = _fiber_family_pair_key(left, right)
    if keyed is None:
        return None
    key, orientation = keyed
    cached = CACHE.family_relationships.get(key)
    if cached is not None:
        CACHE.stats.family_relationship_hits += 1
        return cached, orientation
    CACHE.stats.family_relationship_misses += 1
    left_context = left.fiber_context
    right_context = right.fiber_context
    if left_context is None or right_context is None:
        raise RuntimeError("fiber-family cache key requires two fiber contexts")
    from .signs import sign_at_sample

    fiber = left_context.fiber_variable
    parents = left_context.variables[:-1]
    resultant_zero = False
    try:
        domain = sp.QQ.frac_field(*parents) if parents else sp.QQ
        left_compile_key = (
            "fiber-family-field-poly",
            poly_key(left_context.polynomial),
            fiber,
            parents,
        )
        right_compile_key = (
            "fiber-family-field-poly",
            poly_key(right_context.polynomial),
            fiber,
            parents,
        )
        left_poly = CACHE.compiled_polys.get(left_compile_key)
        if not isinstance(left_poly, sp.Poly):
            left_poly = sp.Poly(left_context.polynomial.as_expr(), fiber, domain=domain)
            CACHE.compiled_polys.put(left_compile_key, left_poly)
        right_poly = CACHE.compiled_polys.get(right_compile_key)
        if not isinstance(right_poly, sp.Poly):
            right_poly = sp.Poly(right_context.polynomial.as_expr(), fiber, domain=domain)
            CACHE.compiled_polys.put(right_compile_key, right_poly)
        # The resultant is computed once for the family pair.  Normalize its
        # rational-function representation once before signing it at the parent.
        resultant = sp.cancel(sp.resultant(left_poly.as_expr(), right_poly.as_expr(), fiber))
        numerator, denominator = sp.fraction(resultant)

        # The ordering relationship has already paid for the same degree-zero
        # subresultant used by fiber-zero certification.  When its field-domain
        # normalization has only a rational unit denominator, publish the exact
        # parent polynomial into the first-class coprimality cache so the zero
        # path never recomputes this symbolic resultant.
        if denominator.is_Rational:
            left_family_key = _common_factor_family_key(left_context.polynomial)
            right_family_key = _common_factor_family_key(right_context.polynomial)
            family_pair = tuple(sorted((left_family_key, right_family_key), key=repr))
            candidate_key = ("fiber-degree-zero-subresultant", fiber, family_pair)
            if CACHE.compiled_polys.get(candidate_key) is None:
                candidate_expr = sp.expand(numerator / denominator)
                parent_poly = (
                    sp.Poly(candidate_expr, *parents, domain=sp.QQ)
                    if parents
                    else sp.Poly(candidate_expr)
                )
                CACHE.compiled_polys.put(
                    candidate_key, FiberCoprimalityCandidate(candidate_expr, parent_poly)
                )

        numerator_sign = sign_at_sample(
            sp.Poly(sp.expand(numerator), *parents, domain=sp.QQ), left_context.parent_samples
        )
        denominator_sign = sign_at_sample(
            sp.Poly(sp.expand(denominator), *parents, domain=sp.QQ), left_context.parent_samples
        )
        if denominator_sign == 0:
            # The generic resultant lives over QQ(parents).  At an exceptional
            # parent section its normalized denominator can vanish even though
            # the *specialized* fiber polynomials are perfectly regular (degree
            # drops are the usual cause).  Do not evaluate through that pole.
            # Defer the zero/nonzero decision to the specialization-safe
            # subresultant path below.
            resultant_zero = None
        else:
            resultant_zero = numerator_sign == 0
    except (
        sp.PolynomialError,
        sp.polys.polyerrors.CoercionFailed,
        sp.polys.polyerrors.DomainError,
    ):
        resultant_zero = None
    relationship = FiberFamilyRelationship(resultant_zero, set(), {})
    CACHE.family_relationships.put(key, relationship)
    if resultant_zero is None:
        # A pole in the generic coefficient-field resultant is not a
        # mathematical indeterminacy.  Specialize the subresultant sequence at
        # the parent sample and use the first nonzero positive-degree member as
        # the exact common-factor certificate.  This keeps exceptional sections
        # inside the parent tower instead of constructing global algebraic
        # substitutions.
        common = _fiber_common_gcd_context(left, right, relationship)
        relationship.resultant_zero = common is not None
    return relationship, orientation


def _fiber_root_order_from_family_cert(left: AlgebraicRoot, right: AlgebraicRoot) -> int | None:
    """Order fiber roots through one multi-family ranked rational partition.

    Every separator inserted for any family pair is evaluated against all fiber
    families already seen over the same parent cell.  Later pair comparisons can
    therefore consume global partition topology instead of starting another
    family-pair binary Sturm search.
    """
    related = _fiber_family_relationship(left, right)
    if related is None:
        return None
    relationship, orientation = related
    if relationship.resultant_zero is None:
        return None
    pair = (
        (left.root_index, right.root_index)
        if orientation == 1
        else (right.root_index, left.root_index)
    )
    if pair in relationship.common_root_pairs:
        return 0
    cached = relationship.root_order_pairs.get(pair)
    if cached is not None:
        return orientation * cached

    left_context = left.fiber_context
    right_context = right.fiber_context
    if left_context is None or right_context is None:
        return None
    partition = _shared_fiber_partition(left)
    if partition is None:
        return None
    _shared_fiber_partition(right)  # register the second family
    for root in (left, right):
        context = root.fiber_context
        if context is not None and context._sturm_identity not in partition.roots:
            partition.roots[context._sturm_identity] = (root,)

    known = _order_from_shared_boundaries(partition, left, right)
    if known is None:
        known = _one_sided_shared_boundary(partition, left, right)
    if known is not None:
        result, separator = known
        relationship.root_order_pairs[pair] = orientation * result
        lower, upper = (left, right) if result < 0 else (right, left)
        CACHE.sector_separators.put(
            ("sector-between", sample_identity_key(lower), sample_identity_key(upper)), separator
        )
        return result

    search_left = min(left.interval.left, right.interval.left)
    search_right = max(left.interval.right, right.interval.right)
    for _ in range(128):
        if search_left >= search_right:
            return None
        separator = sp.Rational(search_left + search_right, 2)
        # Classify all families that topology can answer for free, then ask
        # only about the two selected roots.  A sibling interval crossing this
        # separator must not force a whole-family Sturm rank when the selected
        # root is already known to one side.
        _classify_shared_boundary(partition, separator)
        left_below = _root_below_shared_boundary(partition, left, separator)
        right_below = _root_below_shared_boundary(partition, right, separator)
        if left_below != right_below:
            result = -1 if left_below else 1
            relationship.root_order_pairs[pair] = orientation * result
            lower, upper = (left, right) if result < 0 else (right, left)
            CACHE.sector_separators.put(
                ("sector-between", sample_identity_key(lower), sample_identity_key(upper)),
                separator,
            )
            return result
        if left_below:
            search_right = separator
        else:
            search_left = separator
    return None


def _common_factor_family_key(poly: sp.Poly) -> object:
    """Canonical zero-set key for a QQ fiber family, ignoring rational units."""
    if poly.domain == sp.QQ and not poly.is_zero:
        terms = tuple(poly.terms())
        leading = terms[0][1]
        return (
            tuple(poly.gens),
            poly.domain,
            tuple((powers, coeff / leading) for powers, coeff in terms),
        )
    return poly_key(poly)


def _fiber_common_factor_cache_key(left_context, right_context) -> tuple[object, ...]:
    family_pair = tuple(
        sorted(
            (
                _common_factor_family_key(left_context.polynomial),
                _common_factor_family_key(right_context.polynomial),
            ),
            key=repr,
        )
    )
    return (
        "fiber-common-gcd",
        tuple(sample_identity_key(sample) for sample in left_context.parent_samples),
        left_context.fiber_variable,
        family_pair,
    )


def _fiber_degree_zero_coprimality_candidate(
    left_context, right_context, family_pair: tuple[object, ...]
) -> FiberCoprimalityCandidate | None:
    """Return the compiled degree-zero subresultant for a family pair.

    The resultant of the two QQ fiber families is the degree-zero
    subresultant.  Its nonvanishing after specialization is sufficient to
    prove that the specialized families are coprime.  Vanishing is only
    inconclusive: a parent-section degree drop can make the generic
    resultant vanish even when the specialized families are coprime, so the
    higher subresultant scan must then take over.
    """
    fiber = left_context.fiber_variable
    parents = left_context.variables[:-1]
    key = ("fiber-degree-zero-subresultant", fiber, family_pair)
    cached = CACHE.compiled_polys.get(key)
    if cached is False:
        return None
    if isinstance(cached, FiberCoprimalityCandidate):
        return cached
    try:
        resultant = sp.expand(
            sp.resultant(
                left_context.polynomial.as_expr(),
                right_context.polynomial.as_expr(),
                fiber,
            )
        )
        parent_poly = sp.Poly(resultant, *parents, domain=sp.QQ) if parents else sp.Poly(resultant)
    except (
        sp.PolynomialError,
        sp.polys.polyerrors.CoercionFailed,
        sp.polys.polyerrors.DomainError,
    ):
        CACHE.compiled_polys.put(key, False)
        return None
    candidate = FiberCoprimalityCandidate(resultant, parent_poly)
    CACHE.compiled_polys.put(key, candidate)
    return candidate


def _fiber_degree_zero_coprime_at_section(
    left_context, right_context, family_pair: tuple[object, ...], gcd_cache_key: object
) -> bool:
    """Certify coprimality from only the degree-zero subresultant."""
    section_key = ("fiber-degree-zero-subresultant-sign", gcd_cache_key)
    cached = CACHE.fiber_subresultant_section_scans.get(section_key)
    if isinstance(cached, int):
        return cached != 0
    candidate = _fiber_degree_zero_coprimality_candidate(left_context, right_context, family_pair)
    if candidate is None:
        return False
    if candidate.parent_poly.total_degree() == 0:
        sign = int(sp.sign(sp.Rational(candidate.parent_poly.as_expr())))
    else:
        from .signs import sign_at_sample

        sign = int(sign_at_sample(candidate.parent_poly, left_context.parent_samples))
    CACHE.fiber_subresultant_section_scans.put(section_key, sign)
    return sign != 0


def _fiber_common_gcd_context(
    left: AlgebraicRoot,
    right: AlgebraicRoot,
    relationship: FiberFamilyRelationship,
    *,
    instantiate_context: bool = True,
):
    """Return the specialized common-factor family from a subresultant certificate."""
    left_context = left.fiber_context
    right_context = right.fiber_context
    if left_context is None or right_context is None:
        return None

    def materialize(gcd_poly: sp.Poly):
        if not instantiate_context:
            return gcd_poly
        from .samples import FiberRootContext

        return FiberRootContext(
            gcd_poly,
            left_context.variables,
            left_context.parent_samples,
            left_context.fiber_variable,
        )

    if relationship.gcd_computed:
        if relationship.common_gcd is None:
            return None
        return materialize(relationship.common_gcd)
    relationship.gcd_computed = True

    left_family_key = _common_factor_family_key(left_context.polynomial)
    right_family_key = _common_factor_family_key(right_context.polynomial)
    family_pair = tuple(sorted((left_family_key, right_family_key), key=repr))
    gcd_cache_key = _fiber_common_factor_cache_key(left_context, right_context)
    cached_certificate = CACHE.fiber_common_gcds.get(gcd_cache_key)
    if cached_certificate is False:
        # A nonzero degree-zero specialized subresultant certified this
        # section-specific family pair as coprime.
        relationship.common_gcd = None
        return None
    if isinstance(cached_certificate, FiberCommonFactorCertificate):
        relationship.common_gcd = cached_certificate.gcd_poly
        if cached_certificate.gcd_poly is None:
            return None
        return materialize(cached_certificate.gcd_poly)
    fiber = left_context.fiber_variable
    parents = left_context.variables[:-1]

    # Most zero queries are negative.  Test the degree-zero subresultant
    # first, without constructing the full PRS.  A nonzero value at this
    # exact parent section is a complete coprimality certificate.  Only a
    # zero value pays for the higher-degree subresultants below.
    try:
        if _fiber_degree_zero_coprime_at_section(
            left_context, right_context, family_pair, gcd_cache_key
        ):
            relationship.common_gcd = None
            CACHE.fiber_common_gcds.put(gcd_cache_key, False)
            return None
    except (ArithmeticError, ValueError, TypeError):
        # An inconclusive section sign must not become a coprimality proof.
        # The specialization-safe PRS below remains authoritative.
        pass

    from .signs import sign_at_sample

    try:
        subresultant_key = ("fiber-subresultants", fiber, family_pair)
        scan_plan = CACHE.fiber_subresultants.get(subresultant_key)
        if not (
            isinstance(scan_plan, tuple)
            and all(isinstance(item, FiberSubresultantCandidate) for item in scan_plan)
        ):
            candidates: list[FiberSubresultantCandidate] = []
            for candidate in reversed(
                sp.subresultants(
                    left_context.polynomial.as_expr(),
                    right_context.polynomial.as_expr(),
                    fiber,
                )
            ):
                candidate_expr = sp.expand(candidate)
                candidate_key = (
                    "fiber-subresultant-univariate",
                    expr_key(candidate_expr),
                    fiber,
                )
                candidate_poly = CACHE.compiled_polys.get(candidate_key)
                if not isinstance(candidate_poly, sp.Poly):
                    candidate_poly = sp.Poly(candidate_expr, fiber)
                    CACHE.compiled_polys.put(candidate_key, candidate_poly)
                # Keep degree-zero subresultants in the section scan.  A
                # nonzero specialized constant is the strongest possible
                # certificate here: the two specialized fiber families are
                # coprime, so no common-factor context (or its Sturm
                # machinery) needs to be constructed.  Constants that vanish
                # at the selected parent section are skipped naturally by the
                # coefficient-sign scan, allowing the next positive-degree
                # subresultant to describe the section-specific gcd after a
                # degree drop.
                # Keep coefficient expressions lazy.  Most CMXY:14 family
                # pairs are coprime, so the reversed PRS stops after signing
                # the degree-zero subresultant.  Compiling every coefficient
                # Poly for all earlier PRS members here wasted work on the
                # overwhelmingly common negative zero query.
                coefficient_expressions = tuple(
                    sp.expand(coefficient) for coefficient in candidate_poly.all_coeffs()
                )
                candidates.append(
                    FiberSubresultantCandidate(
                        candidate_expr, candidate_poly.degree(), coefficient_expressions
                    )
                )
            scan_plan = tuple(candidates)
            CACHE.fiber_subresultants.put(subresultant_key, scan_plan)

        section_scan_key = ("fiber-subresultant-section-scan", gcd_cache_key)
        section_scan = CACHE.fiber_subresultant_section_scans.get(section_scan_key)
        if not isinstance(section_scan, FiberSubresultantSectionScan):
            section_scan = FiberSubresultantSectionScan({})
            # Publish the mutable scan before doing section sign work.  If a
            # later coefficient is temporarily inconclusive, a retry resumes
            # after every coefficient whose sign was already certified.
            CACHE.fiber_subresultant_section_scans.put(section_scan_key, section_scan)
        if section_scan.candidate_index is None and not section_scan.complete:
            for candidate_index, candidate in enumerate(scan_plan):
                for coefficient_index, coefficient_expr in enumerate(
                    candidate.coefficient_expressions
                ):
                    coefficient_key = (candidate_index, coefficient_index)
                    coefficient_sign = section_scan.coefficient_signs.get(coefficient_key)
                    if coefficient_sign is None:
                        compiled_key = (
                            "fiber-subresultant-coefficient",
                            expr_key(coefficient_expr),
                            parents,
                        )
                        coefficient_poly = CACHE.compiled_polys.get(compiled_key)
                        if not isinstance(coefficient_poly, sp.Poly):
                            coefficient_poly = sp.Poly(coefficient_expr, *parents, domain=sp.QQ)
                            CACHE.compiled_polys.put(compiled_key, coefficient_poly)
                        if coefficient_poly.total_degree() == 0:
                            constant = sp.Rational(coefficient_poly.as_expr())
                            coefficient_sign = int(sp.sign(constant))
                        else:
                            coefficient_sign = int(
                                sign_at_sample(coefficient_poly, left_context.parent_samples)
                            )
                        section_scan.coefficient_signs[coefficient_key] = coefficient_sign
                    if coefficient_sign != 0:
                        section_scan.candidate_index = candidate_index
                        section_scan.coefficient_index = coefficient_index
                        section_scan.coefficient_sign = coefficient_sign
                        break
                if section_scan.candidate_index is not None:
                    break
            section_scan.complete = True

        if section_scan.candidate_index is not None:
            candidate = scan_plan[section_scan.candidate_index]
            # The first nonzero member of the reversed subresultant PRS is a
            # unit exactly when the specialized families are coprime.  Publish
            # that negative certificate immediately; materializing a
            # FiberRootContext for a unit gcd would only turn a cheap
            # nonvanishing proof into a second root-isolation problem.
            if candidate.degree == 0:
                relationship.common_gcd = None
                CACHE.fiber_common_gcds.put(gcd_cache_key, False)
                return None
            gcd_key = (
                "fiber-subresultant-gcd-poly",
                expr_key(candidate.expression),
                tuple(left_context.variables),
            )
            gcd_poly = CACHE.compiled_polys.get(gcd_key)
            if not isinstance(gcd_poly, sp.Poly):
                gcd_poly = sp.Poly(candidate.expression, *left_context.variables, domain=sp.QQ)
                CACHE.compiled_polys.put(gcd_key, gcd_poly)
            relationship.common_gcd = gcd_poly
            CACHE.fiber_common_gcds.put(gcd_cache_key, FiberCommonFactorCertificate(gcd_poly))
            return materialize(gcd_poly)
    except (
        ArithmeticError,
        ValueError,
        TypeError,
        sp.PolynomialError,
        sp.polys.polyerrors.DomainError,
    ):
        return None
    CACHE.fiber_common_gcds.put(gcd_cache_key, FiberCommonFactorCertificate(None))
    return None


def _linear_common_factor_zero_at_root(
    context, certificate: FiberCommonFactorCertificate, root: AlgebraicRoot
) -> bool | None:
    """Answer linear common-factor membership without a Sturm context.

    A linear specialized subresultant has at most one root.  Test only the two
    boundaries of the requested ranked section.  A positive answer identifies
    the unique member and therefore certifies every sibling rank at once;
    negative answers are retained incrementally without eagerly signing unused
    partition boundaries.
    """
    gcd_poly = certificate.gcd_poly
    if gcd_poly is None or gcd_poly.degree(context.fiber_variable) != 1:
        return None
    family = context._sturm_identity
    complete = certificate.root_memberships.get(family)
    if complete is not None:
        return root.root_index in complete
    partial = certificate.root_results.setdefault(family, {})
    if root.root_index in partial:
        return partial[root.root_index]

    from .signs import sign_at_sample

    fiber = context.fiber_variable
    parents = context.variables[:-1]

    def endpoint_sign(boundary: sp.Rational) -> int:
        specialized = gcd_poly.eval(fiber, sp.Rational(boundary))
        if isinstance(specialized, sp.Poly):
            parent_poly = specialized
            if tuple(parent_poly.gens) != parents or parent_poly.domain != sp.QQ:
                parent_poly = sp.Poly(parent_poly.as_expr(), *parents, domain=sp.QQ)
        else:
            parent_poly = sp.Poly(sp.expand(specialized), *parents, domain=sp.QQ)
        return int(sign_at_sample(parent_poly, context.parent_samples))

    left = sp.Rational(root.interval.left)
    right = sp.Rational(root.interval.right)
    left_sign = endpoint_sign(left)
    right_sign = endpoint_sign(right)
    member_index: int | None = None
    if left_sign * right_sign < 0:
        member_index = root.root_index
    else:
        for boundary, sign in ((left, left_sign), (right, right_sign)):
            if sign != 0 or context.squarefree_sign_at_rational(boundary) != 0:
                continue
            rank = context.roots_left_of(boundary)
            if rank == root.root_index:
                member_index = rank
                break
            # The unique linear root is a different ranked boundary.  Its rank
            # completely classifies this defining family as well.
            if 0 <= rank < len(context.isolate_all_roots()):
                certificate.root_memberships[family] = frozenset({rank})
                partial.update(
                    {index: index == rank for index in range(len(context.isolate_all_roots()))}
                )
                return False
    if member_index is not None:
        count = len(context.isolate_all_roots())
        members = frozenset({member_index})
        certificate.root_memberships[family] = members
        partial.update({index: index == member_index for index in range(count)})
        return True
    partial[root.root_index] = False
    return False


def _higher_degree_common_factor_positive_at_root(
    context, certificate: FiberCommonFactorCertificate, root: AlgebraicRoot
) -> bool | None:
    """Return True from a direct degree-2/3 sign-change certificate.

    This is deliberately a positive-only certificate.  A sign change across a
    defining-family isolating interval proves that the specialized common
    factor has a root there, hence (because it divides the defining family)
    that the selected ranked root belongs to it.  No sign change is not a
    negative certificate when multiplicities may specialize, so callers fall
    back to the existing specialization-safe Sturm context.
    """
    gcd_poly = certificate.gcd_poly
    if gcd_poly is None:
        return None
    fiber = context.fiber_variable
    if gcd_poly.degree(fiber) not in (2, 3):
        return None
    family = context._sturm_identity
    signs = certificate.boundary_signs.setdefault(family, {})
    from .signs import sign_at_sample

    def endpoint_sign(boundary: sp.Rational) -> int:
        q = sp.Rational(boundary)
        cached = signs.get(q)
        if cached is not None:
            return cached
        specialized = gcd_poly.eval(fiber, q)
        parents = context.variables[:-1]
        if isinstance(specialized, sp.Poly):
            parent_poly = specialized
            if tuple(parent_poly.gens) != parents or parent_poly.domain != sp.QQ:
                parent_poly = sp.Poly(parent_poly.as_expr(), *parents, domain=sp.QQ)
        else:
            parent_poly = sp.Poly(sp.expand(specialized), *parents, domain=sp.QQ)
        result = int(sign_at_sample(parent_poly, context.parent_samples))
        signs[q] = result
        return result

    left_sign = endpoint_sign(root.interval.left)
    right_sign = endpoint_sign(root.interval.right)
    return True if left_sign * right_sign < 0 else None


def _common_factor_zero_at_rank(
    context, gcd_context, certificate: FiberCommonFactorCertificate, root: AlgebraicRoot
) -> bool:
    """Consume or extend the section-wide ranked common-factor certificate.

    The first two distinct roots are answered individually.  A third query is
    evidence that this family pair is being scanned by a sign table, so the
    complete defining-family partition is classified at once.  All sibling
    answers then become O(1), while one-off pairs avoid eager isolation work.
    """
    family = context._sturm_identity
    direct = _linear_common_factor_zero_at_root(context, certificate, root)
    if direct is not None:
        return direct
    complete = certificate.root_memberships.get(family)
    if complete is not None:
        return root.root_index in complete
    partial = certificate.root_results.setdefault(family, {})
    if root.root_index in partial:
        return partial[root.root_index]
    degree = certificate.gcd_poly.degree(context.fiber_variable) if certificate.gcd_poly else -1
    # One-off and two-off family pairs keep the already-cheap ranked-boundary
    # root_count path.  Only repeated degree-2/3 scans pay for direct endpoint
    # signs, where those signs are shared by adjacent ranked intervals.
    if len(partial) >= 2 and degree in (2, 3):
        direct = _higher_degree_common_factor_positive_at_root(context, certificate, root)
        if direct:
            partial[root.root_index] = True
            members = frozenset(index for index, value in partial.items() if value)
            if len(members) == degree:
                count = len(context.isolate_all_roots())
                certificate.root_memberships[family] = members
                partial.update({index: index in members for index in range(count)})
            return True
    if len(partial) < 5:
        result = gcd_context.root_count(root.interval.left, root.interval.right) > 0
        partial[root.root_index] = result
        if result and degree in (2, 3):
            members = frozenset(index for index, value in partial.items() if value)
            if len(members) == degree:
                count = len(context.isolate_all_roots())
                certificate.root_memberships[family] = members
                partial.update({index: index in members for index in range(count)})
        return result

    intervals = context.isolate_all_roots()
    boundaries = sorted(
        {sp.Rational(x) for interval in intervals for x in (interval.left, interval.right)}
    )
    ranks = {boundary: gcd_context.roots_left_of(boundary) for boundary in boundaries}
    members = frozenset(
        index
        for index, interval in enumerate(intervals)
        if ranks[sp.Rational(interval.right)] > ranks[sp.Rational(interval.left)]
    )
    certificate.root_memberships[family] = members
    partial.update({index: index in members for index in range(len(intervals))})
    return root.root_index in members


def _fiber_sign_cert_key(target_poly: sp.Poly, root: AlgebraicRoot) -> tuple[object, ...]:
    """Stable key for a polynomial sign on a ranked fiber section."""
    context = root.fiber_context
    partition_identity: object = sample_identity_key(root)
    if context is not None:
        cert = context.root_rank_cert(root.root_index)
        if cert is not None:
            partition_identity = (
                "ranked-section",
                context._sturm_identity,
                cert.root_index,
            )
    return ("fiber-sign", poly_key(target_poly), partition_identity)


def _fiber_cached_sign(target_poly: sp.Poly, root: AlgebraicRoot) -> int | None:
    cached = CACHE.fiber_sign_certs.get(_fiber_sign_cert_key(target_poly, root))
    return int(cached) if cached in (-1, 0, 1) else None


def _publish_fiber_sign(target_poly: sp.Poly, root: AlgebraicRoot, sign: int) -> None:
    """Publish a certified section sign for all consumers of the ranked partition."""
    sign = int(sign)
    if sign not in (-1, 0, 1):
        raise ValueError("fiber sign certificate must be -1, 0, or 1")
    CACHE.fiber_sign_certs.put(_fiber_sign_cert_key(target_poly, root), sign)


def _fiber_poly_zero_at_root(target_poly: sp.Poly, root: AlgebraicRoot) -> bool:
    """Certify a fiber-family zero once and share its common-factor certificate."""
    context = root.fiber_context
    if context is None or tuple(target_poly.gens) != context.variables:
        return False
    target_key = poly_key(target_poly)
    cert_key = (
        "fiber-poly-zero",
        target_key,
        sample_identity_key(root),
    )
    # A fiber section is, by construction, a zero of its own defining family.
    # Recognize this before constructing a subresultant certificate; exceptional
    # degree-drop sections otherwise force unnecessary coefficient-field work.
    defining_poly = context.polynomial
    if tuple(defining_poly.gens) != context.variables or defining_poly.domain != sp.QQ:
        defining_poly = sp.Poly(defining_poly.as_expr(), *context.variables, domain=sp.QQ)
    defining_key = poly_key(defining_poly)
    if target_key == defining_key:
        CACHE.fiber_zero_certs.put(cert_key, True)
        _publish_fiber_sign(target_poly, root, 0)
        return True
    # Polynomial ideal membership is stronger than a section-specific GCD
    # query and is independent of the selected parent section.  Consume it
    # before constructing any subresultant/common-factor certificate.
    try:
        if target_poly.rem(defining_poly).is_zero:
            CACHE.fiber_zero_certs.put(cert_key, True)
            _publish_fiber_sign(target_poly, root, 0)
            return True
    except (sp.PolynomialError, sp.polys.polyerrors.DomainError):
        pass
    sign_cert = _fiber_cached_sign(target_poly, root)
    if sign_cert is not None:
        return sign_cert == 0
    cached = CACHE.fiber_zero_certs.get(cert_key)
    if cached is not None:
        return bool(cached)

    from .samples import FiberRootContext

    target_context = FiberRootContext(
        target_poly, context.variables, context.parent_samples, context.fiber_variable
    )
    target_root = AlgebraicRoot(
        root.polynomial,
        root.interval,
        root.root_index,
        1,
        None,
        target_context,
    )
    # This certificate asks whether the specialized target family vanishes at
    # this particular root.  Do not infer a negative result from the generic
    # family resultant: specialization can create a common section.
    relationship = FiberFamilyRelationship(True, set(), {})
    gcd_poly = _fiber_common_gcd_context(root, target_root, relationship, instantiate_context=False)
    gcd_cache_key = _fiber_common_factor_cache_key(context, target_context)
    if not isinstance(gcd_poly, sp.Poly):
        # A nonzero constant subresultant is an exact section-specific
        # coprimality certificate.  Cache this root query as negative without
        # constructing a common-factor certificate or FiberRootContext.
        if CACHE.fiber_common_gcds.get(gcd_cache_key) is False:
            CACHE.fiber_zero_certs.put(cert_key, False)
        # Other failures remain unknown and must not become negative proofs.
        return False
    common_factor = CACHE.fiber_common_gcds.get(gcd_cache_key)
    if not isinstance(common_factor, FiberCommonFactorCertificate):
        return False
    try:
        direct = _linear_common_factor_zero_at_root(context, common_factor, root)
        if direct is not None:
            result = direct
        else:
            gcd_context = FiberRootContext(
                gcd_poly, context.variables, context.parent_samples, context.fiber_variable
            )
            result = _common_factor_zero_at_rank(context, gcd_context, common_factor, root)
    except (ArithmeticError, ValueError):
        return False
    CACHE.fiber_zero_certs.put(cert_key, result)
    if result:
        _publish_fiber_sign(target_poly, root, 0)
    return result


def _fiber_common_root_is_cached_or_certified(left: AlgebraicRoot, right: AlgebraicRoot) -> bool:
    related = _fiber_family_relationship(left, right)
    if related is None:
        return False
    relationship, orientation = related
    if not relationship.resultant_zero:
        return False
    pair = (
        (left.root_index, right.root_index)
        if orientation == 1
        else (right.root_index, left.root_index)
    )
    if pair in relationship.common_root_pairs:
        return True
    overlap_left = max(left.interval.left, right.interval.left)
    overlap_right = min(left.interval.right, right.interval.right)
    if overlap_left >= overlap_right:
        return False
    left_context = left.fiber_context
    right_context = right.fiber_context
    if left_context is None or right_context is None:
        return False
    gcd_context = _fiber_common_gcd_context(left, right, relationship)
    if gcd_context is not None:
        try:
            if gcd_context.root_count(overlap_left, overlap_right) > 0:
                relationship.common_root_pairs.add(pair)
                return True
        except (ArithmeticError, ValueError):
            pass
    try:
        # A vanishing family resultant only says that *some* sections meet.
        # Certify this particular root pair by evaluating one unspecialized
        # family polynomial at the other's complete tower sample.  This is
        # tower-native and avoids the unsound inference that two one-root
        # intervals inside the same overlap necessarily contain the same root.
        from .signs import sign_at_sample

        if (
            left_context.root_count(overlap_left, overlap_right) == 1
            and right_context.root_count(overlap_left, overlap_right) == 1
            and sign_at_sample(
                left_context.polynomial,
                (*right_context.parent_samples, right),
            )
            == 0
        ):
            relationship.common_root_pairs.add(pair)
            return True
    except (ArithmeticError, ValueError):
        return False
    return False


def _sample_is_root_at_rational(sample: Sample, point: sp.Rational) -> bool:
    if isinstance(sample, RationalSample):
        return sample.value == point
    if not isinstance(sample, AlgebraicRoot):
        return False
    if sample.interval.is_point():
        return sample.interval.left == point
    if point < sample.interval.left or point > sample.interval.right:
        return False
    if sample.fiber_context is not None:
        try:
            # A fiber polynomial can vanish at an interval endpoint because a
            # *different* root of the same family lies there.  Equality with a
            # ranked section therefore requires both endpoint vanishing and the
            # section's absolute root rank at that boundary.
            return (
                sample.fiber_context.sign_at_rational(point) == 0
                and sample.fiber_context.roots_left_of(point) == sample.root_index + 1
            )
        except (ArithmeticError, ValueError):
            return False
    if certified_zero(sample.polynomial.eval(point)) is not True:
        return False
    # Ordinary isolating intervals may likewise use another root as an endpoint.
    # Ask for the selected root, rather than equating polynomial vanishing with
    # equality of this root object.
    try:
        return bool(sp.CRootOf(sample.polynomial.as_expr(), sample.root_index) == point)
    except (NotImplementedError, ValueError, TypeError):
        return False


def compare_samples(left: Sample, right: Sample) -> int:
    """Compare two explicit sample objects exactly enough for CAD ordering."""

    if not isinstance(left, (RationalSample, AlgebraicRoot)) or not isinstance(
        right, (RationalSample, AlgebraicRoot)
    ):
        raise TypeError("compare_samples requires explicit semialg Sample objects")
    canonical, orientation = _structural_cache_pair(left, right)
    cached = CACHE.comparisons.get(canonical)
    if cached is not None:
        CACHE.stats.comparison_hits += 1
        return orientation * cached
    CACHE.stats.comparison_misses += 1
    if isinstance(left, AlgebraicRoot) and isinstance(right, AlgebraicRoot):
        left_context = left.fiber_context
        right_context = right.fiber_context
        same_family = (
            left_context is None and right_context is None and left.polynomial == right.polynomial
        ) or (
            left_context is not None
            and right_context is not None
            and left_context.polynomial == right_context.polynomial
            and tuple(sample_identity_key(s) for s in left_context.parent_samples)
            == tuple(sample_identity_key(s) for s in right_context.parent_samples)
            and left_context.fiber_variable == right_context.fiber_variable
        )
        if same_family:
            result = int(left.root_index > right.root_index) - int(
                left.root_index < right.root_index
            )
            CACHE.comparisons.put(canonical, orientation * result)
            return result
    if isinstance(left, RationalSample) and isinstance(right, RationalSample):
        result = int(bool(left.value > right.value)) - int(bool(left.value < right.value))
        CACHE.comparisons.put(canonical, orientation * result)
        return result

    # Root objects can carry an older certificate than the persistent root-index
    # partition.  Refresh their interval state before deciding whether another
    # comparison refinement is necessary.
    if isinstance(left, AlgebraicRoot) and left.fiber_context is not None:
        cached = left.fiber_context.best_root_intvl(left.root_index)
        if cached is not None and cached.width < left.interval.width:
            left = AlgebraicRoot(
                left.polynomial,
                cached,
                left.root_index,
                left.multiplicity,
                left.root_expr,
                left.fiber_context,
            )
    if isinstance(right, AlgebraicRoot) and right.fiber_context is not None:
        cached = right.fiber_context.best_root_intvl(right.root_index)
        if cached is not None and cached.width < right.interval.width:
            right = AlgebraicRoot(
                right.polynomial,
                cached,
                right.root_index,
                right.multiplicity,
                right.root_expr,
                right.fiber_context,
            )

    interval_order = left.interval.strict_order(right.interval)
    if interval_order is not None:
        CACHE.comparisons.put(canonical, orientation * interval_order)
        return interval_order
    if left.interval.right == right.interval.left:
        point = left.interval.right
        left_at_point = _sample_is_root_at_rational(left, point)
        right_at_point = _sample_is_root_at_rational(right, point)
        result = 0 if left_at_point and right_at_point else -1
        CACHE.comparisons.put(canonical, orientation * result)
        return result
    if right.interval.right == left.interval.left:
        point = right.interval.right
        left_at_point = _sample_is_root_at_rational(left, point)
        right_at_point = _sample_is_root_at_rational(right, point)
        result = 0 if left_at_point and right_at_point else 1
        CACHE.comparisons.put(canonical, orientation * result)
        return result

    # A separator established by an earlier cross-family comparison is part of
    # the persistent partition topology.  Consume it before resultant/Sturm work.
    if isinstance(left, AlgebraicRoot) and isinstance(right, AlgebraicRoot):
        for lower, upper, result in ((left, right, -1), (right, left, 1)):
            separator = CACHE.sector_separators.get(
                ("sector-between", sample_identity_key(lower), sample_identity_key(upper))
            )
            if separator is not None:
                CACHE.comparisons.put(canonical, orientation * result)
                return result

    # Resolve cross-family sections at family/root-index granularity.  If the
    # family resultant vanishes, first certify whether this particular pair is
    # one of the common sections.  Otherwise (and for all remaining distinct
    # pairs) order the root ranks by rational Sturm separators, without
    # repeatedly contracting the individual root certificates.
    if isinstance(left, AlgebraicRoot) and isinstance(right, AlgebraicRoot):
        related = _fiber_family_relationship(left, right)
        if related is not None:
            relationship, _ = related
            if (
                relationship.resultant_zero is True
                and max(left.interval.left, right.interval.left)
                < min(left.interval.right, right.interval.right)
                and _fiber_common_root_is_cached_or_certified(left, right)
            ):
                CACHE.comparisons.put(canonical, 0)
                return 0
            family_order = _fiber_root_order_from_family_cert(left, right)
            if family_order is not None:
                CACHE.comparisons.put(canonical, orientation * family_order)
                return family_order

    # Non-tower projection polynomials can also encode the same section.
    if isinstance(left, AlgebraicRoot) and isinstance(right, AlgebraicRoot):
        try:
            var = left.variable
            right_expr = right.polynomial.as_expr().subs(right.variable, var)
            left_poly = sp.Poly(left.polynomial.as_expr(), var, extension=True)
            right_poly = sp.Poly(right_expr, var, extension=True)
            common = sp.gcd(left_poly, right_poly)
            overlap_left = max(left.interval.left, right.interval.left)
            overlap_right = min(left.interval.right, right.interval.right)
            if (
                common.degree() > 0
                and overlap_left <= overlap_right
                and _ordinary_root_count(common, overlap_left, overlap_right) > 0
            ):
                CACHE.comparisons.put(canonical, 0)
                return 0
        except (
            NotImplementedError,
            sp.PolynomialError,
            sp.polys.polyerrors.DomainError,
            ValueError,
            TypeError,
        ):
            pass
    # Refine certified isolating intervals before reconstructing algebraic
    # expressions. CAD roots normally separate after a few rational bisections.
    refined_left = left
    refined_right = right
    for _ in range(64):
        candidates = []
        if isinstance(refined_left, AlgebraicRoot):
            candidates.append((refined_left.interval.width, "left"))
        if isinstance(refined_right, AlgebraicRoot):
            candidates.append((refined_right.interval.width, "right"))
        if not candidates:
            break
        # Refine only the interval currently obstructing separation most.  If
        # that family cannot refine tower-natively, try the other side before
        # declaring the order unresolved.
        candidates.sort(reverse=True, key=lambda item: item[0])
        changed = False
        for _, side in candidates:
            try:
                if side == "left":
                    updated = refine_root_once(refined_left, reason="separation")
                    changed = updated.interval != refined_left.interval
                    refined_left = updated
                else:
                    updated = refine_root_once(refined_right, reason="separation")
                    changed = updated.interval != refined_right.interval
                    refined_right = updated
            except ArithmeticError:
                continue
            if changed:
                break
        interval_order = refined_left.interval.strict_order(refined_right.interval)
        if interval_order is not None:
            if changed:
                mark_refinement_decisive("separation")
            CACHE.comparisons.put(canonical, orientation * interval_order)
            return interval_order
        if refined_left.interval.right == refined_right.interval.left:
            point = refined_left.interval.right
            if _sample_is_root_at_rational(refined_left, point) and _sample_is_root_at_rational(
                refined_right, point
            ):
                if changed:
                    mark_refinement_decisive("separation")
                CACHE.comparisons.put(canonical, 0)
                return 0
        if refined_right.interval.right == refined_left.interval.left:
            point = refined_right.interval.right
            if _sample_is_root_at_rational(refined_left, point) and _sample_is_root_at_rational(
                refined_right, point
            ):
                if changed:
                    mark_refinement_decisive("separation")
                CACHE.comparisons.put(canonical, 0)
                return 0
        CACHE.stats.comparison_refinements += 1
        if not changed:
            break

    # Parent-aware CAD roots must remain inside their algebraic tower.  Falling
    # back to expression reconstruction here recreates a global primitive
    # element and defeats the fiber certificate.
    has_fiber_root = any(
        isinstance(sample, AlgebraicRoot) and sample.fiber_context is not None
        for sample in (refined_left, refined_right)
    )
    if has_fiber_root:
        raise ArithmeticError(
            "fiber-root order remained unresolved after certified interval refinement"
        )

    diff = sample_to_expr(refined_left) - sample_to_expr(refined_right)
    if certified_zero(diff) is True:
        CACHE.comparisons.put(canonical, 0)
        return 0
    sign = sp.sign(diff)
    if sign in (-1, 1):
        result = int(sign)
        CACHE.comparisons.put(canonical, orientation * result)
        return result
    result = exact_algebraic_sign(diff)
    if result is None:
        raise ValueError("could not compare algebraic samples exactly")
    CACHE.comparisons.put(canonical, orientation * result)
    return result


def sort_samples(samples: list[Sample] | tuple[Sample, ...]) -> tuple[Sample, ...]:
    out = list(samples)
    for i in range(1, len(out)):
        item = out[i]
        pos = i
        while pos > 0 and compare_samples(out[pos - 1], item) > 0:
            out[pos] = out[pos - 1]
            pos -= 1
        out[pos] = item
    unique: list[Sample] = []
    for sample in out:
        if not unique or compare_samples(unique[-1], sample) != 0:
            unique.append(sample)
    return tuple(unique)
