"""Certified correspondence of ordered real branches at a collision fiber."""

from dataclasses import dataclass

import sympy as sp

from ..algebraic.cache import CACHE, poly_key, sample_identity_key
from ..algebraic.groebner_utils import algebraic_coefficient_field
from ..algebraic.roots import isolate_real_roots, root_multiplicity
from ..algebraic.sample_points import choose_sector_sample
from ..algebraic.samples import FiberRootContext, sample_to_expr
from ..algebraic.signs import sign_at_sample
from ..cache_utils import cache_disabled
from ..cad_algorithms.decomposition import _roots_in_parent_tower
from ..cad_algorithms.resource_limits import current_budget
from ..errors import ExactComputationUnavailable
from ..exact_arithmetic import exact_sign
from ..internal_symbols import fresh_dummy


@dataclass(frozen=True)
class FiberSpectrum:
    """Actual degree and ordered real-root multiplicities on a selected fiber."""

    degree: int
    multiplicities: tuple[int, ...]


@dataclass(frozen=True)
class BranchCorrespondenceCertificate:
    """An ordered branch limit proved by exact spectra and one-sided Sturm counts."""

    polynomial: sp.Expr
    variable: sp.Symbol
    base_variables: tuple[sp.Symbol, ...]
    source_samples: tuple[object, ...]
    target_values: tuple[sp.Expr, ...]
    source_index: int
    source_spectrum: FiberSpectrum
    target_roots: tuple[sp.Expr, ...]
    target_multiplicities: tuple[int, ...]
    target_index: int
    barriers: tuple[sp.Rational, ...] = ()
    continuation_ranks: tuple[int, ...] = ()

    @property
    def value(self):
        """Return the finite limit, or a signed point at infinity.

        The indices -1 and len(target_roots) select the two oriented ends of
        the affine line. Finite endpoint multiplicities stay in affine order.
        """
        if self.target_index == -1:
            return -sp.oo
        if self.target_index == len(self.target_roots):
            return sp.oo
        return self.target_roots[self.target_index]


def _source_spectrum(polynomial, variables, samples):
    """Compute multiplicities in the selected fiber, not over the generic field."""
    poly = sp.Poly(polynomial, *variables)
    key = (
        "closure-source-spectrum",
        poly_key(poly),
        tuple(sample_identity_key(s) for s in samples),
    )
    cached = CACHE.compiled_polys.get(key)
    if isinstance(cached, FiberSpectrum):
        return cached
    fiber = variables[-1]
    if not samples:
        roots = isolate_real_roots(sp.Poly(polynomial, fiber, extension=True))
        result = FiberSpectrum(poly.degree(fiber), tuple(root.multiplicity for root in roots))
    else:
        context = FiberRootContext(poly, variables, samples, fiber)
        if context.is_nullified_at_parent():
            raise ExactComputationUnavailable("branch selector nullifies on its source fiber")
        coefficients = sp.Poly(polynomial, fiber).all_coeffs()
        degree = len(coefficients) - 1
        for coefficient in coefficients:
            if sign_at_sample(sp.Poly(coefficient, *variables[:-1], domain="EX"), samples):
                break
            degree -= 1
        roots = _roots_in_parent_tower(context)
        multiplicities = [root_multiplicity(root) for root in roots]
        result = FiberSpectrum(degree, tuple(multiplicities))
    CACHE.compiled_polys.put(key, result)
    return result


def _possible_targets(source_weights, target_weights, selected):
    """Locate a branch in every admissible ordered collision allocation.

    Each endpoint cluster contains at most its multiplicity in real roots;
    the missing roots come in complex-conjugate pairs. Source branches retain
    their multiplicities and order on the source CAD cell. Dynamic programming
    enumerates allocations, not paths or numerical samples near the boundary.
    """
    count = len(source_weights)
    totals = [0]
    for weight in source_weights:
        totals.append(totals[-1] + weight)

    def ends(cluster, start):
        maximum = target_weights[cluster]
        for end in range(start, count + 1):
            real_weight = totals[end] - totals[start]
            if real_weight > maximum:
                break
            if (maximum - real_weight) % 2 == 0:
                yield end

    memo = {}

    def feasible(cluster, start):
        key = (cluster, start)
        if key in memo:
            return memo[key]
        if cluster == len(target_weights):
            return start == count
        memo[key] = any(feasible(cluster + 1, end) for end in ends(cluster, start))
        return memo[key]

    frontier = {0}
    targets = set()
    for cluster in range(len(target_weights)):
        next_frontier = set()
        for start in frontier:
            for end in ends(cluster, start):
                if feasible(cluster + 1, end):
                    next_frontier.add(end)
                    if start <= selected < end:
                        targets.add(cluster)
        frontier = next_frontier
    return targets


def _sturm_continuation(polynomial, variable, base, source, target, endpoint, source_count, index):
    """Count endpoint clusters by exact Sturm signs as a parameter tends to zero.

    The segment from a one-dimensional base interval's endpoint to its sample
    stays inside that interval. Signs over K(t) have certified eventual values
    given by their first nonzero Taylor coefficients. No numerical perturbation
    or epsilon grid is used to choose a branch.
    """
    parameter = fresh_dummy("branch_parameter", real=True)
    path = sp.expand(polynomial.subs(base, target + parameter * (sample_to_expr(source) - target)))
    ground = sp.Poly(path, variable, parameter, extension=True).domain
    family = sp.Poly(path, variable, domain=ground.frac_field(parameter))
    # The p,p' Sturm chain counts distinct roots even with persistent multiple
    # factors. Avoid repeated square-free/monic conversions over K(t), which
    # can inflate algebraic units without changing the variation count.
    sequence = [family, family.diff()]
    budget = current_budget()
    while not sequence[-1].is_zero:
        if budget is not None:
            budget.check()
        remainder = -sequence[-2].rem(sequence[-1])
        if remainder.is_zero:
            break
        sequence.append(remainder)

    def germ_sign(coefficient):
        # Keep K(t) coefficients native: converting a zero coefficient through
        # an expanded radical expression can create an artificial 0/0 before
        # the common number field has reduced it.
        signs = []
        for part in (coefficient.numer, coefficient.denom):
            if not part:
                signs.append(0)
            else:
                first = min(part, key=lambda monomial: monomial[0])
                signs.append(exact_sign(ground.to_sympy(part[first])))
        return signs[0] * signs[1]

    def variation(value=None, negative=False):
        signs = []
        for member in sequence:
            if value is None:
                sign = germ_sign(member.rep.LC())
                if negative and member.degree() % 2:
                    sign = -sign
            else:
                sign = germ_sign(member.rep.eval(family.domain.convert(value)))
            if sign:
                signs.append(sign)
        return sum(a != b for a, b in zip(signs, signs[1:], strict=False))

    minus = variation(negative=True)
    if minus - variation() != source_count:
        raise ExactComputationUnavailable("source spectrum does not match the continuation germ")
    endpoint_roots = isolate_real_roots(endpoint)
    barriers = tuple(
        choose_sector_sample(left, right).value
        for left, right in zip((None, *endpoint_roots), (*endpoint_roots, None), strict=True)
    )
    ranks = tuple(minus - variation(barrier) for barrier in barriers)
    if index < ranks[0] or index >= ranks[-1]:
        # In the reciprocal chart u=1/y, these branches meet u=0. The
        # nonzero endpoint polynomial excludes any other bounded limit; the
        # exterior Sturm counts distinguish the two real orientations.
        reciprocal = (
            endpoint.as_expr().subs(variable, 1 / parameter)
            * parameter ** sp.Poly(path, variable).degree()
        )
        if sp.cancel(reciprocal).subs(parameter, 0) != 0:
            raise ExactComputationUnavailable("escaping branch has no reciprocal endpoint")
        return (-1 if index < ranks[0] else len(endpoint_roots)), barriers, ranks
    targets = [j for j in range(len(endpoint_roots)) if ranks[j] <= index < ranks[j + 1]]
    if len(targets) != 1:
        raise ExactComputationUnavailable("Sturm continuation did not isolate the selected cluster")
    return targets[0], barriers, ranks


def certify_branch_correspondence(
    polynomial, variable, index, base_variables, source_samples, target_values
):
    """Certify an ordered branch limit from a delineable source CAD cell.

    Source multiplicities and endpoint degree are computed exactly. Ordered
    real/complex allocations certify unique targets; ambiguous one-dimensional
    base cases use exact Sturm continuation germs. The source cell must be
    delineable and the target in its closure. One base parameter supports
    degree loss and signed infinite limits through reciprocal charts. Unresolved
    higher-dimensional collisions raise ``ExactComputationUnavailable``.
    """
    base_variables = tuple(base_variables)
    source_samples = tuple(source_samples)
    target_values = tuple(map(sp.sympify, target_values))
    if len(source_samples) != len(base_variables) or len(target_values) != len(base_variables):
        raise ValueError("branch correspondence needs complete source and target prefixes")
    ground = algebraic_coefficient_field((polynomial,), (*base_variables, variable))
    coefficients = ground.poly_ring(*base_variables) if base_variables else ground
    polynomial = sp.Poly(polynomial, variable, domain=coefficients).primitive()[1].as_expr()
    source = _source_spectrum(polynomial, (*base_variables, variable), source_samples)
    if index < 0 or index >= len(source.multiplicities):
        raise ValueError("source branch index is outside its real-root spectrum")
    specialized = sp.Poly(
        sp.expand(polynomial.subs(dict(zip(base_variables, target_values, strict=True)))),
        variable,
        extension=True,
    )
    degree_loss = specialized.degree() != source.degree
    if specialized.is_zero or (degree_loss and len(base_variables) != 1):
        raise ExactComputationUnavailable("branch continuation needs refinement at degree loss")
    key = ("closure-endpoint-spectrum", poly_key(specialized))
    roots = CACHE.compiled_polys.get(key)
    if roots is None:
        roots = tuple(specialized.real_roots(multiple=False))
        CACHE.compiled_polys.put(key, roots)
    target_roots = tuple(value for value, _ in roots)
    target_weights = tuple(int(weight) for _, weight in roots)
    possible = (
        set() if degree_loss else _possible_targets(source.multiplicities, target_weights, index)
    )
    barriers = ranks = ()
    if len(possible) == 1:
        target_index = possible.pop()
    elif len(base_variables) == 1:
        target_index, barriers, ranks = _sturm_continuation(
            polynomial,
            variable,
            base_variables[0],
            source_samples[0],
            target_values[0],
            specialized,
            len(source.multiplicities),
            index,
        )
    else:
        raise ExactComputationUnavailable(
            "collision needs additional branch-continuation refinement"
        )
    return BranchCorrespondenceCertificate(
        polynomial,
        variable,
        base_variables,
        source_samples,
        target_values,
        int(index),
        source,
        target_roots,
        target_weights,
        target_index,
        barriers,
        ranks,
    )


def replay_branch_correspondence(certificate):
    """Recompute the spectra and conservation proof for the recorded statement."""
    try:
        with cache_disabled():
            expected = certify_branch_correspondence(
                certificate.polynomial,
                certificate.variable,
                certificate.source_index,
                certificate.base_variables,
                certificate.source_samples,
                certificate.target_values,
            )
    except (ExactComputationUnavailable, ValueError, TypeError, sp.PolynomialError):
        return False
    return expected == certificate
