from __future__ import annotations

import itertools
from dataclasses import dataclass
from functools import lru_cache

import sympy as sp
from sympy import S

from ..exactness import certify_equal, certify_zero
from ..exceptions import SYMBOLIC_QUERY_ERRORS
from ..results import Monotonicity, TruthValue
from .common import condition_to_set, exact_order_real_points


@dataclass(frozen=True)
class MonotonicDecomposition:
    expression: sp.Basic
    variable: sp.Symbol
    inner: sp.Basic
    outer: sp.Lambda
    direction: int
    condition: sp.Basic = S.true
    inverse: sp.Lambda | None = None
    strict: bool = True
    continuous: bool = True
    source: str = "structural-monotonic-decomposition"

    def reconstruct(self):
        return sp.simplify(self.outer(self.inner))


@dataclass(frozen=True)
class MonotonicInterval:
    interval: sp.Set
    direction: Monotonicity
    left_limit: sp.Basic | None = None
    right_limit: sp.Basic | None = None
    strict: bool = True
    continuous: bool = True
    left_attained: bool | None = None
    right_attained: bool | None = None
    evidence: tuple = ()

    @property
    def increasing(self):
        return self.direction in (
            Monotonicity.STRICTLY_INCREASING,
            Monotonicity.NONDECREASING,
        )

    @property
    def decreasing(self):
        return self.direction in (
            Monotonicity.STRICTLY_DECREASING,
            Monotonicity.NONINCREASING,
        )


@dataclass(frozen=True)
class IndexedMonotonicPartition:
    domain: sp.Set
    boundary_geometry: tuple
    # Empty means the geometry proves alternation but not its phase.  This is
    # deliberately different from assuming that the leftmost cell increases.
    direction_pattern: tuple[Monotonicity, ...] = ()
    complete: bool = False
    first_direction: Monotonicity | None = None
    source: str = "indexed-derivative-zero-geometry"


@dataclass(frozen=True)
class MonotonicPartition:
    domain: sp.Set
    cells: tuple[MonotonicInterval, ...] = ()
    indexed: IndexedMonotonicPartition | None = None
    complete: bool = False
    source: str = "monotonic-partition"
    derivative: sp.Basic | None = None
    critical_points: tuple = ()
    property_decomposition: object = None
    critical_points_complete: bool = False

    @property
    def finite(self):
        return self.indexed is None

    @property
    def direction_changes(self):
        vals = [c.direction for c in self.cells]
        n = 0
        for a, b in itertools.pairwise(vals):
            if a != b:
                n += 1
        return n

    @property
    def has_indexed_oscillation(self):
        return bool(
            self.indexed is not None
            and self.indexed.complete
            and self.indexed.boundary_geometry
        )

    def same_strict_direction(self):
        if not self.complete or not self.cells:
            return None
        dirs = {c.direction for c in self.cells}
        if dirs == {Monotonicity.STRICTLY_INCREASING}:
            return Monotonicity.STRICTLY_INCREASING
        if dirs == {Monotonicity.STRICTLY_DECREASING}:
            return Monotonicity.STRICTLY_DECREASING
        return None


@dataclass(frozen=True)
class NonMonotonicityWitness:
    points: tuple
    values: tuple
    exact: bool
    source: str = "exact-sample-ordering-witness"


@dataclass(frozen=True)
class InverseSolutionCertificate:
    equation: sp.Basic
    domain: sp.Set
    solutions: sp.Set
    complete: bool
    partition: MonotonicPartition | None = None
    source: str = "monotonic-partition-inverse"
    evidence: tuple = ()
    isolated_intervals: tuple = ()


@dataclass(frozen=True)
class PartitionRange:
    image: sp.Set | None
    complete: bool
    partition: MonotonicPartition
    source: str = "monotonic-partition-range"


def _monotone_head(expr, X):
    """Return (direction, condition, inverse) for high-value primitive heads."""
    f = expr.func
    args = expr.args
    slot = sp.Dummy("u", real=True)
    if len(args) != 1 or args[0] != X:
        return None
    if f in (sp.exp, sp.log, sp.sinh, sp.tanh, sp.atan, sp.asinh, sp.erf, sp.erfi):
        cond = X > 0 if f is sp.log else S.true
        inv = {
            sp.exp: sp.log,
            sp.log: sp.exp,
            sp.sinh: sp.asinh,
            sp.tanh: sp.atanh,
            sp.atan: sp.tan,
            sp.asinh: sp.sinh,
            sp.erf: sp.erfinv,
        }.get(f)
        invlam = sp.Lambda(slot, inv(slot)) if inv is not None else None
        return 1, cond, invlam
    if f in (sp.acos, sp.erfc):
        cond = sp.And(X >= -1, X <= 1) if f is sp.acos else S.true
        inv = sp.cos if f is sp.acos else getattr(sp, "erfcinv", None)
        return -1, cond, (sp.Lambda(slot, inv(slot)) if inv else None)
    if f is sp.asin:
        return 1, sp.And(X >= -1, X <= 1), sp.Lambda(slot, sp.sin(slot))
    if f is sp.acosh:
        return 1, X >= 1, sp.Lambda(slot, sp.cosh(slot))
    if f is sp.atanh:
        return 1, sp.And(X > -1, X < 1), sp.Lambda(slot, sp.tanh(slot))
    if f is sp.LambertW and len(args) == 1:
        return 1, X >= -1 / sp.E, sp.Lambda(slot, slot * sp.exp(slot))
    return None


def _monotonic_decomposition_result(expr, var, *, strict=True, continuous=True):
    """Find ``expr = F(inner)`` with a certified monotone outer ``F``.

    By default return the mathematical decomposition
    ``(inner, outer, direction, condition)``.  With ``return_result=True``
    return the full :class:`MonotonicDecomposition` certificate object.
    """
    expr = sp.sympify(expr)
    var = sp.sympify(var)
    if not expr.has(var):
        return MonotonicDecomposition(
            expr,
            var,
            S.Zero,
            sp.Lambda(sp.Dummy("u"), expr),
            1,
            S.true,
            None,
            False,
            True,
            "constant",
        )
    u = sp.Dummy("u", real=True)
    if expr.is_Add:
        c, dep = expr.as_independent(var, as_Add=True)
        if c != 0 and dep != expr:
            return MonotonicDecomposition(
                expr,
                var,
                dep,
                sp.Lambda(u, u + c),
                1,
                S.true,
                sp.Lambda(u, u - c),
                strict,
                continuous,
                "affine-shift",
            )
    if expr.is_Mul:
        c, dep = expr.as_independent(var, as_Add=False)
        if dep != expr and c.is_real is True and certify_zero(c) is False:
            if c.is_positive is True:
                return MonotonicDecomposition(
                    expr,
                    var,
                    dep,
                    sp.Lambda(u, c * u),
                    1,
                    S.true,
                    sp.Lambda(u, u / c),
                    strict,
                    continuous,
                    "positive-scale",
                )
            if c.is_negative is True:
                return MonotonicDecomposition(
                    expr,
                    var,
                    dep,
                    sp.Lambda(u, c * u),
                    -1,
                    S.true,
                    sp.Lambda(u, u / c),
                    strict,
                    continuous,
                    "negative-scale",
                )
    if isinstance(expr, sp.Pow):
        base, exponent = expr.args
        if base.has(var) and not exponent.has(var):
            if exponent.is_integer and exponent.is_positive:
                if exponent.is_odd:
                    return MonotonicDecomposition(
                        expr,
                        var,
                        base,
                        sp.Lambda(u, u**exponent),
                        1,
                        S.true,
                        sp.Lambda(u, sp.real_root(u, exponent)),
                        strict,
                        continuous,
                        "odd-positive-power",
                    )
                if exponent.is_even and exponent != 2:
                    return MonotonicDecomposition(
                        expr,
                        var,
                        base**2,
                        sp.Lambda(u, u ** (exponent / 2)),
                        1,
                        S.true,
                        None,
                        strict,
                        continuous,
                        "even-power-square-reduction",
                    )
            if exponent.is_Rational:
                d = 1 if exponent.is_positive else -1
                cond = base >= 0 if exponent.is_positive else base > 0
                return MonotonicDecomposition(
                    expr,
                    var,
                    base,
                    sp.Lambda(u, u**exponent),
                    d,
                    cond,
                    None,
                    strict,
                    continuous,
                    "rational-power",
                )
        if (
            not base.has(var)
            and exponent.has(var)
            and base.is_positive is True
            and base != 1
        ):
            d = 1 if (base - 1).is_positive else (-1 if (base - 1).is_negative else 0)
            if d:
                return MonotonicDecomposition(
                    expr,
                    var,
                    exponent,
                    sp.Lambda(u, base**u),
                    d,
                    S.true,
                    sp.Lambda(u, sp.log(u, base)),
                    strict,
                    continuous,
                    "positive-base-exponential",
                )
    varying = [
        (i, a)
        for i, a in enumerate(getattr(expr, "args", ()))
        if getattr(a, "has", lambda *_: False)(var)
    ]
    if len(varying) == 1 and expr.func not in (sp.Add, sp.Mul, sp.Pow, sp.Piecewise):
        i, inner = varying[0]
        if inner != var:
            X = sp.Dummy("X", real=True)
            aa = list(expr.args)
            aa[i] = X
            outer_expr = expr.func(*aa)
            info = _monotone_head(outer_expr, X)
            if info:
                d, cond, inv = info
                return MonotonicDecomposition(
                    expr,
                    var,
                    inner,
                    sp.Lambda(X, outer_expr),
                    d,
                    cond.xreplace({X: inner}),
                    inv,
                    strict,
                    continuous,
                    "elementary-monotone-head",
                )
            try:
                from .monotonicity import _direct_registry

                m = _direct_registry(outer_expr, X)
                inc = (Monotonicity.STRICTLY_INCREASING, Monotonicity.NONDECREASING)
                dec = (Monotonicity.STRICTLY_DECREASING, Monotonicity.NONINCREASING)
                if m in inc + dec:
                    d = 1 if m in inc else -1
                    from .domain import domain

                    cond = domain(outer_expr, X)
                    cond = (
                        cond.xreplace({X: inner})
                        if isinstance(cond, sp.Basic)
                        else S.true
                    )
                    return MonotonicDecomposition(
                        expr,
                        var,
                        inner,
                        sp.Lambda(X, outer_expr),
                        d,
                        cond,
                        None,
                        m in (inc[0], dec[0]),
                        continuous,
                        "registry-monotone-head",
                    )
            except SYMBOLIC_QUERY_ERRORS:
                pass
    return None


def monotonic_decomposition(
    expr, var, *, strict=True, continuous=True, return_result=False
):
    """Find ``expr = F(inner)`` with a certified monotone outer ``F``.

    The ordinary result is the mathematical quadruple
    ``(inner, outer, direction, condition)``.  ``direction`` is ``1`` for an
    increasing outer map and ``-1`` for a decreasing one.  Use
    ``return_result=True`` for inverse information, provenance, and replayable
    certification metadata.
    """
    result = _monotonic_decomposition_result(
        expr, var, strict=strict, continuous=continuous
    )
    if return_result or result is None:
        return result
    return (result.inner, result.outer, result.direction, result.condition)


def _endpoint_limit(expr, var, v, left, open_, domain=S.true, assumptions=S.true):
    from ..proof_cache import property_proof_context

    return property_proof_context(
        expr, var, domain=domain, assumptions=assumptions
    ).endpoint_limit(v, direction=("+" if left else "-"), attained=not open_)


def _phase_for_indexed(expr, var, ds, pd, derivative=None):
    """Determine an indexed alternating phase only when an exact theorem/sample does."""
    # For an infinite family the zero geometry itself proves alternation but
    # usually does not identify which sign occurs on the leftmost unbounded
    # cell.  A finite endpoint cell can sometimes be sampled exactly.
    if not isinstance(ds, sp.Interval) or not pd.boundary_geometry:
        return None
    try:
        from ..proof_cache import cached_derivative

        d = (
            sp.factor(derivative)
            if derivative is not None
            else cached_derivative(expr, var, order=1, factored=True)
        )
    except SYMBOLIC_QUERY_ERRORS:
        return None
    # Airy phase is known at x=0 exactly.  The cell immediately left of zero
    # therefore fixes the parity without numerical root evaluation.
    if expr.func in (sp.airyai, sp.airybi) and ds.end == 0:
        val = sp.simplify(d.subs(var, 0))
        if val.is_positive is True:
            return Monotonicity.STRICTLY_INCREASING
        if val.is_negative is True:
            return Monotonicity.STRICTLY_DECREASING
    # Otherwise do not guess a phase from structural zero ordering.
    return None


@lru_cache(maxsize=512)
def _monotonic_partition_on_set(expr, var, domain_set):
    """Build an exact monotonic partition on an already-normalized real set.

    The routine prefers derivative zero geometry, then falls back to finite
    critical-point analysis. Completeness is retained explicitly so consumers do
    not mistake a partial decomposition for a global proof.
    """
    ds = domain_set
    expr = sp.sympify(expr)
    var = sp.sympify(var)
    if not isinstance(ds, sp.Set):
        return MonotonicPartition(S.Reals, complete=False, source="unknown-domain")
    derivative = None
    try:
        from ..proof_cache import property_proof_context

        derivative = property_proof_context(expr, var, domain=ds).derivative(
            1, factored=True
        )
    except SYMBOLIC_QUERY_ERRORS:
        pass
    from ..transcendental.zero_geometry import property_interval_decomposition

    pd = property_interval_decomposition(
        expr,
        var,
        derivative_order=1,
        domain=ds,
        positive_value=Monotonicity.STRICTLY_INCREASING,
        negative_value=Monotonicity.STRICTLY_DECREASING,
        derivative=derivative,
    )
    if pd.cells:
        cells = []
        crit = ()
        for g in pd.boundary_geometry:
            if getattr(g, "family", None) == "finite-derivative-zeros" and isinstance(
                getattr(g, "support", None), sp.FiniteSet
            ):
                q = exact_order_real_points(tuple(g.support))
                if q is not None:
                    crit = tuple(
                        z
                        for z in q
                        if not (isinstance(ds, sp.Interval) and z in (ds.start, ds.end))
                    )
                break
        for j, c in enumerate(pd.cells):
            it = getattr(c, "interval", None) or getattr(c, "domain", None)
            value = getattr(c, "value", None) or getattr(c, "property", None)
            if not isinstance(it, sp.Interval) or value is None:
                continue
            left_outer = j == 0 and isinstance(ds, sp.Interval) and it.start == ds.start
            right_outer = (
                j == len(pd.cells) - 1
                and isinstance(ds, sp.Interval)
                and it.end == ds.end
            )
            la = (
                (not ds.left_open)
                if left_outer and isinstance(ds, sp.Interval)
                else True
            )
            ra = (
                (not ds.right_open)
                if right_outer and isinstance(ds, sp.Interval)
                else True
            )
            cells.append(
                MonotonicInterval(
                    it,
                    value,
                    _endpoint_limit(expr, var, it.start, True, not la, ds),
                    _endpoint_limit(expr, var, it.end, False, not ra, ds),
                    True,
                    True,
                    la,
                    ra,
                    (c,),
                )
            )
        crit_complete = bool(
            pd.boundary_geometry
            and any(
                getattr(g, "family", None) == "finite-derivative-zeros"
                and getattr(g, "complete_on_support", False)
                for g in pd.boundary_geometry
            )
        )
        return MonotonicPartition(
            ds,
            tuple(cells),
            complete=bool(pd.complete and len(cells) == len(pd.cells)),
            source="finite-derivative-zero-partition",
            derivative=derivative,
            critical_points=crit,
            property_decomposition=pd,
            critical_points_complete=crit_complete,
        )
    if pd.boundary_geometry and pd.alternating:
        phase = _phase_for_indexed(expr, var, ds, pd, derivative)
        pattern = (
            ()
            if phase is None
            else (
                phase,
                Monotonicity.STRICTLY_DECREASING
                if phase is Monotonicity.STRICTLY_INCREASING
                else Monotonicity.STRICTLY_INCREASING,
            )
        )
        idx = IndexedMonotonicPartition(
            ds, pd.boundary_geometry, pattern, pd.complete, phase
        )
        return MonotonicPartition(
            ds,
            indexed=idx,
            complete=pd.complete,
            source="indexed-derivative-zero-geometry",
            derivative=derivative,
            property_decomposition=pd,
            critical_points_complete=False,
        )
    return MonotonicPartition(
        ds,
        complete=False,
        derivative=derivative,
        property_decomposition=pd,
        critical_points_complete=False,
    )


def _monotonic_partition_result(expr, var, *, domain=True):
    """Return the rich cached exact finite or indexed monotonic partition."""
    expr = sp.sympify(expr)
    var = sp.sympify(var)
    if isinstance(domain, sp.Set):
        ds = domain
    else:
        if domain is True:
            from .domain import domain

            cond = domain(expr, var)
        else:
            cond = sp.sympify(domain)
        ds = condition_to_set(cond, var)
    if not isinstance(ds, sp.Set):
        return MonotonicPartition(S.Reals, complete=False, source="unknown-domain")
    try:
        return _monotonic_partition_on_set(expr, var, ds)
    except TypeError:
        return _monotonic_partition_on_set.__wrapped__(expr, var, ds)


def monotonic_partition(expr, var, *, domain=True, return_result=False):
    """Partition a real domain into intervals of certified monotonicity.

    Complete finite partitions are returned by default as ordered
    ``(interval, direction)`` pairs.  If the result is incomplete or uses an
    indexed infinite description, the structured :class:`MonotonicPartition`
    is retained so completeness information is not lost.  Set
    ``return_result=True`` to always receive the structured result.
    """
    result = _monotonic_partition_result(expr, var, domain=domain)
    if return_result:
        return result
    if result.complete and result.indexed is None:
        return tuple((cell.interval, cell.direction) for cell in result.cells)
    return result


def clear_monotonic_partition_cache():
    # Preserve the public spelling while clearing all proof primitives that can
    # feed a partition; otherwise a "cleared" partition cache can immediately
    # repopulate from stale lower-level proofs.
    from ..proof_cache import clear_property_proof_cache

    clear_property_proof_cache()


def _cell_image(cell: MonotonicInterval):
    if not isinstance(cell.interval, sp.Interval):
        return None
    a, b = cell.left_limit, cell.right_limit
    if a is None or b is None:
        return None
    if cell.increasing:
        lo, hi = a, b
        la, ua = cell.left_attained, cell.right_attained
    elif cell.decreasing:
        lo, hi = b, a
        la, ua = cell.right_attained, cell.left_attained
    else:
        return None
    if la is None or ua is None:
        return None
    if lo is S.NegativeInfinity:
        la = False
    if hi is S.Infinity:
        ua = False
    try:
        return sp.Interval(lo, hi, left_open=not la, right_open=not ua)
    except SYMBOLIC_QUERY_ERRORS:
        return None


def range_from_monotonic_partition(expr, var, partition=None, *, domain=True):
    """Exact image set when a finite monotonic partition is complete."""
    part = partition
    if not isinstance(part, MonotonicPartition):
        part = monotonic_partition(expr, var, domain=domain, return_result=True)
    if not part.complete or part.indexed is not None or not part.cells:
        return PartitionRange(None, False, part)
    images = []
    for c in part.cells:
        im = _cell_image(c)
        if im is None:
            return PartitionRange(None, False, part)
        images.append(im)
    try:
        image = sp.Union(*images)
    except SYMBOLIC_QUERY_ERRORS:
        return PartitionRange(None, False, part)
    return PartitionRange(image, True, part)


def injectivity_from_partition(partition):
    """Three-valued injectivity conclusion from a compact or rich partition."""
    if isinstance(partition, tuple):
        directions = [
            item[1] for item in partition if isinstance(item, tuple) and len(item) == 2
        ]
        if len(directions) != len(partition) or not directions:
            return TruthValue.UNKNOWN
        if set(directions) in (
            {Monotonicity.STRICTLY_INCREASING},
            {Monotonicity.STRICTLY_DECREASING},
        ):
            return TruthValue.TRUE
        return TruthValue.UNKNOWN
    if not isinstance(partition, MonotonicPartition) or not partition.complete:
        return TruthValue.UNKNOWN
    if partition.indexed is not None:
        # Infinite simple alternating extrema force collisions on an interval.
        return (
            TruthValue.FALSE
            if partition.indexed.boundary_geometry
            else TruthValue.UNKNOWN
        )
    if not partition.cells:
        return TruthValue.UNKNOWN
    same = partition.same_strict_direction()
    if same is not None:
        return TruthValue.TRUE
    # Direction reversal at an attained continuous interior critical point
    # creates a strict local extremum and hence duplicate nearby values.
    for a, b in zip(partition.cells, partition.cells[1:], strict=False):
        if (
            a.direction != b.direction
            and a.continuous
            and b.continuous
            and a.right_attained
            and b.left_attained
        ):
            return TruthValue.FALSE
    return TruthValue.UNKNOWN


def _certified_inverse_cell(expr, var, value, cell):
    """Resolve one strict monotonic cell exactly or by a rational isolating interval."""
    from ..transcendental.real_sign import robust_real_sign
    from ..transcendental.roots import CertifiedRootInterval, RootCertificate

    interval = cell.interval
    if not isinstance(interval, sp.Interval) or not cell.strict:
        return None
    left = None if cell.left_limit is None else sp.simplify(cell.left_limit - value)
    right = None if cell.right_limit is None else sp.simplify(cell.right_limit - value)
    if left is None or right is None:
        return None
    sl = robust_real_sign(left)
    sr = robust_real_sign(right)
    if sl is None or sr is None:
        return None
    exact = []
    if sl == 0 and cell.left_attained:
        exact.append(interval.start)
    if sr == 0 and cell.right_attained:
        exact.append(interval.end)
    if exact:
        return sp.FiniteSet(*exact), (), True
    if sl == sr and sl in (-1, 1):
        return S.EmptySet, (), True
    if sl * sr != -1:
        return None
    if interval.start in (S.NegativeInfinity, S.Infinity) or interval.end in (
        S.NegativeInfinity,
        S.Infinity,
    ):
        return None

    h = sp.simplify(expr - value)
    a, b = interval.start, interval.end
    # Move to exact interior dyadic points when an endpoint is open or non-rational.
    aa = a
    bb = b
    left_sign = sl
    right_sign = sr
    for level in range(4, 28):
        if robust_real_sign(sp.simplify(h.subs(var, aa))) != left_sign:
            aa = sp.simplify(a + (b - a) / (2**level))
        if robust_real_sign(sp.simplify(h.subs(var, bb))) != right_sign:
            bb = sp.simplify(b - (b - a) / (2**level))
        if (
            robust_real_sign(sp.simplify(h.subs(var, aa))) == left_sign
            and robust_real_sign(sp.simplify(h.subs(var, bb))) == right_sign
        ):
            break
    else:
        return None
    for _ in range(32):
        m = sp.simplify((aa + bb) / 2)
        sm = robust_real_sign(sp.simplify(h.subs(var, m)))
        if sm == 0:
            return sp.FiniteSet(m), (), True
        if sm not in (-1, 1):
            break
        if sm == left_sign:
            aa = m
        else:
            bb = m
    iv = sp.Interval(aa, bb)
    cert = RootCertificate(
        interval,
        1,
        "certified_monotonic_inverse",
        (("target", value), ("endpoint_signs", sl, sr), ("isolating_interval", iv)),
        True,
        1,
    )
    return S.EmptySet, (CertifiedRootInterval(iv, cert),), True


def certified_monotonic_inverse(expr, var, value, *, domain=True, partition=None):
    """Solve ``expr == value`` using a shared monotonic partition certificate.

    Finite cells are solved independently, with uniqueness guaranteed on each
    strict cell.  Indexed partitions remain exact ``ConditionSet`` objects;
    this routine never numerically truncates an infinite family.
    """
    expr = sp.sympify(expr)
    var = sp.sympify(var)
    value = sp.sympify(value)
    part = partition
    if not isinstance(part, MonotonicPartition):
        part = monotonic_partition(expr, var, domain=domain, return_result=True)
    eq = sp.Eq(expr, value)
    if not part.complete:
        ds = part.domain if isinstance(part.domain, sp.Set) else S.Reals
        return InverseSolutionCertificate(
            eq,
            ds,
            sp.ConditionSet(var, eq, ds),
            False,
            part,
            source="incomplete-monotonic-partition",
        )
    if part.indexed is not None:
        return InverseSolutionCertificate(
            eq,
            part.domain,
            sp.ConditionSet(var, eq, part.domain),
            False,
            part,
            source="indexed-monotonic-partition",
        )
    # Prefer one global exact solve; splitting by monotonic cells is a fallback,
    # not a reason to multiply solver work when SymPy already has a closed form.
    try:
        global_sol = sp.solveset(eq, var, domain=part.domain)
    except SYMBOLIC_QUERY_ERRORS:
        global_sol = sp.ConditionSet(var, eq, part.domain)
    if not isinstance(global_sol, sp.ConditionSet) and not getattr(
        global_sol, "has", lambda *_: False
    )(sp.ConditionSet):
        return InverseSolutionCertificate(
            eq,
            part.domain,
            global_sol,
            True,
            part,
            evidence=(part,),
            source="global-solve+monotonic-partition",
        )
    sols = []
    isolated = []
    complete = True
    for cell in part.cells:
        if not isinstance(cell.interval, sp.Interval):
            complete = False
            continue
        image = _cell_image(cell)
        if image is not None:
            try:
                contained = image.contains(value)
                if contained is S.false or contained is False:
                    continue
            except SYMBOLIC_QUERY_ERRORS:
                pass
        try:
            s = sp.solveset(eq, var, domain=cell.interval)
        except SYMBOLIC_QUERY_ERRORS:
            s = sp.ConditionSet(var, eq, cell.interval)
        if isinstance(s, sp.ConditionSet) or getattr(s, "has", lambda *_: False)(
            sp.ConditionSet
        ):
            certified = _certified_inverse_cell(expr, var, value, cell)
            if certified is None:
                complete = False
                sols.append(s)
            else:
                exact_part, isolated_part, cell_complete = certified
                sols.append(exact_part)
                isolated.extend(isolated_part)
                complete = complete and cell_complete
        else:
            sols.append(s)
    # Include attained shared/end boundaries exactly once; solveset on open
    # cells intentionally excludes them.
    boundaries = []
    if isinstance(part.domain, sp.Interval):
        if not part.domain.left_open:
            boundaries.append(part.domain.start)
        if not part.domain.right_open:
            boundaries.append(part.domain.end)
    boundaries.extend(part.critical_points)
    for q in boundaries:
        try:
            if certify_equal(expr.subs(var, q), value) is True:
                sols.append(sp.FiniteSet(q))
        except SYMBOLIC_QUERY_ERRORS:
            pass
    sol = sp.Union(*sols) if sols else S.EmptySet
    return InverseSolutionCertificate(
        eq,
        part.domain,
        sol,
        complete,
        part,
        evidence=(part,),
        isolated_intervals=tuple(isolated),
    )


def quick_nonmonotonicity_witness(expr, var, domain):
    """Try exact rational points and return only an exactly verified violation."""
    ds = domain if isinstance(domain, sp.Set) else condition_to_set(domain, var)
    if (
        not isinstance(ds, sp.Interval)
        or ds.start in (S.NegativeInfinity, S.Infinity)
        or ds.end in (S.NegativeInfinity, S.Infinity)
    ):
        return None
    a, b = ds.start, ds.end
    pts = tuple(
        sp.simplify(a + (b - a) * q)
        for q in (
            sp.Rational(1, 5),
            sp.Rational(2, 5),
            sp.Rational(3, 5),
            sp.Rational(4, 5),
        )
    )
    vals = tuple(sp.simplify(expr.subs(var, p)) for p in pts)
    signs = []
    for u, v in itertools.pairwise(vals):
        d = sp.simplify(v - u)
        if d.is_positive is True:
            signs.append(1)
        elif d.is_negative is True:
            signs.append(-1)
        elif certify_zero(d) is True:
            signs.append(0)
        else:
            return None
    nz = {s for s in signs if s}
    if len(nz) > 1:
        return NonMonotonicityWitness(pts, vals, True)
    return None
