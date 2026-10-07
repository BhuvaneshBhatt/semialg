from __future__ import annotations

from dataclasses import dataclass
from math import gcd, lcm

import sympy as sp
from sympy import S

from ..exactness import certify_identity_equal
from ..exceptions import SYMBOLIC_QUERY_ERRORS


def _canonical_period(period):
    p = sp.simplify(sp.sympify(period))
    if p.is_negative is True:
        p = -p
    return p


def _integer_ratio(a, b):
    r = sp.simplify(a / b)
    return int(r) if r.is_Integer is True and int(r) > 0 else None


def _translate_set(s, delta):
    delta = sp.sympify(delta)
    if delta == 0:
        return s
    if isinstance(s, sp.Interval):
        return sp.Interval(s.start + delta, s.end + delta, s.left_open, s.right_open)
    if isinstance(s, sp.FiniteSet):
        return sp.FiniteSet(*(sp.simplify(v + delta) for v in s))
    if isinstance(s, sp.Union):
        return sp.Union(*(_translate_set(a, delta) for a in s.args))
    return sp.imageset(lambda t: t + delta, s)


@dataclass(frozen=True)
class PeriodicSet:
    """Set obtained by translating a base set by integer multiples of a period."""

    base: sp.Set
    period: sp.Expr
    variable: sp.Symbol | None = None
    origin: sp.Expr = S.Zero

    def __post_init__(self):
        object.__setattr__(self, "base", sp.sympify(self.base))
        object.__setattr__(self, "period", _canonical_period(self.period))
        object.__setattr__(self, "origin", sp.sympify(self.origin))
        if self.period.is_positive is not True:
            raise ValueError("PeriodicSet period must be provably positive")
        object.__setattr__(
            self, "base", sp.Intersection(self.base, self.fundamental_interval)
        )

    @property
    def fundamental_interval(self):
        return sp.Interval.Ropen(self.origin, self.origin + self.period)

    def normalize_point(self, value):
        return sp.simplify(
            self.origin + sp.Mod(sp.sympify(value) - self.origin, self.period)
        )

    def contains(self, value):
        return sp.Contains(self.normalize_point(value), self.base)

    def as_condition(self, variable=None):
        v = sp.sympify(variable or self.variable)
        if v is None:
            raise ValueError("a variable is required")
        return self.contains(v)

    def intersect_fundamental(self, condition_or_set, variable=None):
        if isinstance(condition_or_set, sp.Set):
            s = condition_or_set
        else:
            sp.sympify(variable or self.variable)
            try:
                s = sp.sympify(condition_or_set).as_set()
            except (AttributeError, *SYMBOLIC_QUERY_ERRORS):
                return None
        return sp.Intersection(s, self.fundamental_interval)

    def with_period(self, new_period):
        """Represent the same periodic set on a larger commensurate period."""
        new_period = _canonical_period(new_period)
        n = _integer_ratio(new_period, self.period)
        if n is None:
            raise ValueError("new period must be a positive integer multiple")
        pieces = [_translate_set(self.base, k * self.period) for k in range(n)]
        return PeriodicSet(sp.Union(*pieces), new_period, self.variable, self.origin)

    def _align(self, other):
        if not isinstance(other, PeriodicSet):
            return None
        if certify_identity_equal(self.origin, other.origin) is False:
            return None
        p = common_period_values(self.period, other.period)
        if p is None:
            return None
        return self.with_period(p), other.with_period(p)

    def union(self, other):
        aligned = self._align(other)
        if aligned is None:
            return sp.Union(self.as_sympy_set(), other.as_sympy_set())
        a, b = aligned
        return PeriodicSet(
            sp.Union(a.base, b.base),
            a.period,
            self.variable or other.variable,
            a.origin,
        )

    def intersection(self, other):
        aligned = self._align(other)
        if aligned is None:
            return sp.Intersection(self.as_sympy_set(), other.as_sympy_set())
        a, b = aligned
        return PeriodicSet(
            sp.Intersection(a.base, b.base),
            a.period,
            self.variable or other.variable,
            a.origin,
        )

    def complement(self):
        return PeriodicSet(
            sp.Complement(self.fundamental_interval, self.base),
            self.period,
            self.variable,
            self.origin,
        )

    def affine_pullback(self, a, b=0, variable=None):
        """Exact pullback under ``t = a*x + b`` for nonzero real ``a``.

        The modular condition is exact for arbitrary real translations and
        does not depend on inspecting a finite window of period copies.
        """
        a = sp.sympify(a)
        b = sp.sympify(b)
        if a == 0 or a.is_real is False:
            raise ValueError("a must be nonzero real")
        if a.is_positive is not True and a.is_negative is not True:
            raise ValueError("sign of a must be known")
        p = sp.simplify(self.period / sp.Abs(a))
        x = sp.sympify(variable or self.variable or sp.Dummy("x", real=True))
        fi = sp.Interval.Ropen(S.Zero, p)
        # self.contains performs exact modular normalization in the source
        # period, so no finite translation window is needed.
        cond = self.contains(sp.simplify(a * x + b))
        base = sp.ConditionSet(x, cond, fi)
        # Simplify common finite-set cases without weakening the exact fallback.
        if isinstance(self.base, sp.FiniteSet):
            pts = []
            for q in self.base:
                pts.append(sp.simplify(sp.Mod((q - b) / a, p)))
            base = sp.FiniteSet(*pts)
        return PeriodicSet(base, p, x, S.Zero)

    def as_sympy_set(self):
        sp.Dummy("k", integer=True)
        return sp.imageset(lambda t, k: t + k * self.period, self.base, S.Integers)


@dataclass(frozen=True)
class PeriodicCondition:
    set: PeriodicSet

    def as_expr(self, variable=None):
        return self.set.as_condition(variable)

    def __sympy__(self):
        return self.as_expr()


def positive_period(expr, var):
    try:
        from sympy.calculus.util import periodicity

        p = periodicity(sp.sympify(expr), sp.sympify(var))
    except SYMBOLIC_QUERY_ERRORS:
        return None
    if p in (None, S.Zero, sp.oo, -sp.oo):
        return None
    p = _canonical_period(p)
    return p if p.is_positive is True else None


def common_period_values(*periods):
    ps = [_canonical_period(p) for p in periods]
    if not ps:
        return None
    p0 = ps[0]
    ratios = []
    for p in ps:
        r = sp.simplify(p / p0)
        if r.is_Rational is not True:
            return None
        ratios.append(sp.Rational(r))
    nums = [abs(int(r.p)) for r in ratios]
    dens = [int(r.q) for r in ratios]
    num = 1
    for n in nums:
        num = lcm(num, n)
    den = dens[0]
    for d in dens[1:]:
        den = gcd(den, d)
    return sp.simplify(p0 * sp.Rational(num, den))


def common_positive_period(expressions, var):
    ps = []
    for e in expressions:
        p = positive_period(e, var)
        if p is None:
            return None
        ps.append(p)
    return common_period_values(*ps)
