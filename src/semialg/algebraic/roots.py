from __future__ import annotations

from collections.abc import Iterable
from functools import cmp_to_key

import sympy as sp

from ..exact_arithmetic import compare_exact_reals
from .cache import CACHE, poly_key
from .intervals import RationalInterval
from .samples import AlgebraicRoot


def _as_univariate_poly(poly: sp.Poly | sp.Expr, var: sp.Symbol | None = None) -> sp.Poly:
    if isinstance(poly, sp.Poly):
        if poly.is_multivariate:
            if var is None:
                raise ValueError("a variable is required for multivariate expressions")
            return sp.Poly(poly.as_expr(), var, domain="EX")
        return poly
    if var is None:
        symbols = sorted(poly.free_symbols, key=lambda s: s.name)
        if len(symbols) != 1:
            raise ValueError("a variable is required unless the expression is univariate")
        var = symbols[0]
    return sp.Poly(poly, var, domain="EX")


def _root_multiplicities(poly: sp.Poly) -> dict[sp.Expr, int]:
    out: dict[sp.Expr, int] = {}
    try:
        _, factors = sp.factor_list(poly.as_expr())
    except (NotImplementedError, sp.PolynomialError, ValueError, TypeError):
        return out
    var = poly.gens[0]
    for factor_expr, mult in factors:
        factor_poly = sp.Poly(
            factor_expr, var, domain=poly.domain if poly.domain != sp.EX else "EX"
        )
        try:
            roots = sp.real_roots(factor_poly.as_expr())
        except (NotImplementedError, sp.PolynomialError, ValueError):
            roots = []
        for root in roots:
            out[root] = mult
    return out


def rational_intv_around(root: sp.Expr) -> RationalInterval:
    """Return a certified rational isolating interval for a real algebraic root.

    The enclosure is derived from exact polynomial isolation.  Decimal
    approximations are never used to manufacture certificate metadata.
    """

    root = sp.sympify(root)
    if root.is_Rational:
        value = sp.Rational(root)
        return RationalInterval(value, value)
    # ``CRootOf`` already carries a certified rational isolating interval.
    # Recomputing a minimal polynomial and isolating all of its roots is much
    # more expensive, especially for algebraic-coefficient fibers whose
    # elimination polynomial has higher degree than the original fiber.
    if isinstance(root, sp.polys.rootoftools.ComplexRootOf) and root.is_real is True:
        interval = root._get_interval()
        left, right = interval.as_tuple() if hasattr(interval, "as_tuple") else tuple(interval)
        return RationalInterval(sp.Rational(left), sp.Rational(right))
    t = sp.Dummy("_alpha", real=True)
    try:
        minimal = sp.Poly(sp.minpoly(root, t), t, domain=sp.QQ)
        intervals = minimal.intervals(eps=sp.Rational(1, 10**30))
    except (NotImplementedError, sp.PolynomialError, ValueError, TypeError) as exc:
        raise ValueError(
            "could not construct an exact isolating interval for algebraic root"
        ) from exc
    for bounds, _multiplicity in intervals:
        left, right = map(sp.Rational, bounds)
        if compare_exact_reals(root, left) >= 0 and compare_exact_reals(root, right) <= 0:
            return RationalInterval(left, right)
    raise ValueError("minimal-polynomial intervals did not isolate the requested algebraic root")


def _poly_intervals(poly: sp.Poly) -> list[tuple[RationalInterval, int]]:
    try:
        exact_poly = sp.Poly(poly.as_expr(), poly.gens[0], extension=True)
        raw = exact_poly.intervals(eps=sp.Rational(1, 10**30))
    except (
        NotImplementedError,
        sp.PolynomialError,
        sp.polys.polyerrors.DomainError,
        ValueError,
        TypeError,
    ):
        return []
    intervals: list[tuple[RationalInterval, int]] = []
    for bounds, mult in raw:
        left, right = bounds
        intervals.append((RationalInterval(sp.Rational(left), sp.Rational(right)), int(mult)))
    return intervals


def _count_roots_in_interval(
    poly: sp.Poly,
    left: sp.Rational,
    right: sp.Rational,
    *,
    exclude_left: bool = False,
    exclude_right: bool = False,
) -> int:
    """Count roots exactly while optionally excluding root endpoints."""

    count = int(poly.count_roots(left, right))
    if exclude_left and sp.simplify(poly.eval(left)) == 0:
        count -= 1
    if exclude_right and sp.simplify(poly.eval(right)) == 0:
        count -= 1
    return count


def _trim_excluded_endpoint(
    poly: sp.Poly,
    left: sp.Rational,
    right: sp.Rational,
    *,
    exclude_left: bool,
    exclude_right: bool,
) -> tuple[sp.Rational, sp.Rational]:
    """Move an excluded root endpoint away from a one-root interval."""

    for _ in range(96):
        if not exclude_left and not exclude_right:
            break
        mid = sp.Rational(left + right, 2)
        mid_is_root = sp.simplify(poly.eval(mid)) == 0
        if mid_is_root:
            # The interval's unique included root is the midpoint unless that
            # root is excluded by construction, which cannot occur for a
            # strict interior midpoint.
            return mid, mid
        if exclude_left:
            right_count = _count_roots_in_interval(poly, mid, right, exclude_right=exclude_right)
            if right_count == 1:
                left = mid
                exclude_left = False
                continue
            right = mid
        elif exclude_right:
            left_count = _count_roots_in_interval(poly, left, mid, exclude_left=exclude_left)
            if left_count == 1:
                right = mid
                exclude_right = False
                continue
            left = mid
    return left, right


def _isolate_algebraic_coeff_roots(poly: sp.Poly) -> list[tuple[RationalInterval, int]]:
    """Isolate real roots over an exact algebraic coefficient field.

    SymPy can count real roots over ``QQ<alpha>`` efficiently even where its
    general ``RootOf`` association machinery is expensive.  This routine uses
    those exact counts to construct pairwise ordered rational isolating
    intervals.  It never uses a floating-point approximation to choose a root.
    """

    exact_poly = sp.Poly(poly.as_expr(), poly.gens[0], extension=True)
    try:
        factors = exact_poly.sqf_list()[1]
    except (NotImplementedError, sp.PolynomialError, ValueError, TypeError):
        factors = [(exact_poly, 1)]
    isolated: list[tuple[RationalInterval, int]] = []
    for factor, multiplicity in factors:
        total = int(factor.count_roots(-sp.oo, sp.oo))
        if total == 0:
            continue
        bound = sp.Integer(1)
        for _ in range(128):
            if (
                sp.simplify(factor.eval(-bound)) != 0
                and sp.simplify(factor.eval(bound)) != 0
                and int(factor.count_roots(-bound, bound)) == total
            ):
                break
            bound *= 2
        else:
            raise NotImplementedError("could not bound all exact real roots")

        def recurse(
            left: sp.Rational,
            right: sp.Rational,
            count: int,
            *,
            exclude_left: bool = False,
            exclude_right: bool = False,
            _factor: sp.Poly = factor,
            _multiplicity: int = multiplicity,
        ) -> None:
            if count <= 0:
                return
            if count == 1:
                left, right = _trim_excluded_endpoint(
                    _factor,
                    left,
                    right,
                    exclude_left=exclude_left,
                    exclude_right=exclude_right,
                )
                isolated.append((RationalInterval(left, right), int(_multiplicity)))
                return
            mid = sp.Rational(left + right, 2)
            mid_is_root = sp.simplify(_factor.eval(mid)) == 0
            if mid_is_root:
                isolated.append((RationalInterval(mid, mid), int(_multiplicity)))
            left_count = _count_roots_in_interval(
                _factor,
                left,
                mid,
                exclude_left=exclude_left,
                exclude_right=mid_is_root,
            )
            right_count = _count_roots_in_interval(
                _factor,
                mid,
                right,
                exclude_left=mid_is_root,
                exclude_right=exclude_right,
            )
            recurse(
                left,
                mid,
                left_count,
                exclude_left=exclude_left,
                exclude_right=mid_is_root,
            )
            recurse(
                mid,
                right,
                right_count,
                exclude_left=mid_is_root,
                exclude_right=exclude_right,
            )

        recurse(-sp.Rational(bound), sp.Rational(bound), total)
    isolated.sort(key=lambda item: (item[0].left, item[0].right))
    return isolated


def _dedupe_sorted_roots(roots: Iterable[sp.Expr]) -> list[sp.Expr]:
    """Sort and deduplicate exact real algebraic roots without numerics."""

    ordered = sorted(roots, key=cmp_to_key(compare_exact_reals))
    unique: list[sp.Expr] = []
    for root in ordered:
        if not unique or compare_exact_reals(unique[-1], root) != 0:
            unique.append(root)
    return unique


def isolate_real_roots(
    poly: sp.Poly | sp.Expr, var: sp.Symbol | None = None
) -> tuple[AlgebraicRoot, ...]:
    """Return explicit algebraic samples for the real roots of a univariate polynomial."""

    univar = _as_univariate_poly(poly, var)
    if univar.degree() <= 0:
        return ()
    key = (poly_key(univar), str(univar.gens[0]))
    CACHE.stats.calls += 1
    cached = CACHE.roots.get(key)
    if cached is not None:
        CACHE.stats.cache_hits += 1
        return cached  # type: ignore[return-value]
    CACHE.stats.cache_misses += 1

    expr = sp.expand(univar.as_expr())
    roots_preordered = False
    try:
        raw_roots = sp.real_roots(expr)
        roots_preordered = True
    except (NotImplementedError, sp.PolynomialError, ValueError):
        # Linear and quadratic algebraic-coefficient fibers have simple exact
        # radical representations that remain usable as coefficients at later
        # CAD levels.  Prefer those to opaque interval-root selectors.
        raw_roots = []
        solved_simple = False
        if univar.degree() <= 2:
            try:
                roots_with_mult = sp.roots(expr, univar.gens[0])
                solved_simple = (
                    sum(int(mult) for mult in roots_with_mult.values()) == univar.degree()
                )
                if solved_simple:
                    for root, multiplicity in roots_with_mult.items():
                        reality = root.is_real
                        if reality is None:
                            try:
                                reality = compare_exact_reals(sp.im(root), 0) == 0
                            except ValueError:
                                solved_simple = False
                                break
                        if reality is True:
                            raw_roots.extend([sp.simplify(root)] * int(multiplicity))
            except (NotImplementedError, sp.PolynomialError, ValueError, TypeError):
                solved_simple = False

        if not solved_simple:
            try:
                exact_intervals = _isolate_algebraic_coeff_roots(univar)
            except (
                NotImplementedError,
                sp.PolynomialError,
                sp.polys.polyerrors.DomainError,
                ValueError,
                TypeError,
            ):
                exact_intervals = []
            if exact_intervals:
                samples = tuple(
                    AlgebraicRoot(
                        polynomial=univar,
                        interval=interval,
                        root_index=index,
                        multiplicity=multiplicity,
                        root_expr=None,
                    )
                    for index, (interval, multiplicity) in enumerate(exact_intervals)
                )
                CACHE.roots.put(key, samples)
                return samples

            # Final exact fallback for low-degree domains where root counting is
            # not available.  Radical roots are accepted only if their reality
            # and ordering can themselves be certified exactly.
            solved_low_degree = False
            if univar.degree() <= 4:
                try:
                    roots_with_mult = sp.roots(expr, univar.gens[0])
                    solved_low_degree = (
                        sum(int(mult) for mult in roots_with_mult.values()) == univar.degree()
                    )
                    if solved_low_degree:
                        for root, multiplicity in roots_with_mult.items():
                            reality = root.is_real
                            if reality is None:
                                try:
                                    reality = compare_exact_reals(sp.im(root), 0) == 0
                                except ValueError:
                                    solved_low_degree = False
                                    break
                            if reality is True:
                                raw_roots.extend([sp.simplify(root)] * int(multiplicity))
                except (
                    NotImplementedError,
                    sp.PolynomialError,
                    ValueError,
                    TypeError,
                ):
                    solved_low_degree = False
            if not solved_low_degree:
                raise NotImplementedError(
                    "exact real-root isolation is unavailable for this polynomial"
                ) from None
    if roots_preordered:
        roots = []
        for root in raw_roots:
            if not roots or root != roots[-1]:
                roots.append(root)
    else:
        roots = _dedupe_sorted_roots(raw_roots)
    intervals = _poly_intervals(univar)
    multiplicities = _root_multiplicities(univar)
    samples: list[AlgebraicRoot] = []
    for index, root in enumerate(roots):
        mult = next(
            (m for candidate, m in multiplicities.items() if sp.simplify(candidate - root) == 0), 1
        )
        interval = rational_intv_around(root)
        if index < len(intervals):
            interval, interval_mult = intervals[index]
            mult = interval_mult or mult
        samples.append(
            AlgebraicRoot(
                polynomial=univar,
                interval=interval,
                root_index=index,
                multiplicity=mult,
                root_expr=root,
            )
        )
    result = tuple(samples)
    CACHE.roots.put(key, result)
    return result


def refine_isol_intv(root: AlgebraicRoot, *, steps: int = 4) -> AlgebraicRoot:
    """Return the same root with a narrower certified rational interval.

    Algebraic-coefficient fibers may use an opaque ordered-root expression that
    SymPy cannot convert to a native algebraic number.  In that case refinement
    uses exact Sturm/root counts for the defining polynomial instead of trying
    to compare the opaque expression numerically or symbolically.
    """

    if root.interval.is_point():
        return root
    left, right = root.interval.left, root.interval.right
    if root.root_expr is None:
        try:
            poly = sp.Poly(root.polynomial.as_expr(), root.variable, extension=True)
            for _ in range(max(steps, 0)):
                CACHE.stats.refinements += 1
                mid = sp.Rational(left + right, 2)
                if sp.simplify(poly.eval(mid)) == 0:
                    left = right = mid
                    break
                left_count = int(poly.count_roots(left, mid))
                if left_count >= 1:
                    right = mid
                else:
                    left = mid
            return AlgebraicRoot(
                root.polynomial,
                RationalInterval(left, right),
                root.root_index,
                root.multiplicity,
                root.root_expr,
            )
        except (
            NotImplementedError,
            sp.PolynomialError,
            sp.polys.polyerrors.DomainError,
            ValueError,
            TypeError,
        ):
            pass
    expr = root.as_expr()
    for _ in range(max(steps, 0)):
        CACHE.stats.refinements += 1
        mid = sp.Rational(left + right, 2)
        if compare_exact_reals(expr, mid) <= 0:
            right = mid
        else:
            left = mid
    return AlgebraicRoot(
        root.polynomial,
        RationalInterval(left, right),
        root.root_index,
        root.multiplicity,
        root.root_expr,
    )


def root_multiplicity(root: AlgebraicRoot) -> int:
    return root.multiplicity
