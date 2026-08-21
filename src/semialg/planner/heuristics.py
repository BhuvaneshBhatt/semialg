from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass
from itertools import permutations

import sympy as sp

from ..formula import Formula, equational_constraints, formula_polynomials
from .features import ProblemFeatures


@dataclass(frozen=True)
class OrderScore:
    order: tuple[sp.Symbol, ...]
    score: int
    reason: str = ""
    projection_poly_count: int | None = None
    projection_sotd: int | None = None
    estimated_lifting_roots: int | None = None
    estimated_cell_count: int | None = None
    estimated_alg_degree: int | None = None
    coefficient_height_bits: int | None = None
    pilot_lifting_roots: int | None = None
    pilot_cell_count: int | None = None


def _poly(poly: sp.Expr, order: Sequence[sp.Symbol]) -> sp.Poly | None:
    try:
        return sp.Poly(sp.expand(poly), *order)
    except (sp.PolynomialError, TypeError, ValueError):
        return None


def _degree_in(poly: sp.Expr, sym: sp.Symbol) -> int:
    try:
        return int(sp.Poly(sp.expand(poly), sym).degree())
    except (sp.PolynomialError, TypeError, ValueError):
        return 0


def _total_degree(poly: sp.Expr, order: Sequence[sp.Symbol]) -> int:
    pobj = _poly(poly, order)
    return int(pobj.total_degree()) if pobj is not None else 0


def brown_variable_order(
    polys: Sequence[sp.Expr], variables: Sequence[sp.Symbol]
) -> tuple[sp.Symbol, ...]:
    """Brown-style CAD variable ordering heuristic.

    Brown's heuristic eliminates variables using low maximum degree, then low
    total degree in terms containing that variable, then low occurrence count.
    The returned order is the corresponding CAD variable order.
    """

    vars_tuple = tuple(variables)

    def key(sym: sp.Symbol):
        max_degree = max((_degree_in(poly, sym) for poly in polys), default=0)
        degree_sum = sum(_degree_in(poly, sym) for poly in polys)
        occurrence = sum(1 for poly in polys if sym in poly.free_symbols)
        return (max_degree, degree_sum, occurrence, sym.name)

    return tuple(sorted(vars_tuple, key=key))


def sotd_score(order: Sequence[sp.Symbol], polys: Sequence[sp.Expr]) -> int:
    """Approximate sum-of-total-degrees score for an order.

    This is intentionally cheap: it scores the input family under the proposed
    order without constructing a full projection tower. It is suitable for
    ranking candidate orders before CAD construction.
    """

    return sum(_total_degree(poly, order) for poly in polys)


def ndrr_score(order: Sequence[sp.Symbol], polys: Sequence[sp.Expr]) -> int:
    """Cheap NDRR proxy: count distinct univariate real roots after projection to first variable.

    The full NDRR heuristic requires projection. This proxy keeps decomposition selection
    lightweight by looking at input polynomials that are univariate in the first
    variable under the proposed order.
    """

    if not order:
        return 0
    first = order[0]
    count = 0
    seen: set[str] = set()
    for poly in polys:
        if poly.free_symbols and poly.free_symbols <= {first}:
            try:
                roots = sp.Poly(poly, first).real_roots()
            except (sp.PolynomialError, TypeError, ValueError, NotImplementedError):
                roots = ()
            for root in roots:
                key = sp.sstr(root)
                if key not in seen:
                    seen.add(key)
                    count += 1
    return count


def _projection_complexity(
    order: Sequence[sp.Symbol], polys: Sequence[sp.Expr]
) -> tuple[int, int] | None:
    """Return exact Collins projection size/SOTD for a candidate order.

    Used only for a small shortlist: projection cost is a much better predictor
    than input-only proxies, and the tower is cached for later CAD use.
    """
    if not order or len(order) > 4:
        return None
    try:
        from ..cad.projection.collins import build_collins_proj_set

        tower = build_collins_proj_set(polys, order)
    except (sp.PolynomialError, TypeError, ValueError, NotImplementedError, ArithmeticError):
        return None
    count = 0
    total_degree = 0
    for level in tower.levels:
        count += len(level.polynomials)
        for poly in level.polynomials:
            try:
                total_degree += int(poly.total_degree())
            except (TypeError, ValueError, AttributeError):
                total_degree += 1
    return count, total_degree


def _probe_root_count(polys: Sequence[sp.Poly], order: Sequence[sp.Symbol], level: int) -> int:
    """Estimate distinct fiber roots at a few cheap exact rational probes.

    This is intentionally only an ordering heuristic.  It never participates in
    CAD correctness and is restricted to the shortlist already paying projection
    cost.
    """
    if level < 1 or level > len(order):
        return 0
    var = order[level - 1]
    lower = tuple(order[: level - 1])
    probes = (sp.Integer(0), sp.Integer(1), sp.Integer(-1))
    best = 0
    for probe in probes:
        subs = {sym: probe for sym in lower}
        roots: list[object] = []
        for poly in polys:
            expr = sp.expand(poly.as_expr().subs(subs))
            if expr == 0:
                continue
            try:
                univar = sp.Poly(expr, var, domain="EX")
            except (sp.PolynomialError, TypeError, ValueError):
                continue
            if univar.degree() <= 0:
                continue
            try:
                from ..algebraic.roots import isolate_real_roots

                roots.extend(isolate_real_roots(univar))
            except (sp.PolynomialError, TypeError, ValueError, NotImplementedError):
                # Degree remains a safe cheap proxy if exact probing fails.
                best = max(best, int(univar.degree()))
        if roots:
            from ..algebraic.comparison import sort_samples

            try:
                best = max(best, len(sort_samples(tuple(roots))))
            except (TypeError, ValueError):
                best = max(best, len(roots))
    return best


def _lifting_complexity(
    order: Sequence[sp.Symbol], polys: Sequence[sp.Expr]
) -> tuple[int, int] | None:
    """Estimate lifting roots/cells from a cached projection tower and exact probes."""
    if not order or len(order) > 4:
        return None
    try:
        from ..cad.projection.collins import build_collins_proj_set

        tower = build_collins_proj_set(polys, order)
    except (sp.PolynomialError, TypeError, ValueError, NotImplementedError, ArithmeticError):
        return None
    total_roots = 0
    estimated_cells = 1
    for level in range(1, len(order) + 1):
        level_polys = tower.level(level).polynomials
        root_count = _probe_root_count(level_polys, order, level)
        if root_count == 0:
            # Avoid incorrectly predicting a trivial lift merely because the
            # rational probes hit a degenerate point.
            degree_upper = sum(max(0, int(poly.degree(order[level - 1]))) for poly in level_polys)
            root_count = min(degree_upper, 8)
        total_roots += root_count
        estimated_cells *= max(1, 2 * root_count + 1)
        estimated_cells = min(estimated_cells, 10**7)
    return total_roots, estimated_cells


def _integer_height_bits(value: sp.Expr) -> int:
    value = sp.sympify(value)
    if value.is_Rational:
        rat = sp.Rational(value)
        return max(abs(int(rat.p)).bit_length(), abs(int(rat.q)).bit_length())
    if value.is_Integer:
        return abs(int(value)).bit_length()
    return max(1, len(sp.srepr(value)).bit_length())


def _projection_arithmetic_complexity(
    order: Sequence[sp.Symbol], polys: Sequence[sp.Expr]
) -> tuple[int, int] | None:
    """Estimate algebraic degree and coefficient-height growth in a tower."""
    if not order or len(order) > 4:
        return None
    try:
        from ..cad.projection.collins import build_collins_proj_set

        tower = build_collins_proj_set(polys, order)
    except (sp.PolynomialError, TypeError, ValueError, NotImplementedError):
        return None
    max_degree = 1
    max_height = 1
    for level in tower.levels:
        for poly in level.polynomials:
            try:
                max_degree = max(max_degree, int(poly.total_degree()))
                for coeff in poly.coeffs():
                    max_height = max(max_height, _integer_height_bits(coeff))
            except (TypeError, ValueError, AttributeError):
                continue
    return max_degree, max_height


def _sample_algebraic_degree(sample: object) -> int:
    try:
        from ..algebraic.samples import AlgebraicRoot

        if isinstance(sample, AlgebraicRoot):
            return max(1, int(sample.polynomial.degree()))
    except (TypeError, ValueError, AttributeError):
        pass
    return 1


def _pilot_lifting_complexity(
    order: Sequence[sp.Symbol], polys: Sequence[sp.Expr], *, max_parents: int = 3
) -> tuple[int, int, int] | None:
    """Perform a bounded real lift on representative stacks for order scoring.

    Unlike the rational-probe estimate, this uses the actual CAD lifting/root
    isolation primitives.  Only a few parents per level are followed, so it is
    a cost predictor rather than a correctness-producing decomposition.
    """
    if not order or len(order) > 4:
        return None
    try:
        from ..algebraic.comparison import sort_samples
        from ..algebraic.roots import isolate_real_roots
        from ..cad.decomposition import _build_stack, _stack_roots_over_point
        from ..cad.projection.collins import build_collins_proj_set

        tower = build_collins_proj_set(polys, order)
        roots = []
        for poly in tower.level(1).polynomials:
            roots.extend(isolate_real_roots(poly))
        roots = list(sort_samples(tuple(roots)))
        cells = list(_build_stack(None, roots, 1, tower))
        root_total = len(roots)
        cell_total = len(cells)
        algebraic_degree = max((_sample_algebraic_degree(root) for root in roots), default=1)

        def choose_parents(candidates):
            # Prefer rational sample paths for the bounded pilot.  They still
            # exercise the real lifting/root machinery but avoid turning a
            # cheap cost probe into a full algebraic-coefficient root-isolation
            # problem.  If no rational path exists, retain representative
            # candidates and let the exact pilot decline conservatively.
            from ..algebraic.samples import RationalSample

            rational = [
                cell
                for cell in candidates
                if all(isinstance(sample, RationalSample) for sample in cell.sample)
            ]
            pool = rational or list(candidates)
            if len(pool) <= max_parents:
                return pool
            indexes = {0, len(pool) // 2, len(pool) - 1}
            return [pool[i] for i in sorted(indexes)[:max_parents]]

        parents = choose_parents(cells)
        for level in range(2, len(order) + 1):
            next_parents = []
            for parent in parents:
                stack_roots = _stack_roots_over_point(
                    tower.level(level).polynomials,
                    order,
                    parent.sample,
                    order[level - 1],
                )
                root_total += len(stack_roots)
                algebraic_degree = max(
                    algebraic_degree,
                    max((_sample_algebraic_degree(root) for root in stack_roots), default=1),
                )
                stack = list(_build_stack(parent, stack_roots, level, tower))
                cell_total += len(stack)
                next_parents.extend(stack)
            if not next_parents:
                break
            parents = choose_parents(next_parents)
        return root_total, cell_total, algebraic_degree
    except (sp.PolynomialError, TypeError, ValueError, NotImplementedError, ArithmeticError):
        return None


def score_variable_order(
    order: Sequence[sp.Symbol],
    polys: Sequence[sp.Expr],
    *,
    ec_exprs: Sequence[sp.Expr] = (),
    include_projection: bool = False,
    include_lifting: bool = False,
    include_pilot: bool = False,
) -> OrderScore:
    """Score one CAD variable order using symbolic and measured cost proxies.

    Cheap Brown, SOTD, NDRR, and equational-constraint terms are always used.
    Projection, lifting, arithmetic-growth, and bounded pilot measurements are
    optional because they are progressively more expensive to obtain.
    """

    order = tuple(order)
    position = {sym: i for i, sym in enumerate(order)}
    score = 0
    brown_rank = {sym: idx for idx, sym in enumerate(brown_variable_order(polys, order))}
    score += 40 * sum(abs(position[sym] - brown_rank[sym]) for sym in order)
    score += 8 * sotd_score(order, polys)
    score += 20 * ndrr_score(order, polys)
    for poly in polys:
        syms = [sym for sym in order if sym in poly.free_symbols]
        if not syms:
            continue
        width = position[syms[-1]] - position[syms[0]]
        score += width + len(syms)

    # Reduced projection can exploit an EC when its main variable is high in
    # the CAD order (high variables are projected first). Prefer orders that
    # place a low-degree designated EC at such a level.
    for ec in ec_exprs:
        ec_syms = [sym for sym in order if sym in ec.free_symbols]
        if not ec_syms:
            continue
        main_pos = max(position[sym] for sym in ec_syms)
        try:
            degree = int(sp.Poly(sp.expand(ec), *order).total_degree())
        except (sp.PolynomialError, TypeError, ValueError):
            degree = 1
        score += 15 * (len(order) - 1 - main_pos) + 2 * degree

    projection_count = None
    projection_degree = None
    if include_projection:
        complexity = _projection_complexity(order, polys)
        if complexity is not None:
            projection_count, projection_degree = complexity
            score += 75 * projection_count + 5 * projection_degree
    lifting_roots = None
    cell_count = None
    algebraic_degree = None
    coefficient_height = None
    pilot_roots = None
    pilot_cells = None
    if include_lifting:
        lifting = _lifting_complexity(order, polys)
        if lifting is not None:
            lifting_roots, cell_count = lifting
            # Lifting dominates once stacks branch.  Log-like integer scaling
            # keeps this estimate from overwhelming exact projection metrics.
            score += 35 * lifting_roots + 12 * int(max(0, cell_count).bit_length())
    arithmetic = (
        _projection_arithmetic_complexity(order, polys)
        if (include_projection or include_lifting)
        else None
    )
    if arithmetic is not None:
        algebraic_degree, coefficient_height = arithmetic
        score += 10 * algebraic_degree + 2 * coefficient_height
    if include_pilot:
        pilot = _pilot_lifting_complexity(order, polys)
        if pilot is not None:
            pilot_roots, pilot_cells, pilot_degree = pilot
            # Pilot lifting refines, rather than duplicates, the cheap lifting
            # estimate so piloted and unpiloted orders remain comparable.
            baseline_roots = lifting_roots or 0
            baseline_cell_bits = int(max(1, cell_count or 1).bit_length())
            pilot_cell_bits = int(max(1, pilot_cells).bit_length())
            baseline_degree = algebraic_degree or 1
            score += 20 * (pilot_roots - baseline_roots)
            score += 8 * (pilot_cell_bits - baseline_cell_bits)
            score += 12 * (pilot_degree - baseline_degree)
            algebraic_degree = max(baseline_degree, pilot_degree)
    return OrderScore(
        order=order,
        score=score,
        reason="Brown/input-SOTD/NDRR/EC plus shortlisted projection and lifting/root-count estimates; smaller is better",
        projection_poly_count=projection_count,
        projection_sotd=projection_degree,
        estimated_lifting_roots=lifting_roots,
        estimated_cell_count=cell_count,
        estimated_alg_degree=algebraic_degree,
        coefficient_height_bits=coefficient_height,
        pilot_lifting_roots=pilot_roots,
        pilot_cell_count=pilot_cells,
    )


def candidate_variable_orders(
    features: ProblemFeatures,
    polys: Sequence[sp.Expr],
    *,
    equational_constraints: Sequence[sp.Expr] = (),
    limit: int = 12,
) -> tuple[OrderScore, ...]:
    """Return the best candidate CAD orders under a staged cost model.

    The planner generates structural candidates, ranks them cheaply, measures
    projection and lifting cost for a short list, and runs a bounded pilot lift
    only for the strongest remaining candidates.
    """

    vars_ = tuple(features.variables)
    if len(vars_) <= 1:
        return (OrderScore(order=vars_, score=0, reason="single variable"),)
    ecs = tuple(sp.expand(ec) for ec in equational_constraints)
    sorted_vars = tuple(sorted(vars_, key=lambda s: s.name))
    brown = brown_variable_order(polys, vars_)
    candidates: set[tuple[sp.Symbol, ...]] = {
        tuple(vars_),
        sorted_vars,
        tuple(reversed(sorted_vars)),
        brown,
        tuple(reversed(brown)),
        tuple(
            sorted(
                vars_, key=lambda sym: (sum(sym in poly.free_symbols for poly in polys), sym.name)
            )
        ),
        tuple(
            sorted(
                vars_, key=lambda sym: (-sum(sym in poly.free_symbols for poly in polys), sym.name)
            )
        ),
    }
    # EC-aware seeds: put variables serving as main variables of simple ECs
    # toward the high end of the CAD order.
    if ecs:
        ec_occurrence = {sym: sum(sym in ec.free_symbols for ec in ecs) for sym in vars_}
        candidates.add(tuple(sorted(vars_, key=lambda sym: (ec_occurrence[sym], sym.name))))
        candidates.add(tuple(sorted(vars_, key=lambda sym: (-ec_occurrence[sym], sym.name))))
    if len(vars_) <= 5:
        all_orders = list(permutations(vars_))
        prelim = sorted(
            (score_variable_order(order, polys, ec_exprs=ecs) for order in all_orders),
            key=lambda item: item.score,
        )
        for item in prelim[: max(limit, 8)]:
            candidates.add(item.order)

    prelim_scores = sorted(
        (score_variable_order(cand, polys, ec_exprs=ecs) for cand in candidates),
        key=lambda item: (item.score, tuple(sym.name for sym in item.order)),
    )
    # Exact projection/lifting scoring is deliberately limited to a shortlist
    # and at most four variables.  Pilot lifting is a second stage over the two
    # best projection-scored orders, so its measured adjustment is comparable
    # with the remaining estimated candidates.
    shortlist = (
        {item.order for item in prelim_scores[: min(6, len(prelim_scores))]}
        if len(vars_) <= 4
        else set()
    )
    baseline_scores = [
        score_variable_order(
            item.order,
            polys,
            ec_exprs=ecs,
            include_projection=item.order in shortlist,
            include_lifting=item.order in shortlist,
        )
        for item in prelim_scores
    ]
    baseline_scores.sort(key=lambda item: (item.score, tuple(sym.name for sym in item.order)))
    pilot_shortlist = (
        {item.order for item in baseline_scores[: min(2, len(baseline_scores))]}
        if len(vars_) <= 4
        else set()
    )
    scores = [
        score_variable_order(
            item.order,
            polys,
            ec_exprs=ecs,
            include_projection=item.order in shortlist,
            include_lifting=item.order in shortlist,
            include_pilot=item.order in pilot_shortlist,
        )
        if item.order in pilot_shortlist
        else item
        for item in baseline_scores
    ]
    scores.sort(key=lambda item: (item.score, tuple(sym.name for sym in item.order)))
    return tuple(scores[:limit])


def choose_best_variable_order(
    features: ProblemFeatures,
    polys: Sequence[sp.Expr],
    *,
    equational_constraints: Sequence[sp.Expr] = (),
) -> tuple[sp.Symbol, ...]:
    return candidate_variable_orders(
        features, polys, equational_constraints=equational_constraints
    )[0].order


def choose_formula_variable_order(
    formula: Formula, features: ProblemFeatures
) -> tuple[sp.Symbol, ...]:
    return choose_best_variable_order(
        features,
        tuple(formula_polynomials(formula)),
        equational_constraints=tuple(equational_constraints(formula)),
    )


__all__ = [
    "OrderScore",
    "brown_variable_order",
    "candidate_variable_orders",
    "choose_formula_variable_order",
    "choose_best_variable_order",
    "ndrr_score",
    "score_variable_order",
    "sotd_score",
]
