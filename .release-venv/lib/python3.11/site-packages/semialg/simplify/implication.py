from __future__ import annotations

from collections.abc import Iterable
from dataclasses import dataclass
from functools import lru_cache

import sympy as sp
from sympy.core.relational import Relational
from sympy.logic.boolalg import And as SymAnd

from ..formulas.boolean import bounded_dnf_branches, make_and
from ..structural_keys import symbol_identity_key
from .boolean import simplify_boolean


@dataclass
class ImplicationMinimizationStats:
    """Lightweight counters for profiling implication-based simplification."""

    minimize_calls: int = 0
    input_branches: int = 0
    output_branches: int = 0
    conjunction_checks: int = 0
    disjunction_checks: int = 0
    removed_atoms: int = 0
    removed_branches: int = 0
    cad_requests: int = 0
    cad_cache_misses: int = 0
    guard_rejections: int = 0
    shared_cad_builds: int = 0
    truth_evaluations: int = 0
    pairwise_fallbacks: int = 0
    generated_implicants: int = 0
    cover_candidates: int = 0
    cover_selected: int = 0

    @property
    def cad_cache_hits(self) -> int:
        # Shared-CAD requests bypass the pairwise implication LRU, so they
        # are excluded from these cache-hit statistics.
        return max(0, self.cad_requests - self.shared_cad_builds - self.cad_cache_misses)


_STATS = ImplicationMinimizationStats()


def clear_implication_minimization_stats() -> None:
    """Reset implication-minimization profiling counters and its CAD cache."""

    fresh = ImplicationMinimizationStats()
    _STATS.__dict__.update(fresh.__dict__)
    _is_unsatisfiable_cached.cache_clear()


def implication_minimization_stats() -> ImplicationMinimizationStats:
    """Return a snapshot of implication-minimization profiling counters."""

    return ImplicationMinimizationStats(
        minimize_calls=_STATS.minimize_calls,
        input_branches=_STATS.input_branches,
        output_branches=_STATS.output_branches,
        conjunction_checks=_STATS.conjunction_checks,
        disjunction_checks=_STATS.disjunction_checks,
        removed_atoms=_STATS.removed_atoms,
        removed_branches=_STATS.removed_branches,
        cad_requests=_STATS.cad_requests,
        cad_cache_misses=_STATS.cad_cache_misses,
        guard_rejections=_STATS.guard_rejections,
        shared_cad_builds=_STATS.shared_cad_builds,
        truth_evaluations=_STATS.truth_evaluations,
        pairwise_fallbacks=_STATS.pairwise_fallbacks,
        generated_implicants=_STATS.generated_implicants,
        cover_candidates=_STATS.cover_candidates,
        cover_selected=_STATS.cover_selected,
    )


def _dnf_branches(expr: sp.Expr, *, max_branches: int = 4096) -> tuple[tuple[sp.Expr, ...], ...]:
    """Return a DNF branch representation for a SymPy Boolean expression."""

    simplified = simplify_boolean(expr)
    expansion = bounded_dnf_branches(simplified, max_branches=max_branches)
    if not expansion.complete:
        return (_branch_atoms(simplified),)
    branches = []
    for branch in expansion.branches:
        expr_branch = [item for item in branch if item is not sp.true and item is not True]
        if any(item is sp.false or item is False for item in branch):
            continue
        branches.append(_branch_atoms(make_and(*expr_branch)))
    return tuple(branches)


def _branch_atoms(expr: sp.Expr) -> tuple[sp.Expr, ...]:
    if expr is sp.true:
        return tuple()
    if isinstance(expr, SymAnd):
        atoms = expr.args
    else:
        atoms = (expr,)
    return tuple(sorted(dict.fromkeys(atoms), key=sp.sstr))


def _branch_expr(branch: Iterable[sp.Expr]) -> sp.Expr:
    items = tuple(branch)
    if not items:
        return sp.true
    return sp.And(*items)


def _safe_not(expr: sp.Expr) -> sp.Expr:
    return sp.Not(expr)


def _complexity_ok(expr: sp.Expr, *, max_atoms: int, max_vars: int, max_degree: int) -> bool:
    atoms = tuple(expr.atoms(Relational))
    if len(atoms) > max_atoms:
        return False
    vars_ = sorted(expr.free_symbols, key=symbol_identity_key)
    if len(vars_) > max_vars:
        return False
    for atom in atoms:
        try:
            poly = sp.Poly(sp.expand(atom.lhs - atom.rhs), *vars_)
        except (sp.PolynomialError, TypeError, ValueError):
            return False
        if poly.total_degree() > max_degree:
            return False
    return True


@lru_cache(maxsize=512)
def _is_unsatisfiable_cached(expr: sp.Expr, variables: tuple[sp.Symbol, ...]) -> bool | None:
    """CAD unsatisfiability cache preserving exact SymPy symbol identity."""
    _STATS.cad_cache_misses += 1
    try:
        from ..cad_algorithms.decomposition import decomp_collins_complete
        from ..formula import formula_polynomials, parse_formula
        from ..qe.complete import evaluate_formula_on_cell

        formula = parse_formula(expr)
        polys = tuple(formula_polynomials(formula)) or (sp.Integer(1),)
        cad = decomp_collins_complete(polys, variables)
        level = len(variables)
        return not any(
            evaluate_formula_on_cell(formula, cell, variables)
            for cell in cad.cells_by_level.get(level, ())
        )
    except (
        sp.PolynomialError,
        TypeError,
        ValueError,
        ArithmeticError,
        NotImplementedError,
    ):
        return None


def is_unsatisfiable_by_cad(
    expr: sp.Expr, *, max_atoms: int = 8, max_vars: int = 3, max_degree: int = 4
) -> bool | None:
    """Conservatively decide small real-polynomial unsatisfiability by CAD.

    ``None`` means the implication minimizer should leave the formula alone.
    The size guard prevents the pretty-printer from accidentally launching an
    expensive CAD on every large output branch.
    """

    if expr is sp.false:
        return True
    if expr is sp.true:
        return False
    if not _complexity_ok(expr, max_atoms=max_atoms, max_vars=max_vars, max_degree=max_degree):
        _STATS.guard_rejections += 1
        return None
    variables = tuple(sorted(expr.free_symbols, key=symbol_identity_key))
    _STATS.cad_requests += 1
    return _is_unsatisfiable_cached(expr, variables)


def implies_by_cad(antecedent: sp.Expr, consequent: sp.Expr, **kwargs) -> bool | None:
    """Return whether ``antecedent => consequent`` for guarded small formulas."""

    return is_unsatisfiable_by_cad(sp.And(antecedent, _safe_not(consequent)), **kwargs)


def minimize_conj_by_impl(branch: tuple[sp.Expr, ...], **kwargs) -> tuple[sp.Expr, ...] | None:
    """Remove atoms implied by the other atoms in a conjunction."""

    atoms = list(branch)
    changed = False
    idx = 0
    while idx < len(atoms):
        atom = atoms[idx]
        rest = atoms[:idx] + atoms[idx + 1 :]
        if not rest:
            idx += 1
            continue
        _STATS.conjunction_checks += 1
        result = implies_by_cad(_branch_expr(rest), atom, **kwargs)
        if result is True:
            atoms.pop(idx)
            _STATS.removed_atoms += 1
            changed = True
            continue
        idx += 1
    return tuple(atoms) if changed else None


@dataclass
class _SharedCADTruthContext:
    """Truth-bitset evaluator over one sign-invariant CAD."""

    variables: tuple[sp.Symbol, ...]
    cells: tuple[object, ...]
    truth_cache: dict[sp.Expr, int]

    @property
    def full_mask(self) -> int:
        return (1 << len(self.cells)) - 1

    def truth_mask(self, expr: sp.Expr) -> int:
        expr = sp.sympify(expr)
        cached = self.truth_cache.get(expr)
        if cached is not None:
            return cached
        from ..formula import parse_formula
        from ..qe.complete import evaluate_formula_on_cell

        formula = parse_formula(expr)
        mask = 0
        for pos, cell in enumerate(self.cells):
            _STATS.truth_evaluations += 1
            if evaluate_formula_on_cell(formula, cell, self.variables):
                mask |= 1 << pos
        self.truth_cache[expr] = mask
        return mask

    def conjunction_mask(self, atoms: Iterable[sp.Expr]) -> int:
        mask = self.full_mask
        for atom in atoms:
            mask &= self.truth_mask(atom)
            if mask == 0:
                break
        return mask


def _shared_truth_context(
    expr: sp.Expr, *, max_atoms: int = 8, max_vars: int = 3, max_degree: int = 4
) -> _SharedCADTruthContext | None:
    """Build one guarded CAD containing every atom needed by a minimization pass."""

    if not _complexity_ok(expr, max_atoms=max_atoms, max_vars=max_vars, max_degree=max_degree):
        _STATS.guard_rejections += 1
        return None
    variables = tuple(sorted(expr.free_symbols, key=symbol_identity_key))
    if not variables:
        return None
    try:
        from ..cad_algorithms.decomposition import decomp_collins_complete
        from ..formula import formula_polynomials, parse_formula

        parsed = parse_formula(expr)
        polys = tuple(formula_polynomials(parsed)) or (sp.Integer(1),)
        _STATS.cad_requests += 1
        cad = decomp_collins_complete(polys, variables)
        _STATS.shared_cad_builds += 1
        cells = tuple(cad.cells_by_level.get(len(variables), ()))
        return _SharedCADTruthContext(variables, cells, {})
    except (
        sp.PolynomialError,
        TypeError,
        ValueError,
        ArithmeticError,
        NotImplementedError,
    ):
        return None


def _prime_implicants_from_branches(
    branches: list[tuple[sp.Expr, ...]],
    context: _SharedCADTruthContext,
    target_mask: int,
) -> list[tuple[tuple[sp.Expr, ...], int]]:
    """Generate exact CAD-valid generalized DNF cubes by literal deletion.

    Unlike implication-only minimization, a literal may be deleted even when it
    is not implied by the remaining literals, provided the resulting cube does
    not include any CAD cell outside the truth set of the whole formula.  This
    merges complementary/overlapping branches such as ``(A&B)|(A&~B)`` into
    ``A`` whenever the shared CAD proves the equivalence.
    """

    queue = {tuple(branch) for branch in branches}
    accepted: dict[tuple[sp.Expr, ...], int] = {}
    while queue:
        cube = queue.pop()
        mask = context.conjunction_mask(cube)
        if mask & ~target_mask:
            continue
        generalized = False
        for idx in range(len(cube)):
            _STATS.conjunction_checks += 1
            smaller = cube[:idx] + cube[idx + 1 :]
            smask = context.conjunction_mask(smaller)
            if smask and not (smask & ~target_mask):
                if smaller not in accepted:
                    queue.add(smaller)
                generalized = True
        if not generalized and mask:
            accepted[cube] = mask
    _STATS.generated_implicants += len(accepted)
    return list(accepted.items())


def _select_exact_cover(
    implicants: list[tuple[tuple[sp.Expr, ...], int]], target_mask: int
) -> list[tuple[sp.Expr, ...]]:
    """Choose a small exact implicant cover using essentials then greedy gain."""

    if not implicants:
        return []
    # Remove dominated cubes: when two cubes cover the same/more cells, prefer
    # the one with fewer literals.  The remaining set-cover step is exact in
    # truth semantics even though the size objective uses a deterministic greedy
    # tie-break instead of an exponential minimum-cover search.
    filtered: list[tuple[tuple[sp.Expr, ...], int]] = []
    for cube, mask in sorted(
        implicants, key=lambda item: (len(item[0]), tuple(map(sp.sstr, item[0])))
    ):
        dominated = False
        for other, other_mask in filtered:
            _STATS.disjunction_checks += 1
            if mask & ~other_mask == 0 and len(other) <= len(cube):
                dominated = True
                break
        if dominated:
            continue
        filtered.append((cube, mask))
    _STATS.cover_candidates += len(filtered)

    uncovered = target_mask
    selected: list[tuple[sp.Expr, ...]] = []
    while uncovered:
        best = max(
            filtered,
            key=lambda item: (
                (item[1] & uncovered).bit_count(),
                -len(item[0]),
                tuple(reversed(tuple(map(sp.sstr, item[0])))),
            ),
        )
        gain = best[1] & uncovered
        if not gain:
            break
        selected.append(best[0])
        uncovered &= ~best[1]
        filtered = [item for item in filtered if item[0] != best[0]]
    _STATS.cover_selected += len(selected)
    return selected


def _minimize_disj_with_shared_cad(
    branches: list[tuple[sp.Expr, ...]],
    context: _SharedCADTruthContext,
) -> sp.Expr:
    """Reconstruct a smaller DNF from exact shared-CAD truth sets."""

    target_mask = 0
    for branch in branches:
        target_mask |= context.conjunction_mask(branch)
    if target_mask == 0:
        _STATS.output_branches += 0
        return sp.false
    if target_mask == context.full_mask:
        _STATS.output_branches += 1
        return sp.true

    implicants = _prime_implicants_from_branches(branches, context, target_mask)
    selected = _select_exact_cover(implicants, target_mask)
    if not selected:
        selected = branches
    kept = [_branch_expr(branch) for branch in selected]
    _STATS.removed_branches += max(0, len(branches) - len(kept))
    _STATS.output_branches += len(kept)
    return simplify_boolean(sp.Or(*kept))


def _minimize_disj_pairwise(branches: list[tuple[sp.Expr, ...]], **kwargs) -> sp.Expr:
    """Minimize DNF branches with independent implication checks."""

    _STATS.pairwise_fallbacks += 1
    for pos, branch in enumerate(tuple(branches)):
        minimized = minimize_conj_by_impl(branch, **kwargs)
        if minimized is not None:
            branches[pos] = minimized
    unique: list[tuple[sp.Expr, ...]] = []
    for branch in branches:
        if branch not in unique:
            unique.append(branch)
    branches = unique
    keep = [True] * len(branches)
    for i, left in enumerate(branches):
        if not keep[i]:
            continue
        left_expr = _branch_expr(left)
        for j, right in enumerate(branches):
            if i == j or not keep[j]:
                continue
            right_expr = _branch_expr(right)
            _STATS.disjunction_checks += 1
            result = implies_by_cad(left_expr, right_expr, **kwargs)
            if result is True:
                keep[i] = False
                _STATS.removed_branches += 1
                break
    kept = [_branch_expr(branch) for branch, flag in zip(branches, keep, strict=True) if flag]
    _STATS.output_branches += len(kept)
    if not kept:
        return sp.false
    return simplify_boolean(sp.Or(*kept))


def minimize_disj_by_impl(
    expr: sp.Expr,
    *,
    max_dnf_branches: int = 4096,
    max_atoms: int = 8,
    max_vars: int = 3,
    max_degree: int = 4,
) -> sp.Expr:
    """Remove redundant DNF branches/atoms using one shared sign-invariant CAD.

    ``max_dnf_branches`` bounds Boolean expansion. ``max_atoms``, ``max_vars``,
    and ``max_degree`` bound when semantic CAD implication checks are attempted.
    Exceeding these guards leaves the formula conservatively less simplified;
    it never changes its truth set.
    """

    if max_dnf_branches < 1:
        raise ValueError("max_dnf_branches must be positive")
    if min(max_atoms, max_vars, max_degree) < 0:
        raise ValueError("implication CAD guards must be nonnegative")
    kwargs = {"max_atoms": max_atoms, "max_vars": max_vars, "max_degree": max_degree}
    _STATS.minimize_calls += 1
    branches = [tuple(branch) for branch in _dnf_branches(expr, max_branches=max_dnf_branches)]
    _STATS.input_branches += len(branches)
    if not branches:
        return sp.false
    try:
        context = _shared_truth_context(expr, **kwargs)
        if context is not None:
            return _minimize_disj_with_shared_cad(branches, context)
    except (
        sp.PolynomialError,
        TypeError,
        ValueError,
        ArithmeticError,
        NotImplementedError,
    ):
        pass
    return _minimize_disj_pairwise(branches, **kwargs)


__all__ = [
    "ImplicationMinimizationStats",
    "clear_implication_minimization_stats",
    "implication_minimization_stats",
    "implies_by_cad",
    "is_unsatisfiable_by_cad",
    "minimize_conj_by_impl",
    "minimize_disj_by_impl",
]
