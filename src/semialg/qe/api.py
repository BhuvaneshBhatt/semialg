"""Exact real QE dispatch using specialist solvers & quantifier-guided CAD."""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass, field

import sympy as sp

from ..algebraic.groebner_utils import compute_groebner_basis
from ..formula import ParsedPrenexFormula, parse_formula
from ..normalization import (
    conjuncts,
    normalize_formula,
    normalize_problem_variables,
    normalize_variables,
)
from ..presolve import fourier_motzkin_eliminate, presolve_semialgebraic
from ..quantifiers import split_quantifiers
from ..simplify.boolean import simplify_boolean
from ..structural_keys import symbol_identity_key


@dataclass(frozen=True)
class QuantifierEliminationResult:
    """Certified quantifier-free result together with dispatcher provenance."""

    formula: sp.Expr
    free_variables: tuple[sp.Symbol, ...]
    quantified_variables: tuple[sp.Symbol, ...]
    quantifiers: tuple[tuple[str, sp.Symbol], ...]
    method: str
    certified: bool = True
    presolve_substitutions: tuple[tuple[sp.Symbol, sp.Expr], ...] = ()
    variable_blocks: tuple[tuple[sp.Symbol, ...], ...] = ()
    backend_result: object | None = None
    notes: tuple[str, ...] = ()
    metadata: dict[str, object] = field(default_factory=dict)


def _normalize_prefix(
    formula: object,
    quantifiers: Sequence[tuple[str, sp.Symbol | str]] | None,
) -> tuple[tuple[tuple[str, sp.Symbol], ...], sp.Expr]:
    expr = normalize_formula(formula)
    embedded, matrix = split_quantifiers(expr)
    if quantifiers is not None and embedded:
        raise ValueError(
            "quantifiers must be supplied either in the formula or explicitly, not both"
        )
    raw = tuple(quantifiers) if quantifiers is not None else embedded
    context = (matrix, *(v for _, v in embedded))
    prefix: list[tuple[str, sp.Symbol]] = []
    seen: set[sp.Symbol] = set()
    for raw_name, raw_variable in raw:
        name = str(raw_name).lower()
        if name not in {"exists", "forall"}:
            raise ValueError(f"unsupported quantifier: {raw_name!r}")
        variable = normalize_variables((raw_variable,), *context, append_context_symbols=False)[0]
        if variable in seen:
            raise ValueError(f"variable {variable} is quantified more than once")
        seen.add(variable)
        prefix.append((name, variable))
    return tuple(prefix), normalize_formula(matrix)


def _make_parsed(
    matrix: sp.Expr,
    prefix: tuple[tuple[str, sp.Symbol], ...],
    variables: Sequence[sp.Symbol | str] | None,
) -> ParsedPrenexFormula:
    quantified = tuple(v for _, v in prefix)
    all_variables = normalize_problem_variables(variables, matrix, quantified)
    # Keep free coordinates first and preserve quantifier-prefix order after them.
    qset = set(quantified)
    free = tuple(v for v in all_variables if v not in qset)
    missing_q = tuple(v for v in quantified if v not in all_variables)
    ordered = free + tuple(v for v in quantified if v not in free) + missing_q
    ordered = tuple(dict.fromkeys(ordered))
    return ParsedPrenexFormula(ordered, prefix, parse_formula(matrix), matrix)


def _all_existential(prefix: Sequence[tuple[str, sp.Symbol]]) -> bool:
    return bool(prefix) and all(name == "exists" for name, _ in prefix)


def _conjunctive_block_factorization(
    matrix: sp.Expr,
    blocks: Sequence[tuple[sp.Symbol, ...]],
    prefix: Sequence[tuple[str, sp.Symbol]],
) -> tuple[tuple[sp.Expr, tuple[tuple[str, sp.Symbol], ...]], ...] | None:
    """Factor a conjunction into independent quantified variable blocks.

    This is an exact logical factorization: each conjunct must involve variables
    from at most one incidence block.  Quantifiers whose variables are absent
    from a factor are dropped.  The transformation is valid for both existential
    and universal prefixes because the matrix is a conjunction and the blocks
    are variable-disjoint.
    """
    if len(blocks) <= 1 or not isinstance(matrix, sp.And):
        return None
    factors: list[list[sp.Expr]] = [[] for _ in blocks]
    constants: list[sp.Expr] = []
    for atom in conjuncts(matrix):
        hits = [i for i, block in enumerate(blocks) if set(block) & atom.free_symbols]
        if not hits:
            constants.append(atom)
        elif len(hits) == 1:
            factors[hits[0]].append(atom)
        else:
            return None
    out: list[tuple[sp.Expr, tuple[tuple[str, sp.Symbol], ...]]] = []
    if constants:
        out.append((sp.And(*constants), ()))
    for block, atoms in zip(blocks, factors, strict=True):
        if not atoms:
            continue
        block_set = set(block)
        local_prefix = tuple(item for item in prefix if item[1] in block_set)
        out.append((sp.And(*atoms), local_prefix))
    return tuple(out) if len(out) > 1 else None


def _augment_with_groebner_equalities(
    matrix: sp.Expr, variables: Sequence[sp.Symbol]
) -> tuple[sp.Expr, tuple[sp.Expr, ...]]:
    """Add exact equality-ideal consequences without performing real projection.

    The added Groebner-basis equations generate the same equality ideal, so
    conjoining them is equivalence-preserving over the reals. They are used
    only to expose safe affine substitutions and zero-dimensional structure to
    existing certified preprocessors. No elimination-ideal condition is ever
    treated as an existential projection.
    """
    if not isinstance(matrix, sp.And):
        return matrix, ()
    atoms = conjuncts(matrix)
    equalities = tuple(
        sp.expand(atom.lhs - atom.rhs) for atom in atoms if isinstance(atom, sp.Equality)
    )
    if len(equalities) < 2 or not variables:
        return matrix, ()
    try:
        basis = compute_groebner_basis(equalities, tuple(variables), order="lex")
    except (sp.PolynomialError, ValueError, TypeError):
        return matrix, ()
    consequences = tuple(sp.expand(poly.as_expr()) for poly in basis.polys)
    existing = set(equalities) | {-eq for eq in equalities}
    added = tuple(eq for eq in consequences if eq not in existing and -eq not in existing)
    if not added:
        return matrix, ()
    return sp.And(matrix, *(sp.Eq(eq, 0) for eq in added), evaluate=False), added


def quantifier_eliminate(
    formula,
    quantifiers=None,
    *,
    variables=None,
    strategy="auto",
    return_result=False,
    max_cells=None,
    timeout=None,
):
    """Eliminate real quantifiers using exact specialists and partial CAD.

    ``formula`` is a Boolean matrix with explicit ``quantifiers``, or contains
    ``Exists``/``ForAll`` nodes. ``variables`` supplies the coordinate order;
    free coordinates precede the quantified prefix internally. ``strategy``
    selects automatic specialist dispatch, explicit partial lifting, or a
    specific backend; ``"collins"`` requests full CAD.

    Return the quantifier-free SymPy formula. ``return_result=True`` retains a
    ``QuantifierEliminationResult`` with backend metadata. Free-parameter CAD
    retains truth on every free cell while lifting quantified stacks as needed.

    ``max_cells`` counts cells across attempts and cache returns. ``timeout``
    sets a cooperative deadline shared through preprocessing and CAD work.
    Exhaustion raises ``ResourceLimitError`` without an incomplete formula;
    one symbolic operation can overrun the deadline.
    """
    from ..cad_algorithms.resource_limits import CADBudget, budget_scope, current_budget

    budget = (
        CADBudget(max_cells, timeout)
        if max_cells is not None or timeout is not None
        else current_budget()
    )
    with budget_scope(budget):
        if budget is not None:
            budget.check()
        result = _quantifier_eliminate(
            formula,
            quantifiers,
            variables=variables,
            strategy=strategy,
            return_result=return_result,
        )
        if budget is not None:
            budget.check()
        return result


def _quantifier_eliminate(
    formula: object,
    quantifiers: Sequence[tuple[str, sp.Symbol | str]] | None = None,
    *,
    variables: Sequence[sp.Symbol | str] | None = None,
    strategy: str = "auto",
    return_result: bool = False,
) -> sp.Expr | QuantifierEliminationResult:
    """Dispatch a normalized prefix within the caller's resource scope."""
    prefix, matrix = _normalize_prefix(formula, quantifiers)
    active_symbols = matrix.free_symbols
    prefix = tuple((kind, variable) for kind, variable in prefix if variable in active_symbols)
    parsed = _make_parsed(matrix, prefix, variables)
    quantified = tuple(v for _, v in prefix)
    qset = set(quantified)
    free = tuple(v for v in parsed.vars if v not in qset)
    if not prefix:
        result = QuantifierEliminationResult(
            formula=matrix,
            free_variables=free,
            quantified_variables=(),
            quantifiers=(),
            method="quantifier-free",
            notes=("input contains no quantifiers",),
        )
        return result if return_result else result.formula

    notes: list[str] = []
    presolved = presolve_semialgebraic(
        matrix,
        parsed.vars,
        eliminate=quantified if _all_existential(prefix) else (),
        detect_zero_dimensional=True,
    )
    working_matrix = presolved.formula
    working_variables = presolved.variables
    substitutions = presolved.substitutions
    working_prefix = prefix

    # Exploit the incidence decomposition discovered by presolve instead of
    # sending independent blocks through one larger  CAD.
    factored = _conjunctive_block_factorization(
        working_matrix, presolved.variable_blocks, working_prefix
    )
    if factored is not None:
        pieces: list[sp.Expr] = []
        block_methods: list[str] = []
        for factor, local_prefix in factored:
            if not local_prefix:
                pieces.append(factor)
                block_methods.append("quantifier-free")
                continue
            local = quantifier_eliminate(
                factor, local_prefix, strategy=strategy, return_result=True
            )
            pieces.append(local.formula)
            block_methods.append(local.method)
        out = simplify_boolean(sp.And(*pieces))
        result = QuantifierEliminationResult(
            formula=out,
            free_variables=tuple(sorted(out.free_symbols, key=symbol_identity_key)),
            quantified_variables=quantified,
            quantifiers=prefix,
            method="independent-blocks",
            presolve_substitutions=substitutions,
            variable_blocks=presolved.variable_blocks,
            notes=("eliminated independent incidence blocks separately",),
            metadata={"block_methods": tuple(block_methods)},
        )
        return result if return_result else out
    if presolved.substitutions:
        removed = {v for v, _ in presolved.substitutions}
        working_prefix = tuple(item for item in prefix if item[1] not in removed)
        notes.append("eliminated globally safe affine equalities before QE")
    if not working_prefix:
        out = simplify_boolean(working_matrix)
        result = QuantifierEliminationResult(
            formula=out,
            free_variables=tuple(sorted(out.free_symbols, key=symbol_identity_key)),
            quantified_variables=quantified,
            quantifiers=prefix,
            method="affine-presolve",
            presolve_substitutions=substitutions,
            variable_blocks=presolved.variable_blocks,
            notes=tuple(notes),
        )
        return result if return_result else out

    if _all_existential(working_prefix):
        remaining_q = tuple(v for _, v in working_prefix)
        fm = fourier_motzkin_eliminate(working_matrix, remaining_q)
        if fm is not None and set(fm[1]) == set(remaining_q):
            out = simplify_boolean(fm[0])
            result = QuantifierEliminationResult(
                formula=out,
                free_variables=tuple(sorted(out.free_symbols, key=symbol_identity_key)),
                quantified_variables=quantified,
                quantifiers=prefix,
                method="fourier-motzkin",
                presolve_substitutions=substitutions,
                variable_blocks=presolved.variable_blocks,
                notes=(*notes, "used exact Fourier–Motzkin elimination"),
            )
            return result if return_result else out

    groebner_added: tuple[sp.Expr, ...] = ()
    if _all_existential(working_prefix):
        augmented, groebner_added = _augment_with_groebner_equalities(
            working_matrix, working_variables
        )
        if groebner_added:
            second = presolve_semialgebraic(
                augmented,
                working_variables,
                eliminate=tuple(v for _, v in working_prefix),
                detect_zero_dimensional=True,
            )
            working_matrix = second.formula
            working_variables = second.variables
            substitutions += second.substitutions
            if second.substitutions:
                removed = {v for v, _ in second.substitutions}
                working_prefix = tuple(item for item in working_prefix if item[1] not in removed)
                notes.append("used exact Groebner equality consequences for affine presolve")
            if not working_prefix:
                out = simplify_boolean(working_matrix)
                result = QuantifierEliminationResult(
                    formula=out,
                    free_variables=tuple(sorted(out.free_symbols, key=symbol_identity_key)),
                    quantified_variables=quantified,
                    quantifiers=prefix,
                    method="groebner-affine-presolve",
                    presolve_substitutions=substitutions,
                    variable_blocks=presolved.variable_blocks,
                    notes=tuple(notes),
                    metadata={"groebner_equalities_added": groebner_added},
                )
                return result if return_result else out

    working = _make_parsed(working_matrix, working_prefix, working_variables)
    # Imported lazily so low-level CAD/decomposition modules can import
    # ``semialg.qe.complete`` without creating a package-initialization cycle.
    from ..solve.reduce import reduce_formula

    solved = reduce_formula(working, domain="reals", strategy=strategy, return_result=True)
    out = normalize_formula(solved.result)
    backend = solved.metadata.get("qe_result")
    result = QuantifierEliminationResult(
        formula=out,
        free_variables=tuple(v for v in working.vars if v not in {qv for _, qv in working_prefix}),
        quantified_variables=quantified,
        quantifiers=prefix,
        method=solved.method,
        presolve_substitutions=substitutions,
        variable_blocks=presolved.variable_blocks,
        backend_result=backend,
        notes=tuple(notes),
        metadata={**dict(solved.metadata), "groebner_equalities_added": groebner_added},
    )
    return result if return_result else out


def project_region(
    region: object,
    eliminate: Sequence[sp.Symbol | str],
    *,
    variables: Sequence[sp.Symbol | str] | None = None,
    strategy: str = "auto",
    return_result: bool = False,
) -> sp.Expr | QuantifierEliminationResult:
    """Return the exact real existential projection of ``region``.

    This is the geometric facade over :func:`quantifier_eliminate`. projection
    is existential quantification, not a separate approximate operation.
    """
    matrix = normalize_formula(region)
    eliminated = normalize_variables(eliminate, matrix, append_context_symbols=False)
    return quantifier_eliminate(
        matrix,
        tuple(("exists", variable) for variable in eliminated),
        variables=variables,
        strategy=strategy,
        return_result=return_result,
    )


__all__ = ["QuantifierEliminationResult", "project_region", "quantifier_eliminate"]
