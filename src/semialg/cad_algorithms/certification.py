"""Independent replay of native CAD projection and geometric lifting data."""

from dataclasses import dataclass, replace

import sympy as sp

from ..algebraic.cache import sample_identity_key
from ..algebraic.comparison import compare_samples, sort_samples
from ..algebraic.roots import _ordinary_root_count, isolate_real_roots
from ..algebraic.samples import AlgebraicRoot, RationalSample
from ..algebraic.signs import sign_at_sample
from ..cache_utils import cache_disabled
from ..errors import ExactComputationUnavailable
from ..exact_arithmetic import exact_sign
from ..formula import And, Atom, BoolConst, Not, Or, parse_formula, to_sympy
from .decomposition import CompleteCAD, _build_stack, _stack_roots_over_point
from .lifting.stack import CADCell
from .polynomial_utils import polynomial_family_key
from .projection.collins import build_collins_proj_set
from .projection.ec_hypotheses import ec_parent_admissible, replay_ec_projection
from .projection.lazard_complete import build_lazard_complete_projection


@dataclass(frozen=True)
class CADVerification:
    """Replay outcome with the invariant proved and checked geometric data."""

    valid: bool
    invariant: str
    checked_cells: int
    checked_root_selectors: int
    checked_valuation_claims: int
    failures: tuple[str, ...] = ()


class _InvalidCAD(ValueError):
    pass


def _require(condition, message):
    if not condition:
        raise _InvalidCAD(message)


def _sign(expression, variables, samples):
    """Determine an exact rational-expression sign at a certified sample.

    Separating numerator and denominator detects poles instead of treating an
    undefined rational expression as a zero or an arbitrary sign.
    """
    if not variables:
        return exact_sign(expression)
    numerator, denominator = sp.fraction(sp.cancel(expression))
    denominator_sign = sign_at_sample(sp.Poly(denominator, *variables, domain="EX"), samples)
    _require(denominator_sign != 0, "a root representation has a pole on its parent")
    return sign_at_sample(sp.Poly(numerator, *variables, domain="EX"), samples) * denominator_sign


def _same_fiber_polynomial(left, right, variable, parents, samples):
    """Compare specialized polynomial coefficients up to a nonzero scalar."""
    first = list(sp.Poly(left, variable).all_coeffs())
    second = list(sp.Poly(right, variable).all_coeffs())
    while first and _sign(first[0], parents, samples) == 0:
        first.pop(0)
    while second and _sign(second[0], parents, samples) == 0:
        second.pop(0)
    return (
        bool(first)
        and len(first) == len(second)
        and all(
            _sign(a * second[0] - b * first[0], parents, samples) == 0
            for a, b in zip(first, second, strict=True)
        )
    )


def _verified_root(provided, reference, variables, prefix):
    """Bind a sample to a freshly isolated root, including its cached expression."""
    _require(isinstance(provided, AlgebraicRoot), "section sample is not an algebraic root")
    _require(provided.variable == variables[-1], "root polynomial uses the wrong fiber variable")
    _require(provided.root_index == reference.root_index, "root family rank was changed")
    _require(
        _same_fiber_polynomial(
            provided.polynomial.as_expr(),
            reference.polynomial.as_expr(),
            variables[-1],
            variables[:-1],
            prefix,
        ),
        "root polynomial does not represent its declared family",
    )
    context = reference.fiber_context
    if provided.fiber_context is not None:
        source = provided.fiber_context
        _require(context is not None, "unexpected parent context on a ground root")
        _require(
            source.variables == context.variables
            and source.fiber_variable == context.fiber_variable,
            "root context uses the wrong coordinate order",
        )
        _require(
            source.polynomial.as_expr() == context.polynomial.as_expr(),
            "root selector polynomial was changed",
        )
        _require(
            tuple(map(sample_identity_key, source.parent_samples))
            == tuple(map(sample_identity_key, prefix)),
            "root selector is attached to another parent",
        )
    elif prefix:
        raise _InvalidCAD("a native fiber section is missing its parent context")
    left, right = provided.interval.left, provided.interval.right
    _require(
        compare_samples(reference, RationalSample(left)) >= 0
        and compare_samples(reference, RationalSample(right)) <= 0,
        "isolating interval does not enclose the selected root",
    )
    if context is None:
        ground = sp.Poly(provided.polynomial.as_expr(), variables[-1], extension=True).sqf_part()
        count = (
            int(ground.eval(left) == 0)
            if left == right
            else _ordinary_root_count(ground, left, right)
        )
        _require(count == 1, "ground root interval is not isolating")
    elif left != right:
        _require(
            context.root_count(left, right) == 1
            and context.roots_left_of(left) == provided.root_index
            and context.roots_left_of(right) == provided.root_index + 1,
            "fiber interval has the wrong absolute root rank",
        )
    else:
        _require(
            context.sign_at_rational(left) == 0
            and context.roots_left_of(left) == provided.root_index + 1,
            "point interval is not the declared root",
        )
    expression = provided.root_expr
    if expression is not None:
        _require(
            not expression.free_symbols, "cached sample expression is not a ground algebraic value"
        )
        _require(
            exact_sign(expression - left) >= 0 and exact_sign(right - expression) >= 0,
            "cached root expression lies outside its interval",
        )
        _require(
            _sign(
                reference.polynomial.as_expr().subs(variables[-1], expression),
                variables[:-1],
                prefix,
            )
            == 0,
            "cached root expression does not satisfy its family",
        )
        if context is not None and expression == left:
            _require(
                context.roots_left_of(left) == provided.root_index + 1,
                "cached expression selects the excluded left endpoint",
            )
    return replace(provided, root_expr=None, fiber_context=context)


def _projection_replay(cad, formula, polynomials):
    """Recompute projection closure and provenance without completion flags."""
    tower = cad.tower
    variables = tower.variables
    if polynomials is not None:
        from .projection.collins import squarefree_basis

        expected_input = squarefree_basis(sp.Poly(poly, *variables) for poly in polynomials)
        _require(
            polynomial_family_key(expected_input)
            == polynomial_family_key(
                tuple(sp.Poly(poly, *variables) for poly in tower.original_polynomials)
            ),
            "input polynomial family does not match the proof statement",
        )
    if tower.metadata.get("ec_lifting"):
        _require(formula is not None, "EC replay requires the formula whose ECs are used")
        _require(
            replay_ec_projection(tower, formula),
            "EC implications or projection hypotheses do not replay",
        )
        uses_ec = any(
            plan.equational_constraint is not None for plan in tower.metadata["ec_levels"]
        )
        _require(
            not uses_ec or not tower.metadata.get("complete", False),
            "EC lifting cannot claim global sign invariance",
        )
        _require(
            tower.metadata.get("certified_invariant") in (None, "truth"),
            "EC proof claims an unsupported invariant",
        )
        scope = tower.metadata.get("truth_formulas", ())
        _require(
            all(item in (to_sympy(formula), sp.Not(to_sympy(formula))) for item in scope),
            "truth certificate advertises another formula",
        )
        return "truth" if uses_ec else "sign"
    if tower.metadata.get("projection") == "lazard-complete":
        expected = build_lazard_complete_projection(tower.original_polynomials, variables)
    else:
        expected = build_collins_proj_set(tower.original_polynomials, variables)
    _require(len(tower.levels) == len(expected.levels), "projection tower has missing levels")
    for actual, required in zip(tower.levels, expected.levels, strict=True):
        _require(
            actual.level == required.level and actual.variable == required.variable,
            "projection coordinate order was changed",
        )
        _require(
            polynomial_family_key(actual.polynomials)
            == polynomial_family_key(required.polynomials),
            "projection closure is missing or contains altered factors",
        )
        _require(
            actual.entries == required.entries, "projection operation provenance does not replay"
        )
    return "sign"


def verify_cad(cad, *, formula=None, polynomials=None, return_result=False):
    """Replay native ambient/recursive-EC CAD geometry and projection obligations.

    Memoization is bypassed. Every stack is reconstructed over verified parent
    samples; selector identities, rank intervals, neighboring bounds, coverage
    and recorded signs are checked. Lazard Taylor records and recursive EC
    implications are replayed as algebraic claims. Finite-variety and reduced-projection
    certificates lacking these obligations are outside this verifier.
    Completion flags and diagnostic counters are not used as proof.
    """
    if not isinstance(cad, CompleteCAD):
        raise TypeError("verify_cad requires a CompleteCAD")
    if formula is not None and not isinstance(formula, (And, Atom, BoolConst, Not, Or)):
        formula = parse_formula(formula)
    checked = selectors = valuation_count = 0
    invariant = "unknown"
    try:
        with cache_disabled():
            variables = cad.tower.variables
            if not variables:
                _require(
                    set(cad.cells_by_level) == {0} and len(cad.cells) == 1,
                    "invalid zero-coordinate CAD",
                )
                cell = cad.cells[0]
                _require(
                    not cad.tower.levels and cell == CADCell(level=0, index=(), sample=()),
                    "invalid zero-coordinate cell",
                )
                statements = tuple(cad.tower.original_polynomials)
                if polynomials is not None:
                    statements += tuple(polynomials)
                _require(
                    not any(
                        sp.sympify(
                            item.as_expr() if isinstance(item, sp.Poly) else item
                        ).free_symbols
                        for item in statements
                    ),
                    "zero-coordinate proof contains an unassigned polynomial variable",
                )
                _require(
                    formula is None or not to_sympy(formula).free_symbols,
                    "zero-coordinate formula contains an unassigned variable",
                )
                for item in statements:
                    expression = item.as_expr() if isinstance(item, sp.Poly) else sp.sympify(item)
                    _require(expression.is_real is not False, "nonreal zero-coordinate coefficient")
                    exact_sign(expression)
                result = CADVerification(True, "sign", 1, 0, 0)
                return result if return_result else True
            invariant = _projection_replay(cad, formula, polynomials)
            _require(
                set(cad.cells_by_level) == set(range(1, len(variables) + 1)),
                "CAD tree has missing or extra levels",
            )
            plans = (
                {p.level: p for p in cad.tower.metadata.get("ec_levels", ())}
                if invariant == "truth"
                else {}
            )
            lazard = cad.tower.metadata.get("projection") in {"lazard", "lazard-complete"}
            expected_records = []
            verified = {}
            raw_samples = {}
            parents = [None]
            for level, variable in enumerate(variables, 1):
                actual_cells = cad.cells_by_level[level]
                _require(
                    len({cell.index for cell in actual_cells}) == len(actual_cells),
                    "duplicate CAD cell index",
                )
                consumed = set()
                next_parents = []
                for parent in parents:
                    parent_index = None if parent is None else parent.index
                    prefix = () if parent is None else parent.sample
                    family = cad.tower.level(level).polynomials
                    admissible = ec_parent_admissible(cad.tower, parent) if plans else True
                    if plans and not admissible:
                        family = ()
                    elif level in plans and plans[level].equational_constraint is not None:
                        family = (sp.Poly(plans[level].equational_constraint, *variables[:level]),)
                    if plans and admissible and level in plans:
                        for coefficient in plans[level].degree_coefficients:
                            _require(
                                _sign(coefficient, variables[: level - 1], prefix) != 0,
                                "EC degree-invariance condition fails",
                            )
                    if level == 1:
                        roots = sort_samples(
                            tuple(root for poly in family for root in isolate_real_roots(poly))
                        )
                    else:
                        before_records = len(expected_records)
                        roots = _stack_roots_over_point(
                            family,
                            variables,
                            prefix,
                            variable,
                            use_lazard=lazard,
                            valuation_records=expected_records,
                            delineate_nullification=level < len(variables) and not lazard,
                            base_is_point=all(position % 2 == 1 for position in parent.index),
                        )
                        expected_records[before_records:] = [
                            replace(record, parent_samples=raw_samples[parent.index])
                            for record in expected_records[before_records:]
                        ]
                    expected_stack = _build_stack(parent, roots, level, cad.tower)
                    stack = sorted(
                        (cell for cell in actual_cells if cell.parent_index == parent_index),
                        key=lambda cell: cell.stack_position,
                    )
                    _require(
                        len(stack) == len(expected_stack),
                        "stack does not cover its full real fiber",
                    )
                    for actual, expected in zip(stack, expected_stack, strict=True):
                        _require(
                            actual.index == expected.index
                            and actual.level == level
                            and actual.kind == expected.kind
                            and actual.stack_position == expected.stack_position,
                            "cell index, kind or stack order was changed",
                        )
                        _require(
                            len(actual.sample) == level,
                            "sample coordinate count does not match its level",
                        )
                        _require(
                            actual.sample[:-1]
                            == (() if parent is None else raw_samples[parent.index]),
                            "cell sample has a different parent prefix",
                        )
                        _require(
                            actual.root_index == expected.root_index
                            and actual.defining_polynomial_key == expected.defining_polynomial_key
                            and actual.section_polynomial == expected.section_polynomial,
                            "section selector identity was changed",
                        )
                        if actual.kind == "section":
                            value = _verified_root(
                                actual.sample[-1], expected.sample[-1], variables[:level], prefix
                            )
                            selectors += 1
                        else:
                            value = actual.sample[-1]
                            _require(
                                isinstance(value, RationalSample), "sector sample is not rational"
                            )
                            left, right = expected.interval
                            _require(
                                (left is None or compare_samples(left, value) < 0)
                                and (right is None or compare_samples(value, right) < 0),
                                "sector sample is outside its open interval",
                            )
                        sample = (*prefix, value)
                        _require(
                            set(actual.signs) == set(expected.signs), "sign family was changed"
                        )
                        for poly in cad.tower.level(level).polynomials:
                            from .polynomial_utils import polynomial_key

                            _require(
                                actual.signs[polynomial_key(poly)] == sign_at_sample(poly, sample),
                                "recorded sign is incorrect",
                            )
                        verified[actual.index] = replace(actual, sample=sample)
                        raw_samples[actual.index] = actual.sample
                        next_parents.append(verified[actual.index])
                        consumed.add(actual.index)
                        checked += 1
                    for position, cell in enumerate(stack):
                        if cell.kind == "section":
                            bounds = (cell.sample[-1], cell.sample[-1])
                        else:
                            bounds = (
                                None if position == 0 else stack[position - 1].sample[-1],
                                None
                                if position == len(stack) - 1
                                else stack[position + 1].sample[-1],
                            )
                        _require(
                            cell.interval == bounds, "cell bounds are not its neighboring sections"
                        )
                _require(
                    len(consumed) == len(actual_cells), "CAD contains cells with an unknown parent"
                )
                parents = next_parents
            actual_records = cad.tower.metadata.get("lazard_lifts", ())
            _require(
                tuple(actual_records) == tuple(expected_records),
                "Lazard valuation/residue claims do not replay",
            )
            valuation_count = len(expected_records)
        result = CADVerification(True, invariant, checked, selectors, valuation_count)
    except (
        _InvalidCAD,
        ExactComputationUnavailable,
        ArithmeticError,
        ValueError,
        TypeError,
        KeyError,
        sp.polys.polyerrors.BasePolynomialError,
    ) as exc:
        result = CADVerification(False, invariant, checked, selectors, valuation_count, (str(exc),))
    return result if return_result else result.valid
