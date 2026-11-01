"""Exact existence, uniqueness and composition checks for real semialgebraic graphs."""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass

import sympy as sp

from .decision import equivalent, implies, is_satisfiable
from .decision.solution import EquivalenceResult, ImplicationResult, SatisfiabilityResult
from .function_graph import SemialgebraicFunctionGraph
from .qe import QuantifierEliminationResult, quantifier_eliminate
from .simplify.atoms import normalize_atoms


def _symbols(symbols, name):
    symbols = tuple(symbols)
    if any(not isinstance(s, sp.Symbol) or s.is_real is False for s in symbols) or len(
        set(symbols)
    ) != len(symbols):
        raise ValueError(f"{name} must contain distinct real Symbols")
    return symbols


@dataclass(frozen=True)
class RealCompositionVerification:
    """Exact comparison of real relations on a declared input domain.

    With ``check_functions=True``, totality and single-valuedness are checked
    separately from graph equality. ``eliminations`` records exact QE provenance, not a
    claim that every backend exposes a separately replayable proof object.
    """

    target_formula: sp.Expr
    composition_formula: sp.Expr
    domain: sp.Expr
    eliminations: tuple[QuantifierEliminationResult, ...]
    comparison: EquivalenceResult
    function_checks: tuple[RealFunctionVerification, ...] = ()

    @property
    def equivalent(self) -> bool:
        """Whether the two real relations agree on the declared domain."""
        return self.comparison.equivalent

    @property
    def certified(self) -> bool:
        """Whether the comparison used certified exact elimination paths."""
        return all(result.certified for result in self.eliminations) and all(
            result.certified for result in self.function_checks
        )

    @property
    def function_verified(self) -> bool | None:
        """Whether both graphs are total, single-valued functions; None if unchecked."""
        return (
            all(result.total and result.single_valued for result in self.function_checks)
            if self.function_checks
            else None
        )


def _fresh_graph(graph, expected_target, allowed):
    if isinstance(graph, SemialgebraicFunctionGraph):
        if graph.target != expected_target or not graph.supported:
            raise ValueError("graph must be supported and use the declared target symbol")
        formula = graph.formula
        auxiliary = _symbols(graph.auxiliary_variables, "auxiliary variables")
    else:
        formula = sp.sympify(graph)
        auxiliary = ()
    if set(auxiliary) & set(allowed) or formula.free_symbols - {*allowed, *auxiliary}:
        raise ValueError("graph contains undeclared symbols or captures a declared variable")
    renaming = {v: sp.Dummy("composition_aux", real=True) for v in auxiliary}
    return formula.xreplace(renaming), tuple(renaming.values())


def _quantified_graph_comparison(left, right, variables):
    """Decide both graph differences after existential output presolve."""
    results = tuple(
        quantifier_eliminate(
            difference,
            [("exists", v) for v in variables],
            variables=variables,
            return_result=True,
        )
        for difference in (sp.And(left, sp.Not(right)), sp.And(right, sp.Not(left)))
    )
    decisions = tuple(normalize_atoms(result.formula) for result in results)
    if any(decision not in (sp.true, sp.false) for decision in decisions):
        raise NotImplementedError("closed graph comparison did not yield an exact Boolean decision")
    lnr, rnl = (decision is sp.true for decision in decisions)
    direction = (
        "both" if lnr and rnl else "lhs_implies_rhs" if lnr else "rhs_implies_lhs" if rnl else None
    )
    comparison = EquivalenceResult(
        not lnr and not rnl,
        left,
        right,
        variables,
        failed_direction=direction,
        method="quantified_symmetric_difference",
        diagnostics={"lhs_not_rhs": results[0].method, "rhs_not_lhs": results[1].method},
    )
    return comparison, results


def verify_real_composition(
    target_graph,
    outer_graph,
    inner_graph,
    *,
    input: sp.Symbol,
    output: sp.Symbol,
    intermediate: sp.Symbol,
    domain=sp.true,
    parameters: Sequence[sp.Symbol] = (),
    check_functions: bool = False,
) -> RealCompositionVerification:
    """Verify ``target(output,input) = exists y: outer(output,y)&inner(y,input)``.

    Graphs may be semialgebraic formulas or SemialgebraicFunctionGraph objects.
    Graph-object auxiliaries and the intermediate variable are freshly bound to
    prevent capture. Encode rational poles and branch selectors in the graphs;
    graph builders already preserve these for supported expressions. ``domain``
    may use only the input and declared parameters. Unsupported operations raise
    their ordinary exact-solver errors; grid sampling is never a proof fallback.
    """
    parameters = _symbols(parameters, "parameters")
    _symbols((input, output, intermediate, *parameters), "composition variables and parameters")
    domain = sp.sympify(domain)
    if domain.free_symbols - {input, *parameters}:
        raise ValueError("domain may contain only the input and parameters")
    target, target_aux = _fresh_graph(target_graph, output, (input, output, *parameters))
    outer, outer_aux = _fresh_graph(outer_graph, output, (intermediate, output, *parameters))
    inner, inner_aux = _fresh_graph(inner_graph, intermediate, (input, intermediate, *parameters))
    fresh_middle = sp.Dummy("composition_middle", real=True)
    composition = sp.And(outer, inner).xreplace({intermediate: fresh_middle})
    bound = (fresh_middle, *outer_aux, *inner_aux)
    variables = (*parameters, input, output)
    target_qe = quantifier_eliminate(
        target,
        [("exists", v) for v in target_aux],
        variables=(*variables, *target_aux),
        return_result=True,
    )
    composition_qe = quantifier_eliminate(
        composition,
        [("exists", v) for v in bound],
        variables=(*variables, *bound),
        return_result=True,
    )
    return compare_real_graphs(
        target_qe,
        composition_qe,
        input=input,
        output=output,
        parameters=parameters,
        domain=domain,
        check_functions=check_functions,
        find_counterexample=not parameters,
    )


def compare_real_graphs(
    target_graph,
    composition_graph,
    *,
    input: sp.Symbol,
    output: sp.Symbol,
    parameters: Sequence[sp.Symbol] = (),
    domain=sp.true,
    check_functions: bool = False,
    find_counterexample: bool = False,
) -> RealCompositionVerification:
    """Compare two real graphs on an input domain, optionally checking function obligations.

    Each graph may be a formula, a SemialgebraicFunctionGraph, or an exact
    QuantifierEliminationResult. Projection results retain their elimination
    provenance; other graphs have their auxiliary variables existentially
    eliminated. Parameters are free real inputs quantified by the comparison.
    Equality of relations does not imply totality or single-valuedness.
    Set ``find_counterexample=True`` to request an exact witness when the
    graphs differ; otherwise closed existential decisions establish equality.
    """
    parameters = _symbols(parameters, "parameters")
    _symbols((input, output, *parameters), "graph variables and parameters")
    domain = sp.sympify(domain)
    if domain.free_symbols - {input, *parameters}:
        raise ValueError("domain may contain only the input and parameters")
    variables = (*parameters, input, output)

    def project(graph):
        if isinstance(graph, QuantifierEliminationResult):
            formula = sp.sympify(graph.formula)
            if formula.free_symbols - set(variables):
                raise ValueError("projected graph contains undeclared symbols")
            return graph
        formula, auxiliary = _fresh_graph(graph, output, variables)
        return quantifier_eliminate(
            formula,
            [("exists", v) for v in auxiliary],
            variables=(*variables, *auxiliary),
            return_result=True,
        )

    target_qe, composition_qe = (project(graph) for graph in (target_graph, composition_graph))
    target_formula, composition_formula = (
        normalize_atoms(result.formula) for result in (target_qe, composition_qe)
    )
    left, right = (
        normalize_atoms(sp.And(domain, formula))
        for formula in (target_formula, composition_formula)
    )
    comparison_eliminations = ()
    if left != right and not find_counterexample:
        # Existential output presolve avoids constructing large sample points
        # when graph equalities already determine the output coordinates.
        comparison, comparison_eliminations = _quantified_graph_comparison(left, right, variables)
    else:
        comparison = equivalent(left, right, variables, return_result=True)
    checks = (
        tuple(
            verify_real_function(
                formula, input=input, output=output, parameters=parameters, domain=domain
            )
            for formula in (target_formula, composition_formula)
        )
        if check_functions
        else ()
    )
    return RealCompositionVerification(
        target_formula,
        composition_formula,
        domain,
        (target_qe, composition_qe, *comparison_eliminations),
        comparison,
        checks,
    )


@dataclass(frozen=True)
class RealFunctionVerification:
    """Exact totality and uniqueness checks on a declared input/parameter domain."""

    formula: sp.Expr
    domain: sp.Expr
    totality_condition: sp.Expr
    nonunique_condition: sp.Expr
    totality: ImplicationResult
    nonuniqueness: SatisfiabilityResult
    eliminations: tuple[QuantifierEliminationResult, ...]

    @property
    def total(self):
        """Whether every input in the domain has an output."""
        return bool(self.totality)

    @property
    def single_valued(self):
        """Whether each input has at most one output."""
        return not bool(self.nonuniqueness)

    @property
    def certified(self):
        """Whether the function obligations used certified elimination paths."""
        return all(result.certified for result in self.eliminations)


def verify_real_function(
    graph,
    *,
    input: sp.Symbol,
    output: sp.Symbol,
    parameters: Sequence[sp.Symbol] = (),
    domain=sp.true,
):
    """Prove existence and uniqueness of graph outputs on the requested domain.

    Quantifier elimination computes the input loci with at least one output and
    with two distinct outputs. These conditions and exact decision provenance
    are returned even when totality or single-valuedness fails. A relation with
    no output can be single-valued but is not total. Empty domains are vacuous.
    """
    parameters = _symbols(parameters, "parameters")
    _symbols((input, output, *parameters), "function variables and parameters")
    domain = sp.sympify(domain)
    if domain.free_symbols - {input, *parameters}:
        raise ValueError("domain may contain only the input and parameters")
    formula, aux = _fresh_graph(graph, output, (input, output, *parameters))
    inputs = (*parameters, input)
    projected = quantifier_eliminate(
        formula, [("exists", v) for v in aux], variables=(*inputs, output, *aux), return_result=True
    )
    formula = normalize_atoms(projected.formula)
    coverage = quantifier_eliminate(
        formula, [("exists", output)], variables=(*inputs, output), return_result=True
    )
    other = sp.Dummy("other_function_output", real=True)
    two_outputs = sp.And(formula, formula.xreplace({output: other}), sp.Ne(output, other))
    nonunique = quantifier_eliminate(
        two_outputs,
        [("exists", output), ("exists", other)],
        variables=(*inputs, output, other),
        return_result=True,
    )
    totality = implies(
        normalize_atoms(domain), normalize_atoms(coverage.formula), inputs, return_result=True
    )
    multiple = is_satisfiable(
        normalize_atoms(sp.And(domain, nonunique.formula)), inputs, return_result=True
    )
    return RealFunctionVerification(
        formula,
        domain,
        normalize_atoms(coverage.formula),
        normalize_atoms(nonunique.formula),
        totality,
        multiple,
        (projected, coverage, nonunique),
    )


__all__ = [
    "compare_real_graphs",
    "RealCompositionVerification",
    "RealFunctionVerification",
    "verify_real_composition",
    "verify_real_function",
]
