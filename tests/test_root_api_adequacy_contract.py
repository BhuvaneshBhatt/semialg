"""Semantic-evidence floor for every root-level public function.

This complements the file-level coverage manifest: a function must be exercised by
multiple independently named tests, and critical APIs require a third semantic
contract.  It intentionally does not treat line coverage as mathematical adequacy.
"""

from __future__ import annotations

import ast
import inspect
from collections import defaultdict
from pathlib import Path

import semialg

ROOT = Path(__file__).resolve().parents[1]
CRITICAL = {
    "apply_quantifiers",
    "cad",
    "equivalent",
    "find_negative_point",
    "find_negative_witness_fast",
    "implies",
    "is_satisfiable",
    "is_tautology",
    "polynomial_nonnegative",
    "prove_negative",
    "prove_nonnegative",
    "prove_nonpositive",
    "prove_nonzero",
    "prove_positive",
    "prove_zero",
    "real_algebraic_feasibility",
    "reduce_formula",
    "replay_certificate",
    "resolve_formula",
    "solvability_conditions",
    "solve_real_algebraic_set",
    "solve_semialgebraic",
}


def _named_test_evidence() -> dict[str, set[str]]:
    public = {name for name in semialg.__all__ if inspect.isfunction(getattr(semialg, name))}
    evidence: dict[str, set[str]] = defaultdict(set)
    for path in (ROOT / "tests").rglob("test_*.py"):
        tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
        stack: list[str] = []

        class Visitor(ast.NodeVisitor):
            def __init__(self, current_path: Path, current_stack: list[str]) -> None:
                self.path = current_path
                self.stack = current_stack

            def visit_FunctionDef(self, node: ast.FunctionDef) -> None:
                self.stack.append(node.name)
                self.generic_visit(node)
                self.stack.pop()

            visit_AsyncFunctionDef = visit_FunctionDef

            def visit_Call(self, node: ast.Call) -> None:
                name = None
                if isinstance(node.func, ast.Name):
                    name = node.func.id
                elif (
                    isinstance(node.func, ast.Attribute)
                    and isinstance(node.func.value, ast.Name)
                    and node.func.value.id in {"semialg", "sa"}
                ):
                    name = node.func.attr
                if name in public and self.stack and self.stack[-1].startswith("test_"):
                    relative = self.path.relative_to(ROOT).as_posix()
                    evidence[name].add(f"{relative}::{self.stack[-1]}")
                self.generic_visit(node)

        Visitor(path, stack).visit(tree)
    return evidence


def test_every_root_function_has_three_named_semantic_contracts() -> None:
    evidence = _named_test_evidence()
    public = {name for name in semialg.__all__ if inspect.isfunction(getattr(semialg, name))}
    thin = {name: sorted(evidence[name]) for name in public if len(evidence[name]) < 3}
    assert thin == {}


def test_critical_root_functions_have_a_third_independent_contract() -> None:
    evidence = _named_test_evidence()
    thin = {name: sorted(evidence[name]) for name in CRITICAL if len(evidence[name]) < 3}
    assert thin == {}
