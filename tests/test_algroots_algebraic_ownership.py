from __future__ import annotations

import tomllib
from pathlib import Path

import algroots.border_basis
import algroots.rational_univariate
from packaging.requirements import Requirement

import semialg.algebraic.border_basis as semialg_border
from semialg.algebraic.rational_univariate import RationalUnivariateRepresentation, compute_rur


def test_algroots_owns_core_rur_representation_and_border_basis():
    assert (
        RationalUnivariateRepresentation
        is algroots.rational_univariate.RationalUnivariateRepresentation
    )
    assert semialg_border.BorderBasisResult is algroots.border_basis.BorderBasisResult
    assert semialg_border.compute_border_basis is algroots.border_basis.compute_border_basis


def test_quotient_parent():
    import sympy as sp
    from algroots.quotient import QuotientAlgebra

    from semialg.algebraic.rational_univariate.quotient import QuotientAlgebra as BackendQuotient

    assert BackendQuotient is QuotientAlgebra
    x, y = sp.symbols("x y")
    algebra = BackendQuotient.from_polynomials((x * x, y * y), (x, y))
    assert algebra.dimension == 4 and algebra.geometric_solution_count == 1


def test_algroots_dependency():
    metadata = tomllib.loads((Path(__file__).parents[1] / "pyproject.toml").read_text())
    requirements = [Requirement(value) for value in metadata["project"]["dependencies"]]
    dependency = next(requirement for requirement in requirements if requirement.name == "algroots")
    assert "0.3.0" in dependency.specifier
    assert "0.2.99" not in dependency.specifier


def test_compute_rur_still_returns_algroots_core_type():
    import sympy as sp

    x = sp.Symbol("x")
    result = compute_rur((x**2 - 2,), (x,))
    assert isinstance(result, algroots.rational_univariate.RationalUnivariateRepresentation)


def test_delegation_uses_only_namespaced_algroots_imports():
    import ast

    folder = Path(__file__).parents[1] / "src/semialg/algebraic"
    for path in folder.rglob("*.py"):
        for node in ast.walk(ast.parse(path.read_text())):
            if isinstance(node, ast.ImportFrom):
                assert node.module != "algroots", str(path)


def test_namespaced_backend_preserves_complex_and_multiple_roots():
    import sympy as sp

    from semialg.algebraic.rational_univariate import solve_rur_points, solve_rur_representation

    x = sp.Symbol("x")
    rep = compute_rur(((x**2 + 1) ** 2,), (x,))
    assert rep.dimension == 4 and rep.solution_count == 2
    assert solve_rur_points(rep) == ()
    coords = solve_rur_representation(rep, real=False)
    assert set(coords) == {(-sp.I,), (sp.I,)}
    assert {point.coordinates for point in solve_rur_points(rep, real=False)} == set(coords)
