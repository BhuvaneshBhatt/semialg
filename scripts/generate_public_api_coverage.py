"""Generate the public-API test-coverage inventory.

The generated TOML is committed so reviewers can see coverage ownership and evidence
without running Python.  Tests render it again and require a byte-for-byte match.
"""

from __future__ import annotations

import ast
import inspect
import sys
from collections import defaultdict
from pathlib import Path

try:
    import tomllib
except ModuleNotFoundError:
    import tomli as tomllib

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
MANIFEST = ROOT / "tests" / "public_api_coverage.toml"
DOCUMENTATION_MANIFEST = ROOT / "docs" / "reference" / "primary_api_manifest.toml"

if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

import semialg  # noqa: E402
from semialg._public_api import PUBLIC_EXPORTS  # noqa: E402

ALLOWED_CONTRACTS = {
    "behavioral-call",
    "construction",
    "distribution-metadata",
    "factory-return",
    "inheritance",
    "public-surface",
    "production-path",
    "raising-site",
}
ALLOWED_GAPS = {"raising-site"}

CRITICAL_APIS = {
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
    "root_count_conditions",
    "solvability_conditions",
    "solve_real_algebraic_set",
    "solve_semialgebraic",
}

INDIRECT_TYPES = {
    "AffineBoxClip": ("factory-return", "clip_affine_subspace_to_box"),
    "CADRegion": ("factory-return", "as_cad_region"),
    "Geometry": ("inheritance", "Interval"),
    "StandardRegion": ("inheritance", "Interval"),
    "CriticalValueImage": ("factory-return", "critical_value_image"),
    "ActiveConstraintResult": ("factory-return", "active_constraints"),
    "ComponentConstraintDescription": ("factory-return", "component_constraint_descriptions"),
    "ImpliedPolynomialInequality": ("factory-return", "implied_polynomial_inequality"),
    "IrreducibleAlgebraicComponent": ("factory-return", "irreducible_components"),
    "LocalBranchGeometry": ("factory-return", "local_branch_geometry"),
    "NonnegativeCombinationCertificate": ("factory-return", "nonnegative_combination_certificate"),
    "ParameterizationGeometry": ("factory-return", "parameterization_geometry"),
    "PolynomialConstraintSystem": ("factory-return", "polynomial_constraints"),
    "PolynomialMapImplicitizationResult": ("factory-return", "implicitize_polynomial_map"),
    "ReducedAlgebraicVariety": ("factory-return", "certified_radicalization"),
    "ReducedComponentSingularLocus": ("factory-return", "reduced_component_singular_loci"),
    "RedundantPolynomialInequality": ("factory-return", "redundant_polynomial_inequalities"),
    "StratifiedSingularGeometry": ("factory-return", "stratified_singular_geometry"),
    "ZariskiClosureResult": ("factory-return", "zariski_closure"),
    "RegionConversion": ("factory-return", "convert_region"),
    "BranchTangentGeometry": ("factory-return", "local_branch_geometry"),
    "ComponentIntersectionGeometry": ("factory-return", "local_branch_geometry"),
    "LocalDimensionStratum": ("factory-return", "local_dimension_strata"),
    "MinimalPrimeIntersection": ("factory-return", "minimal_prime_intersections"),
    "SingularGeometryStratum": ("factory-return", "stratified_singular_geometry"),
    "PolynomialConstraint": ("factory-return", "polynomial_constraints"),
    "PolynomialConstraintClause": ("factory-return", "polynomial_constraints"),
    "MixedCellTetrahedralization": ("factory-return", "tetrahedralize_cell"),
    "MixedMeshTetrahedralization": ("factory-return", "tetrahedralize_cells"),
    "PolytopeDecomposition": ("factory-return", "triangulate_polytope"),
    "QuantifierEliminationResult": ("factory-return", "quantifier_eliminate"),
}

EXCEPTION_PROFILES = {
    "ResourceLimitError": {
        "contracts": ["public-surface", "inheritance", "raising-site", "production-path"],
        "tests": ["tests/test_public_type_contracts.py", "tests/test_tunable_limits.py"],
        "production": "cad",
    },
    "SemialgError": {
        "contracts": ["public-surface", "inheritance", "raising-site", "production-path"],
        "tests": ["tests/test_public_type_contracts.py", "tests/test_tunable_limits.py"],
        "production": "ResourceLimitError",
    },
    "SemialgStrategyFailure": {
        "contracts": ["public-surface", "inheritance", "raising-site", "production-path"],
        "tests": ["tests/test_public_type_contracts.py", "tests/test_tunable_limits.py"],
        "production": "ResourceLimitError",
    },
    "UnsupportedFragmentError": {
        "contracts": ["public-surface", "inheritance", "construction"],
        "tests": ["tests/test_public_result_contracts.py"],
        "gaps": ["raising-site"],
    },
}


def _toml_string(value: str) -> str:
    return '"' + value.replace("\\", "\\\\").replace('"', '\\"') + '"'


def _direct_call_evidence(names: set[str]) -> dict[str, list[str]]:
    evidence: dict[str, set[str]] = defaultdict(set)
    for path in sorted((ROOT / "tests").rglob("test_*.py")):
        relative = path.relative_to(ROOT).as_posix()
        tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
        for node in ast.walk(tree):
            if not isinstance(node, ast.Call):
                continue
            name: str | None = None
            if isinstance(node.func, ast.Name):
                name = node.func.id
            elif (
                isinstance(node.func, ast.Attribute)
                and isinstance(node.func.value, ast.Name)
                and node.func.value.id == "semialg"
            ):
                name = node.func.attr
            if name in names:
                evidence[name].add(relative)
    return {name: sorted(paths) for name, paths in evidence.items()}


def _owner(name: str, documentation_target: str) -> str:
    page = documentation_target.split("#", 1)[0]
    if page.endswith("decision_and_qe.md"):
        return "decision"
    if page.endswith("solving_and_sampling.md"):
        return "solving"
    if page.endswith("optimization_and_range.md"):
        return "optimization"
    if page.endswith("integration_and_moments.md"):
        return "integration"
    if page.endswith("algebraic.md"):
        return "algebraic"
    if page.endswith("convexity_backend.md"):
        return "function-analysis"
    if page.endswith("errors_and_failure_modes.md"):
        return "errors"
    if page.endswith("package_metadata.md"):
        return "package"
    if page.endswith("regions.md") or page.endswith("cad.md"):
        return "geometry"
    if page.endswith("root_api_usage.md"):
        module = PUBLIC_EXPORTS[name]
        if any(token in module for token in ("algebraic", "incidence")):
            return "algebraic"
        if "optimization" in module:
            return "optimization"
        if "integrat" in module or "moment" in module or "measure" in module:
            return "integration"
        if any(token in module for token in ("decision", "reasoning", "parameters")):
            return "decision"
        if any(token in module for token in ("solv", "sampling")):
            return "solving"
        if any(token in module for token in ("function", "convex")):
            return "function-analysis"
        return "geometry"
    raise ValueError(f"no coverage owner for documentation target {documentation_target!r}")


def _kind(name: str) -> str:
    if name == "__version__":
        return "metadata"
    obj = getattr(semialg, name)
    if inspect.isfunction(obj):
        return "function"
    if inspect.isclass(obj):
        return "exception" if issubclass(obj, BaseException) else "type"
    raise TypeError(f"unsupported public export {name}: {obj!r}")


def _risk(name: str, kind: str, owner: str) -> str:
    if name in CRITICAL_APIS:
        return "critical"
    if kind == "exception":
        return "high"
    if kind in {"metadata", "type"}:
        return "low" if owner == "package" else "medium"
    if owner in {"algebraic", "decision", "geometry", "integration", "optimization"}:
        return "high"
    return "medium"


def render_manifest() -> str:
    documentation = tomllib.loads(DOCUMENTATION_MANIFEST.read_text(encoding="utf-8"))["apis"]
    exported = set(semialg.__all__)
    if set(documentation) != exported:
        missing = sorted(exported - set(documentation))
        stale = sorted(set(documentation) - exported)
        raise ValueError(f"documentation manifest mismatch: missing={missing}, stale={stale}")

    direct = _direct_call_evidence(exported)
    lines = [
        "# Generated by scripts/generate_public_api_coverage.py; do not edit manually.",
        "schema_version = 1",
        "",
    ]
    for name in sorted(exported, key=lambda item: (item.lower(), item)):
        kind = _kind(name)
        owner = _owner(name, documentation[name])
        risk = _risk(name, kind, owner)
        production: str | None = None
        gaps: list[str] = []
        if name == "__version__":
            contracts = ["public-surface", "distribution-metadata"]
            tests = ["tests/test_package_metadata.py"]
        elif kind == "exception":
            profile = EXCEPTION_PROFILES[name]
            contracts = profile["contracts"]
            tests = profile["tests"]
            production = profile.get("production")
            gaps = profile.get("gaps", [])
        elif name in INDIRECT_TYPES:
            relationship, production = INDIRECT_TYPES[name]
            contracts = ["public-surface", relationship, "production-path"]
            tests = [
                "tests/test_public_api_direct_contracts.py",
                "tests/test_public_type_contracts.py",
            ]
        else:
            contracts = [
                "public-surface",
                "behavioral-call" if kind == "function" else "construction",
            ]
            tests = direct.get(name, [])[:2]
        if not tests:
            raise ValueError(f"public API {name!r} has no test evidence")

        lines.extend(
            [
                f"[apis.{_toml_string(name)}]",
                f"kind = {_toml_string(kind)}",
                f"owner = {_toml_string(owner)}",
                f"risk = {_toml_string(risk)}",
                "contracts = [" + ", ".join(map(_toml_string, contracts)) + "]",
                "tests = [" + ", ".join(map(_toml_string, tests)) + "]",
            ]
        )
        if production is not None:
            lines.append(f"production = {_toml_string(production)}")
        if gaps:
            lines.append("gaps = [" + ", ".join(map(_toml_string, gaps)) + "]")
        if kind == "function":
            lines.append(f"evidence_file_count = {len(direct[name])}")
            lines.append('coverage_status = "expanded"')
        lines.append("")
    return "\n".join(lines)


def main() -> None:
    MANIFEST.write_text(render_manifest(), encoding="utf-8")
    print(f"wrote {MANIFEST.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
