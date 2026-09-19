"""Contracts for the committed public-API test-coverage inventory."""

from __future__ import annotations

import importlib.util
import inspect
from pathlib import Path

import pytest

try:
    import tomllib
except ModuleNotFoundError:
    import tomli as tomllib

import semialg

ROOT = Path(__file__).resolve().parents[1]
MANIFEST_PATH = ROOT / "tests" / "public_api_coverage.toml"
GENERATOR_PATH = ROOT / "scripts" / "generate_public_api_coverage.py"
ALLOWED_KINDS = {"exception", "function", "metadata", "type"}
ALLOWED_RISKS = {"critical", "high", "medium", "low"}
COVERAGE_STATUSES = {"expanded"}
pytestmark = pytest.mark.contract


def _generator_module():
    spec = importlib.util.spec_from_file_location("public_api_coverage_generator", GENERATOR_PATH)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _manifest() -> dict[str, dict[str, object]]:
    return tomllib.loads(MANIFEST_PATH.read_text(encoding="utf-8"))["apis"]


def test_public_api_coverage_manifest_is_current() -> None:
    """A public export addition/removal must update the committed inventory."""

    generator = _generator_module()
    parsed = tomllib.loads(MANIFEST_PATH.read_text(encoding="utf-8"))
    assert parsed["schema_version"] == 1
    assert MANIFEST_PATH.read_text(encoding="utf-8") == generator.render_manifest()
    assert set(_manifest()) == set(semialg.__all__)


def test_every_manifest_entry_has_valid_classification_and_test_evidence() -> None:
    """Every API has a valid owner, risk, contract set, and existing test file."""

    generator = _generator_module()
    for name, entry in _manifest().items():
        kind = entry["kind"]
        assert kind in ALLOWED_KINDS, name
        assert entry["risk"] in ALLOWED_RISKS, name
        assert isinstance(entry["owner"], str) and entry["owner"], name
        contracts = entry["contracts"]
        assert contracts and set(contracts) <= generator.ALLOWED_CONTRACTS, name
        assert set(entry.get("gaps", ())) <= generator.ALLOWED_GAPS, name
        tests = entry["tests"]
        assert tests, name
        for relative in tests:
            path = ROOT / relative
            assert path.is_file(), f"{name}: missing evidence file {relative}"
            assert path.is_relative_to(ROOT / "tests"), f"{name}: evidence is outside tests"

        if name == "__version__":
            assert kind == "metadata"
        elif inspect.isfunction(getattr(semialg, name)):
            assert kind == "function"
            assert "behavioral-call" in contracts
            assert entry["coverage_status"] in COVERAGE_STATUSES
            assert entry["evidence_file_count"] >= 2
            assert len(tests) >= 2
        else:
            assert inspect.isclass(getattr(semialg, name)), name


def test_indirect_types_and_exceptions_declare_production_contracts() -> None:
    """Factory/base/exception types must identify how users encounter them."""

    entries = _manifest()
    assert entries["CADRegion"]["production"] == "as_cad_region"
    assert entries["AffineBoxClip"]["production"] == "clip_affine_subspace_to_box"
    assert entries["Geometry"]["production"] == "IntervalRegion"
    assert entries["StandardRegion"]["production"] == "IntervalRegion"
    for name in ("CADRegion", "AffineBoxClip", "Geometry", "StandardRegion"):
        assert "production-path" in entries[name]["contracts"]
        assert "tests/test_public_type_contracts.py" in entries[name]["tests"]
    for name in ("ResourceLimitError", "SemialgError", "SemialgStrategyFailure"):
        assert entries[name]["kind"] == "exception"
        assert {"inheritance", "raising-site", "production-path"} <= set(entries[name]["contracts"])
        assert "tests/test_public_type_contracts.py" in entries[name]["tests"]
        assert entries[name]["production"]
    unsupported = entries["UnsupportedFragmentError"]
    assert {"inheritance", "construction"} <= set(unsupported["contracts"])
    assert unsupported["gaps"] == ["raising-site"]


def test_manifest_function_evidence_contains_a_direct_call() -> None:
    """Declared function evidence must contain a direct public-function call."""

    generator = _generator_module()
    function_names = {
        name for name in semialg.__all__ if inspect.isfunction(getattr(semialg, name))
    }
    evidence = generator._direct_call_evidence(function_names)
    for name in function_names:
        declared = set(_manifest()[name]["tests"])
        assert declared <= set(evidence[name]), f"{name}: declared evidence has no direct call"
        assert _manifest()[name]["evidence_file_count"] == len(evidence[name])


def test_no_public_function_relies_on_one_test_file() -> None:
    """Require independent file-level evidence for every public function."""

    function_entries = {
        name: entry for name, entry in _manifest().items() if entry["kind"] == "function"
    }
    thin = sorted(
        name for name, entry in function_entries.items() if entry["evidence_file_count"] < 2
    )
    assert thin == []
