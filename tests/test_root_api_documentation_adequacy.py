"""Architecture contracts for root-function documentation adequacy."""

from __future__ import annotations

import inspect
import subprocess
import sys
from pathlib import Path

try:
    import tomllib
except ModuleNotFoundError:
    import tomli as tomllib

import semialg

ROOT = Path(__file__).resolve().parents[1]
REGISTRY = ROOT / "docs" / "reference" / "root_function_documentation_adequacy.toml"
REPORT = ROOT / "docs" / "quality" / "root-api-documentation-adequacy.md"


def _registry():
    return tomllib.loads(REGISTRY.read_text(encoding="utf-8"))["functions"]


def test_documentation_adequacy_registry_covers_every_root_function():
    expected = {name for name in semialg.__all__ if inspect.isfunction(getattr(semialg, name))}
    assert set(_registry()) == expected


def test_documentation_adequacy_registry_has_explainable_statuses():
    allowed = {
        "purpose",
        "signature",
        "parameters",
        "return_semantics",
        "exactness",
        "algorithm",
        "limitations",
        "example",
    }
    for name, entry in _registry().items():
        required = set(entry["required"])
        satisfied = set(entry["satisfied"])
        missing = set(entry["missing"])
        assert required <= allowed, name
        assert missing == required - satisfied, name
        assert entry["status"] == ("adequate" if not missing else "needs-deepening"), name
        assert entry["reference"], name


def test_documentation_adequacy_artifacts_are_reproducible():
    before_registry = REGISTRY.read_bytes()
    before_report = REPORT.read_bytes()
    subprocess.run(
        [sys.executable, str(ROOT / "scripts" / "generate_root_api_documentation_adequacy.py")],
        cwd=ROOT,
        check=True,
    )
    assert REGISTRY.read_bytes() == before_registry
    assert REPORT.read_bytes() == before_report
