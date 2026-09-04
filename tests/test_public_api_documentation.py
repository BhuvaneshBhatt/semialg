"""Documentation coverage and contract regressions for the public API."""

from __future__ import annotations

import inspect
import re
from pathlib import Path

import tomllib

import semialg

ROOT = Path(__file__).resolve().parents[1]
DOCS = ROOT / "docs"
PUBLIC_INDEX = DOCS / "reference" / "public_api.md"

_PRIMARY_REFERENCE_FILES = (
    DOCS / "reference" / "decision_and_qe.md",
    DOCS / "reference" / "cad.md",
    DOCS / "reference" / "solving_and_sampling.md",
    DOCS / "reference" / "optimization_and_range.md",
    DOCS / "reference" / "regions.md",
    DOCS / "reference" / "integration_and_moments.md",
    DOCS / "reference" / "algebraic.md",
    DOCS / "reference" / "parameters.md",
    DOCS / "reference" / "convexity_backend.md",
)


def test_public_api_index_matches_all_exports_exactly():
    text = PUBLIC_INDEX.read_text(encoding="utf-8")
    indexed = re.findall(r"^\| `([^`]+)` \|", text, flags=re.MULTILINE)

    assert len(indexed) == len(set(indexed))
    assert set(indexed) == set(semialg.__all__)


def test_public_functions_and_classes_have_docstrings():
    missing = []
    for name in semialg.__all__:
        obj = getattr(semialg, name)
        if inspect.getdoc(obj):
            continue
        missing.append(name)

    assert missing == []


def test_primary_reference_families_state_the_contract():
    required = (
        "## Family contract",
        "**Mathematical return.**",
        "**Exactness and certification.**",
        "**Algorithm",
        "**Complexity and limitations.**",
    )
    for path in _PRIMARY_REFERENCE_FILES:
        text = path.read_text(encoding="utf-8")
        for marker in required:
            assert marker in text, f"{marker!r} missing from {path.relative_to(ROOT)}"


def test_cad_reference_emphasizes_variable_ordering():
    text = (DOCS / "reference" / "cad.md").read_text(encoding="utf-8")
    assert (
        "## Variable ordering is a first-class performance decision" in text
        or "### Variable ordering is a first-class performance decision" in text
    )
    assert "orders-of-magnitude differences in runtime and memory" in text
    assert "Quantifiers restrict what may be reordered" in text


def test_documentation_relative_links_resolve():
    broken = []
    pattern = re.compile(r"\[[^\]]+\]\(([^)]+)\)")

    for path in DOCS.rglob("*.md"):
        text = path.read_text(encoding="utf-8")
        for target in pattern.findall(text):
            if "://" in target or target.startswith("#") or target.startswith("mailto:"):
                continue

            relative = target.split("#", 1)[0]
            if not relative:
                continue

            resolved = (path.parent / relative).resolve()
            if not resolved.exists():
                broken.append((str(path.relative_to(ROOT)), target))

    assert broken == []


def _section_for_anchor(text: str, anchor: str) -> str | None:
    """Return the H2 section whose GitHub/MkDocs-style slug matches ``anchor``."""

    for match in re.finditer(r"^## (.+)$", text, flags=re.MULTILINE):
        title = match.group(1).strip().lower()
        slug = re.sub(r"[^a-z0-9 -]", "", title).replace(" ", "-")
        if slug != anchor:
            continue
        next_heading = re.search(r"^## ", text[match.end() :], flags=re.MULTILINE)
        end = match.end() + next_heading.start() if next_heading else len(text)
        return text[match.start() : end]
    return None


def test_every_primary_api_maps_to_a_substantive_documentation_section():
    """Every primary root API must map to a maintained topical section."""

    manifest_path = DOCS / "reference" / "primary_api_manifest.toml"
    manifest = tomllib.loads(manifest_path.read_text(encoding="utf-8"))["apis"]
    expected = {name for name in semialg.__all__}
    assert set(manifest) == expected

    failures = []
    for name, target in manifest.items():
        relative, anchor = target.split("#", 1)
        path = DOCS / relative
        if not path.exists():
            failures.append((name, "missing page", target))
            continue
        section = _section_for_anchor(path.read_text(encoding="utf-8"), anchor)
        if section is None:
            failures.append((name, "missing section", target))
            continue
        if len(section.strip()) < 200:
            failures.append((name, "section too short", target))
            continue
        pattern = re.compile(rf"(?<![A-Za-z0-9_]){re.escape(name)}(?![A-Za-z0-9_])")
        if not pattern.search(section):
            failures.append((name, "name absent from mapped section", target))

    assert failures == []
