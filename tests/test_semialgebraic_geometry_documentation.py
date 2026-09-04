"""Contracts for the introductory semialgebraic-geometry documentation section."""

from __future__ import annotations

from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DOCS = ROOT / "docs" / "semialgebraic_geometry"
MKDOCS = ROOT / "mkdocs.yml"

PAGES = {
    "introduction.md": (
        "# Introduction to semialgebraic geometry",
        "Tarski-Seidenberg",
        "SemialgebraicRegion",
    ),
    "theory.md": (
        "# Concepts and theory",
        "quantifier elimination",
        "Euler characteristic",
        "Singular and regular points",
    ),
    "algorithms.md": (
        "# Algorithms and techniques",
        "Cylindrical algebraic decomposition",
        "Virtual substitution",
        "Exact optimization",
    ),
    "applications.md": (
        "# Applications of semialgebraic geometry",
        "Robotics and motion planning",
        "Formal verification",
        "Global optimization",
    ),
}


def test_semialgebraic_geometry_learning_pages_exist_and_cover_core_topics():
    for filename, markers in PAGES.items():
        path = DOCS / filename
        assert path.exists()
        text = path.read_text(encoding="utf-8")
        for marker in markers:
            assert marker in text, f"{marker!r} missing from {filename}"


def test_semialgebraic_geometry_section_is_in_site_navigation():
    nav = MKDOCS.read_text(encoding="utf-8")
    assert "- Semialgebraic geometry:" in nav
    for filename in PAGES:
        assert f"semialgebraic_geometry/{filename}" in nav
