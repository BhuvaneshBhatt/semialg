"""Documentation layers have distinct roles and one canonical solver-scope statement."""

from __future__ import annotations

from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DOCS = ROOT / "docs"


def test_concept_guide_reference_overviews_define_distinct_roles():
    concepts = (DOCS / "concepts" / "index.md").read_text(encoding="utf-8")
    guides = (DOCS / "guides" / "index.md").read_text(encoding="utf-8")
    reference = (DOCS / "reference" / "index.md").read_text(encoding="utf-8")

    assert "what a result means" in concepts
    assert "task-oriented" in guides
    assert "signatures" in reference and "return contracts" in reference


def test_root_and_transcendental_scope_guides_are_navigated():
    nav = (ROOT / "mkdocs.yml").read_text(encoding="utf-8")
    for page in (
        "guides/certified_algebraic_roots.md",
        "guides/transcendental_scope.md",
        "concepts/index.md",
        "guides/index.md",
        "reference/index.md",
    ):
        assert page in nav


def test_reference_delegates_root_workflow_explanation_to_the_guide():
    reference = (DOCS / "reference" / "algebraic.md").read_text(encoding="utf-8")
    assert "Certified algebraic roots" in reference
    assert "For algorithms, endpoint" in reference


def test_transcendental_scope_is_canonical_and_linked_from_solver_docs():
    solver = (DOCS / "transcendental_methods.md").read_text(encoding="utf-8")
    exactness = (DOCS / "concepts" / "exactness_and_certification.md").read_text(encoding="utf-8")
    assert "guides/transcendental_scope.md" in solver
    assert "../guides/transcendental_scope.md" in exactness
