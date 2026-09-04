from __future__ import annotations

from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DOCS = ROOT / "docs"


REQUIRED_REGION_DOCS = (
    "guides/region_representations.md",
    "guides/exact_vs_numerical_regions.md",
    "guides/cad_reuse.md",
    "tutorials/annulus.md",
    "tutorials/projection_3d.md",
    "tutorials/singular_curve.md",
    "tutorials/parameterized_region.md",
    "tutorials/optimization_over_region.md",
    "quality/region_capability_matrix.md",
    "quality/certificate_matrix.md",
    "quality/region_testing.md",
    "architecture/region_computation.md",
)


def test_region_documentation_pages_exist_and_are_navigated():
    nav = (ROOT / "mkdocs.yml").read_text(encoding="utf-8")
    for relative in REQUIRED_REGION_DOCS:
        assert (DOCS / relative).is_file(), relative
        assert relative in nav, f"{relative} is not in mkdocs navigation"


def test_exact_numerical_guide_states_the_important_nonclaims():
    text = (DOCS / "guides/exact_vs_numerical_regions.md").read_text(encoding="utf-8")
    assert "does not prove ambient isotopy" in text
    assert "algebraic-boundary singularity" in text
    assert "does not silently replace" in text


def test_representation_guide_explains_all_public_region_layers():
    text = (DOCS / "guides/region_representations.md").read_text(encoding="utf-8")
    for name in ("StandardRegion", "SemialgebraicRegion", "CADRegion", "StructuredCADCell"):
        assert name in text


def test_capability_and_certificate_matrices_record_planned_topology_gaps():
    capability = (DOCS / "quality/region_capability_matrix.md").read_text(encoding="utf-8")
    certificates = (DOCS / "quality/certificate_matrix.md").read_text(encoding="utf-8")
    assert "Betti numbers/homology" in capability
    assert "Certified isotopic triangulation" in capability
    assert "verified=None" in certificates
