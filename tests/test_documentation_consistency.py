from __future__ import annotations

import re
from pathlib import Path

from semialg._api_policy import PRIMARY_EXPORTS

ROOT = Path(__file__).resolve().parents[1]
DOCS = ROOT / "docs"
MARKDOWN = [ROOT / "README.md", *DOCS.rglob("*.md")]

def test_relative_markdown_links_resolve():
    link_pattern = re.compile(r"\[[^\]]+\]\(([^)]+)\)")
    missing = []
    for path in MARKDOWN:
        text = path.read_text(encoding="utf-8")
        for raw in link_pattern.findall(text):
            target = raw.split("#", 1)[0]
            if not target or "://" in target or target.startswith("mailto:"):
                continue
            resolved = (path.parent / target).resolve()
            if not resolved.exists():
                missing.append((path.relative_to(ROOT), raw))
    assert not missing


def test_new_local_geometry_primary_apis_have_documentation():
    corpus = "\n".join(path.read_text(encoding="utf-8") for path in MARKDOWN)
    names = {
        "local_germ",
        "local_geometry",
        "local_sign_strata",
        "local_components",
        "curve_selection",
        "local_range",
        "local_bound",
        "parameter_strata",
        "vanishing_order",
        "contact_order",
        "blowup_charts",
        "local_image",
        "local_preimage",
        "path_independent",
    }
    assert names <= PRIMARY_EXPORTS.keys()
    for name in names:
        assert f"`{name}`" in corpus
