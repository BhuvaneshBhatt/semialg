from __future__ import annotations

import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def test_all_topic_and_reference_pages_are_in_mkdocs_navigation() -> None:
    config = (ROOT / "mkdocs.yml").read_text(encoding="utf-8")
    nav_refs = set(re.findall(r":\s+([^\s]+\.md)\s*$", config, flags=re.MULTILINE))
    pages = {
        str(path.relative_to(ROOT / "docs"))
        for path in (ROOT / "docs").rglob("*.md")
        if "examples" not in path.parts
    }
    assert pages <= nav_refs
