from __future__ import annotations

import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
README = ROOT / "README.md"
DOCS = ROOT / "docs"

REDUNDANT_ASSERT_OUTPUT = re.compile(
    r"^\s*assert\s+.+(?:\bis\s+True\b|\bis\s+False\b|==\s*(?:sp\.)?(?:true|false)\b).*\n\s*#\s*(?:True|False)\s*$",
    re.MULTILINE,
)


def test_reader_facing_examples_do_not_repeat_boolean_assertions_as_output() -> None:
    failures: list[str] = []
    for path in [README, *sorted(DOCS.rglob("*.md"))]:
        text = path.read_text(encoding="utf-8")
        if REDUNDANT_ASSERT_OUTPUT.search(text):
            failures.append(str(path.relative_to(ROOT)))
    assert not failures, "redundant assert + expected-output examples in: " + ", ".join(failures)


def test_computation_flow_diagrams_are_present_and_linked() -> None:
    page = (DOCS / "concepts" / "computation_flows.md").read_text(encoding="utf-8")
    for name in ("decision-flow.svg", "cad-flow.svg", "analysis-flow.svg"):
        asset = DOCS / "assets" / name
        assert asset.is_file() and asset.stat().st_size > 500
        assert f"../assets/{name}" in page
    assert (
        "raw.githubusercontent.com/BhuvaneshBhatt/semialg/main/docs/assets/decision-flow.svg"
        in README.read_text(encoding="utf-8")
    )
