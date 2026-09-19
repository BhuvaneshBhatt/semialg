from __future__ import annotations

from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parents[1]


def test_all_topic_and_reference_pages_are_in_mkdocs_navigation() -> None:
    config = yaml.safe_load((ROOT / "mkdocs.yml").read_text(encoding="utf-8"))

    def nav_paths(value: object) -> set[str]:
        if isinstance(value, str):
            return {value} if value.endswith(".md") else set()
        if isinstance(value, list):
            return set().union(*(nav_paths(item) for item in value))
        if isinstance(value, dict):
            return set().union(*(nav_paths(item) for item in value.values()))
        return set()

    nav_refs = nav_paths(config["nav"])
    pages = {
        str(path.relative_to(ROOT / "docs"))
        for path in (ROOT / "docs").rglob("*.md")
        if "examples" not in path.parts
    }
    assert pages <= nav_refs


def test_mkdocs_navigation_has_no_missing_or_duplicate_pages() -> None:
    config = yaml.safe_load((ROOT / "mkdocs.yml").read_text(encoding="utf-8"))
    paths: list[str] = []

    def collect(value: object) -> None:
        if isinstance(value, str) and value.endswith(".md"):
            paths.append(value)
        elif isinstance(value, list):
            for item in value:
                collect(item)
        elif isinstance(value, dict):
            for item in value.values():
                collect(item)

    collect(config["nav"])
    assert len(paths) == len(set(paths))
    assert all((ROOT / "docs" / path).is_file() for path in paths)
