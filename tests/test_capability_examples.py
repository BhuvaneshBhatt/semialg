"""Every advertised capability links to an existing executable mathematical check."""

import ast
import json
from pathlib import Path

ROOT = Path(__file__).parents[1]


def test_capability_examples():
    matrix = (ROOT / "docs/feature_matrix.md").read_text(encoding="utf-8")
    titles = [
        line.split("|")[1].strip()
        for line in matrix.splitlines()
        if line.startswith("| ") and not line.startswith("| Capability")
    ]
    manifest = json.loads((ROOT / "docs/capability_examples.json").read_text(encoding="utf-8"))
    assert len(titles) == len(set(titles))
    assert set(titles) == manifest.keys()
    catalog = {
        entry["file"]
        for entry in json.loads((ROOT / "examples/catalog.json").read_text(encoding="utf-8"))
    }
    for title, name in manifest.items():
        path = ROOT / name
        tree = ast.parse(path.read_text(encoding="utf-8"))
        assert any(isinstance(node, ast.Assert) for node in ast.walk(tree)), (title, name)
        if name.startswith("examples/"):
            assert name.removeprefix("examples/") in catalog
        row = next(line for line in matrix.splitlines() if line.startswith("| " + title + " |"))
        assert "https://github.com/BhuvaneshBhatt/semialg/blob/main/" + name in row
