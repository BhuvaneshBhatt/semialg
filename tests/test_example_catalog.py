"""The advertised catalog is complete and its application programs carry assertions."""

import ast
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def test_catalog_scripts():
    entries = json.loads((ROOT / "examples/catalog.json").read_text())
    names = [entry["file"] for entry in entries]
    actual = {
        p.relative_to(ROOT / "examples").as_posix()
        for p in (ROOT / "examples").rglob("*.py")
        if p.name != "run_all.py"
    }
    assert len(names) == len(set(names))
    assert set(names) == actual
    for row in entries:
        tree = ast.parse((ROOT / "examples" / row["file"]).read_text())
        assert any(isinstance(node, ast.Assert) for node in ast.walk(tree)), row["file"]


def test_ownership_manifest():
    import tomllib

    owners = tomllib.loads((ROOT / "docs/ownership.toml").read_text())
    assert owners["algroots"]["requires"] == []
    assert owners["semialg"]["requires"] == ["algroots"]
    assert owners["alggeom"]["requires"] == ["algroots"]
    assert owners["alggeom"]["optional"] == ["semialg"]
    assert "functional decomposition" in owners["alggeom"]["responsibilities"]


def test_backend_boundaries():
    import tomllib

    from packaging.requirements import Requirement

    package = tomllib.loads((ROOT / "pyproject.toml").read_text())["project"]
    name = package["name"]
    dependencies = {Requirement(value).name for value in package.get("dependencies", [])}
    forbidden = {"algroots": {"alggeom", "semialg"}, "semialg": {"alggeom"}, "alggeom": set()}[name]
    assert dependencies.isdisjoint(forbidden)
    for path in (ROOT / "src" / name).rglob("*.py"):
        imports = set()
        for node in ast.walk(ast.parse(path.read_text())):
            if isinstance(node, ast.Import):
                imports.update(alias.name.split(".")[0] for alias in node.names)
            elif isinstance(node, ast.ImportFrom) and node.level == 0 and node.module:
                imports.add(node.module.split(".")[0])
        assert imports.isdisjoint(forbidden), path


def test_tutorial_markers(tmp_path):
    import importlib.util

    spec = importlib.util.spec_from_file_location(
        "example_checker", ROOT / "tools/check_examples.py"
    )
    checker = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(checker)
    path = tmp_path / "tutorial.md"
    path.write_text("<!-- semialg-exec -->\n```python\nassert 2 + 2 == 4\n```\n")
    assert [code.strip() for code in checker.snippets(path)] == ["assert 2 + 2 == 4"]
