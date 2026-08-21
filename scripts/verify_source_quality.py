from __future__ import annotations

import ast
import re
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SOURCE_DIRS = (ROOT / "src", ROOT / "tests")
TEXT_PATHS = (
    *SOURCE_DIRS,
    ROOT / "docs",
    ROOT / "README.md",
    ROOT / "pyproject.toml",
)
MAX_VARIABLE_LENGTH = 24
PHASE_NAME_RE = re.compile(
    r"(?:^|[_-])(phase|milestone|round\d*|pass\d*|hardening)(?:[_-]|$)", re.I
)


def _python_files() -> list[Path]:
    return sorted(path for base in SOURCE_DIRS for path in base.rglob("*.py"))


def _text_files() -> list[Path]:
    files: list[Path] = []
    for path in TEXT_PATHS:
        if path.is_file():
            files.append(path)
        elif path.is_dir():
            files.extend(
                item
                for item in path.rglob("*")
                if item.is_file() and item.suffix in {".py", ".md", ".toml", ".yml", ".yaml"}
            )
    return sorted(set(files))


def _duplicate_bindings(
    body: list[ast.stmt], scope: str = "module"
) -> list[tuple[str, str, list[int]]]:
    """Return repeated definition or assignment names within lexical scopes."""

    bindings: dict[str, list[int]] = defaultdict(list)
    found: list[tuple[str, str, list[int]]] = []
    for node in body:
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)):
            bindings[node.name].append(node.lineno)
        elif isinstance(node, (ast.Assign, ast.AnnAssign)):
            targets = node.targets if isinstance(node, ast.Assign) else [node.target]
            for target in targets:
                if isinstance(target, ast.Name):
                    bindings[target.id].append(node.lineno)
    found.extend((scope, name, lines) for name, lines in bindings.items() if len(lines) > 1)
    for node in body:
        if isinstance(node, ast.ClassDef):
            found.extend(_duplicate_bindings(node.body, f"{scope}.{node.name}"))
    return found


def _long_variables(tree: ast.Module) -> list[tuple[str, int]]:
    """Return bound variable and parameter names longer than the style limit."""

    found: set[tuple[str, int]] = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.arg) and len(node.arg) > MAX_VARIABLE_LENGTH:
            found.add((node.arg, node.lineno))
        elif (
            isinstance(node, ast.Name)
            and isinstance(node.ctx, ast.Store)
            and len(node.id) > MAX_VARIABLE_LENGTH
        ):
            found.add((node.id, node.lineno))
    return sorted(found, key=lambda item: (item[1], item[0]))


def main() -> int:
    """Validate repository-level source hygiene rules not covered by Ruff."""

    failures: list[str] = []

    for path in _text_files():
        relative = path.relative_to(ROOT)
        for line_number, line in enumerate(path.read_text().splitlines(), start=1):
            if line.rstrip() != line:
                failures.append(f"{relative}:{line_number}: trailing whitespace")

    for path in _python_files():
        relative = path.relative_to(ROOT)
        tree = ast.parse(path.read_text(), filename=str(relative))
        for scope, name, lines in _duplicate_bindings(tree.body):
            failures.append(f"{relative}: overridden binding {scope}.{name} at {lines}")
        for name, line in _long_variables(tree):
            failures.append(
                f"{relative}:{line}: variable {name!r} exceeds {MAX_VARIABLE_LENGTH} characters"
            )

    for base in (ROOT / "src", ROOT / "tests", ROOT / "docs"):
        for path in base.rglob("*"):
            if path.is_file() and PHASE_NAME_RE.search(path.stem):
                failures.append(f"{path.relative_to(ROOT)}: phase-shaped filename")

    if failures:
        print("Source quality checks failed:")
        for failure in failures:
            print(f"- {failure}")
        return 1

    print("Source quality checks passed.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
