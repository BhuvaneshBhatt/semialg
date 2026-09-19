"""Execute documentation blocks explicitly marked as semantic contracts."""

from __future__ import annotations

import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DOCS = ROOT / "docs"
EXECUTABLE_BLOCK = re.compile(
    r"<!--\s*semialg-exec\s*-->\s*```python\n(.*?)\n```",
    re.DOTALL,
)
PYTHON_BLOCK = re.compile(r"```python\n(.*?)\n```", re.DOTALL)


def test_marked_documentation_examples_execute_with_their_assertions() -> None:
    executed: list[str] = []
    for path in sorted(DOCS.rglob("*.md")):
        text = path.read_text(encoding="utf-8")
        for index, block in enumerate(EXECUTABLE_BLOCK.findall(text), start=1):
            label = f"{path.relative_to(ROOT)} block {index}"
            namespace = {"__name__": "__semialg_documentation_example__"}
            exec(compile(block, label, "exec"), namespace)
            executed.append(label)

    assert len(executed) >= 6


def test_readme_python_examples_execute_as_one_progressive_session() -> None:
    """The PyPI landing-page examples must remain executable in their written order."""

    readme = (ROOT / "README.md").read_text(encoding="utf-8")
    blocks = PYTHON_BLOCK.findall(readme)
    assert len(blocks) >= 8
    namespace = {"__name__": "__semialg_readme_example__"}
    for index, block in enumerate(blocks, start=1):
        label = f"README.md block {index}"
        exec(compile(block, label, "exec"), namespace)
