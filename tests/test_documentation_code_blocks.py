from __future__ import annotations

import ast
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DOCS = ROOT / "docs"
PYTHON_BLOCK = re.compile(r"```python\n(.*?)\n```", re.DOTALL)


def test_every_python_documentation_block_is_syntactically_valid() -> None:
    failures: list[str] = []
    count = 0
    for path in sorted(DOCS.rglob("*.md")):
        text = path.read_text(encoding="utf-8")
        for index, block in enumerate(PYTHON_BLOCK.findall(text), start=1):
            count += 1
            try:
                ast.parse(block, filename=f"{path}#python-block-{index}")
            except SyntaxError as exc:
                failures.append(f"{path.relative_to(ROOT)} block {index}: {exc}")
    assert count > 0
    assert not failures, "\n".join(failures)
