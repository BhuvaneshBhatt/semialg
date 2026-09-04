"""Contract: every shipped Python example must execute successfully."""

from __future__ import annotations

import runpy
from pathlib import Path

import matplotlib.pyplot as plt
import pytest

ROOT = Path(__file__).resolve().parents[1]
EXAMPLES = tuple(sorted((ROOT / "examples").rglob("*.py")))


def test_examples_directory_is_nonempty() -> None:
    assert EXAMPLES


@pytest.mark.parametrize(
    "script",
    EXAMPLES,
    ids=lambda path: str(path.relative_to(ROOT / "examples")).replace("/", "::"),
)
def test_example_executes(script: Path) -> None:
    """Execute each example as a top-level Python program."""

    try:
        runpy.run_path(str(script), run_name="__main__")
    finally:
        plt.close("all")
