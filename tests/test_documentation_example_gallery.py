"""Regression tests for the executable documentation example gallery."""

from __future__ import annotations

import ast
import re
from pathlib import Path

import sympy as sp

from semialg import (
    integrate_over_region,
    is_equal,
    semialgebraic_measure,
    semialgebraic_projection,
    singular_locus,
    tangent_cone,
    tangent_space,
)

ROOT = Path(__file__).resolve().parents[1]
DOC_EXAMPLES = ROOT / "docs" / "examples"
GALLERY = ROOT / "examples" / "gallery"


def _gallery_stems() -> list[str]:
    return sorted(path.stem for path in GALLERY.glob("*.py"))


def test_gallery_has_documented_executable_examples():
    stems = _gallery_stems()

    assert len(stems) == 18
    assert all((DOC_EXAMPLES / f"{stem}.md").exists() for stem in stems)

    index = (DOC_EXAMPLES / "index.md").read_text(encoding="utf-8")
    for stem in stems:
        assert f"({stem}.md)" in index


def test_all_gallery_scripts_parse_as_python():
    for path in sorted(GALLERY.glob("*.py")):
        ast.parse(path.read_text(encoding="utf-8"), filename=str(path))


def test_documented_code_matches_executable_companion():
    pattern = re.compile(r"```python\n(.*?)\n```", flags=re.DOTALL)

    for script in sorted(GALLERY.glob("*.py")):
        markdown = (DOC_EXAMPLES / f"{script.stem}.md").read_text(encoding="utf-8")
        match = pattern.search(markdown)
        assert match is not None, script.stem

        executable = script.read_text(encoding="utf-8")
        executable = executable.split("\n\n", 1)[1].strip()
        assert match.group(1).strip() == executable


def test_gallery_core_code_paths_smoke():
    x, y, a = sp.symbols("x y a", real=True)

    projection = semialgebraic_projection(
        (y >= x**2) & (y <= 1),
        eliminate=[y],
        variables=[x, y],
    )
    assert is_equal(projection, (x >= -1) & (x <= 1), [x])

    assert semialgebraic_measure(x**2 + y**2 <= 1, [x, y]) == sp.pi

    cusp = y**2 - x**3
    singular = singular_locus([cusp], [x, y])
    assert singular.subs({x: 0, y: 0}) is sp.true

    space = tangent_space([cusp], {x: 0, y: 0}, [x, y])
    cone = tangent_cone([cusp], {x: 0, y: 0}, [x, y])
    assert space.dimension == 2
    assert cone.certified

    result = integrate_over_region(
        1,
        x**2 <= a,
        [x],
        parameters=[a],
        return_stratified=True,
    )
    assert result.select({a: 4}) == 4
