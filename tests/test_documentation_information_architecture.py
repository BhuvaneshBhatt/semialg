"""Navigation and graph behavior contracts without fixed editorial phrasing."""

import re
from pathlib import Path

import sympy as sp

from semialg.decision import equivalent
from semialg.function_verification import verify_real_function

ROOT = Path(__file__).resolve().parents[1]


def test_navigation_targets():
    for page in (
        "concepts/index.md",
        "guides/index.md",
        "reference/index.md",
        "package_ownership.md",
        "worked_applications.md",
        "function_verification.md",
    ):
        assert (ROOT / "docs" / page).is_file()
        assert page in (ROOT / "mkdocs.yml").read_text(encoding="utf-8")


def test_primary_workflow_links():
    for page in ("index.md", "api_overview.md", "guides/choosing_an_api.md"):
        path = ROOT / "docs" / page
        links = re.findall(r"\[[^\]]+\]\(([^)]+)\)", path.read_text(encoding="utf-8"))
        targets = {
            (path.parent / target.split("#", 1)[0]).resolve()
            for target in links
            if "://" not in target
        }
        assert (ROOT / "docs" / "reference/index.md").resolve() in targets or (
            ROOT / "docs" / "guides/choosing_an_api.md"
        ).resolve() in targets
        assert all(target.is_file() for target in targets)


def test_graph_existence_and_uniqueness():
    x, y = sp.symbols("x y", real=True)
    graph = sp.And(sp.Eq(y * y, x), y >= 0)
    result = verify_real_function(graph, input=x, output=y, domain=x >= 0)
    assert result.total and result.single_valued
    assert equivalent(result.totality_condition, x >= 0, (x,))
    assert result.nonunique_condition is sp.false
