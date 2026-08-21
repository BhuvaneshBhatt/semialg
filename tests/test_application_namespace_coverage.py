import ast
import inspect
from pathlib import Path

import semialg.applications as applications


def test_every_application_function_has_a_direct_behavioral_call():
    """Require each exported application function to be exercised directly."""

    public = {
        name for name in applications.__all__ if inspect.isfunction(getattr(applications, name))
    }
    called: set[str] = set()
    for path in Path(__file__).parent.glob("test_*.py"):
        tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
        for node in ast.walk(tree):
            if not isinstance(node, ast.Call):
                continue
            if isinstance(node.func, ast.Name) and node.func.id in public:
                called.add(node.func.id)
            elif isinstance(node.func, ast.Attribute) and node.func.attr in public:
                called.add(node.func.attr)
    assert not sorted(public - called), (
        f"application functions without direct tests: {sorted(public - called)}"
    )


def test_core_math_primitives_are_not_duplicated_in_applications_namespace():
    core_names = {
        "function_range",
        "semialgebraic_measure",
        "integrate_over_region",
        "semialgebraic_minimize",
        "semialgebraic_maximize",
        "cad",
        "qe_by_complete_cad",
    }
    assert not sorted(name for name in core_names if hasattr(applications, name))
