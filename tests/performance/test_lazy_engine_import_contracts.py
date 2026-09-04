import json
import os
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]


def _modules_after(code: str) -> set[str]:
    env = dict(os.environ)
    env["PYTHONPATH"] = str(ROOT / "src")
    completed = subprocess.run(
        [sys.executable, "-c", code],
        cwd=ROOT,
        env=env,
        check=True,
        capture_output=True,
        text=True,
    )
    return set(json.loads(completed.stdout))


def test_function_analysis_facade_does_not_eagerly_import_property_engines():
    modules = _modules_after(
        "import json, sys; import semialg.function_analysis; print(json.dumps(sorted(sys.modules)))"
    )
    assert "semialg._function_analysis_convexity" not in modules
    assert "semialg._function_analysis_monotonicity" not in modules


def test_function_range_core_does_not_eagerly_import_image_engines():
    modules = _modules_after(
        "import json, sys; import semialg._optimization_range; "
        "print(json.dumps(sorted(sys.modules)))"
    )
    assert "semialg._range_special_cases" not in modules
    assert "semialg._function_graph_image" not in modules


def test_public_calls_load_only_the_required_lazy_engines():
    modules = _modules_after(
        "import json, sys, sympy as sp; "
        "from semialg.function_analysis import function_monotonicity; "
        "x=sp.Symbol('x', real=True); function_monotonicity(x, x); "
        "print(json.dumps(sorted(sys.modules)))"
    )
    assert "semialg._function_analysis_monotonicity" in modules
    assert "semialg._function_analysis_convexity" not in modules
