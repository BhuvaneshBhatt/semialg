import os
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]


def test_simplex_dimension_fast_path_does_not_import_cad_decomposition():
    code = r"""
import sys
from pathlib import Path
from semialg.standard_regions import Simplex
from semialg.regions.operations import region_dimension
region = Simplex(((0, 0), (1, 0), (0, 1)))
assert "semialg.decomposition" not in sys.modules
assert region_dimension(region) == 2
assert "semialg.decomposition" not in sys.modules
"""
    env = dict(os.environ)
    env["PYTHONPATH"] = str(ROOT / "src")
    completed = subprocess.run(
        [sys.executable, "-c", code],
        cwd=ROOT,
        env=env,
        capture_output=True,
        text=True,
        check=False,
    )
    assert completed.returncode == 0, completed.stderr
