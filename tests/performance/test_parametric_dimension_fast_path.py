import subprocess
import sys


def test_simplex_dimension_fast_path_does_not_import_cad_decomposition():
    code = r"""
import sys
from semialg.standard_regions import SimplexRegion
from semialg.regions.operations import region_dimension
region = SimplexRegion(((0, 0), (1, 0), (0, 1)))
assert "semialg.decomposition" not in sys.modules
assert region_dimension(region) == 2
assert "semialg.decomposition" not in sys.modules
"""
    completed = subprocess.run(
        [sys.executable, "-c", code], capture_output=True, text=True, check=False
    )
    assert completed.returncode == 0, completed.stderr
