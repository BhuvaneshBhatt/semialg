import subprocess
import sys


def test_affine_box_clip_stays_independent_of_cad_imports():
    code = r"""
import sys
from semialg.polyhedral_clipping import clip_affine_subspace_to_box
assert "semialg.decomposition" not in sys.modules
result = clip_affine_subspace_to_box(
    (0, 0, 0), ((1, 0, 0), (0, 1, 0)), ((-1, 1), (-2, 2), (-3, 3))
)
assert len(result.vertices) == 4
assert "semialg.decomposition" not in sys.modules
"""
    completed = subprocess.run(
        [sys.executable, "-c", code], capture_output=True, text=True, check=False
    )
    assert completed.returncode == 0, completed.stderr
