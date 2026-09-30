"""Documentation smoke tests against an installed wheel, not the source tree."""

from __future__ import annotations

import os
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]


@pytest.mark.slow
def test_gallery_and_readme_public_examples_from_installed_wheel(tmp_path):
    dist = tmp_path / "dist"
    dist.mkdir()
    build_code = f"import setuptools.build_meta as backend; backend.build_wheel({str(dist)!r})"
    subprocess.run(
        [sys.executable, "-c", build_code],
        cwd=ROOT,
        check=True,
        capture_output=True,
        text=True,
    )
    wheel = next(dist.glob("semialg-*.whl"))
    target = tmp_path / "site"
    subprocess.run(
        [sys.executable, "-m", "pip", "install", "--no-deps", "--target", str(target), str(wheel)],
        check=True,
        capture_output=True,
        text=True,
    )
    env = os.environ.copy()
    env["PYTHONPATH"] = str(target)
    code = """
import sympy as sp
from semialg import find_instance, local_germ, local_sign_strata, parameter_strata
x, a = sp.symbols("x a", real=True)
germ = local_germ(sp.Ne(x, 0), (0,), (x,))
assert {s.signs for s in local_sign_strata((x,), germ.formula, germ.point, germ.variables)} == {(-1,), (1,)}
p = sp.Symbol("p", positive=True)
assert find_instance(sp.true, (p,))[p] > 0
strata = parameter_strata(sp.Eq(x**2, a), (x,), (a,))
assert len(strata) == 2
"""
    subprocess.run(
        [sys.executable, "-c", code],
        cwd=tmp_path,
        env=env,
        check=True,
        capture_output=True,
        text=True,
    )
