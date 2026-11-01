"""Coverage, domain poles and independent representative fibers."""

import os
import subprocess
import sys
from pathlib import Path

import sympy as sp

from semialg.decision import equivalent, is_satisfiable
from semialg.function_verification import verify_real_function
from semialg.qe import quantifier_eliminate

x, y, a = sp.symbols("x y a", real=True)


def test_rational_coefficient_domain():
    graph = sp.And(sp.Ne(a, 0), sp.Eq(a * y * y, x), y >= 0)
    result = verify_real_function(
        graph, input=x, output=y, parameters=(a,), domain=sp.And(a > 0, x >= 0)
    )
    assert result.total and result.single_valued
    assert not is_satisfiable(result.formula.subs(a, 0), (x, y))
    projection = quantifier_eliminate(
        sp.And(sp.Ne(a, 0), sp.Eq(a * y * y, 1)), [("exists", y)], variables=(a, y)
    )
    assert equivalent(projection, a > 0, (a,))


def test_fresh_process_replay(tmp_path):
    code = """import sympy as sp
from semialg.algebraic import certify_polynomial_root_interval,verify_root_interval_certificate
from semialg.cache_control import clear_caches
x=sp.Symbol('x')
c=certify_polynomial_root_interval(x*x-2,0,2,var=x)
clear_caches(include_sympy=False,collect=False)
assert verify_root_interval_certificate(c,polynomial=x*x-2,var=x,left=0,right=2,include_left=True,include_right=True)
"""
    bootstrap = (
        "import sys;sys.path[:]= "
        + repr([str(Path(p).resolve()) for p in sys.path if p])
        + ";"
        + code
    )
    result = subprocess.run(
        [sys.executable, "-S", "-B", "-c", bootstrap],
        cwd=tmp_path,
        env={**os.environ, "PYTHONDONTWRITEBYTECODE": "1"},
        capture_output=True,
        text=True,
        timeout=60,
        check=False,
    )
    assert result.returncode == 0, result.stdout + result.stderr
