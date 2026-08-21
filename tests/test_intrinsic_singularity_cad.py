from __future__ import annotations

import pytest
import sympy as sp

from semialg import (
    extract_cylindrical_solution,
    intrinsic_solution_integrals,
    stratify_intrinsic_solution,
)

pytestmark = pytest.mark.slow


def test_cusp_is_explicitly_stratified_and_regular_branches_integrate():
    x, y = sp.symbols("x y", real=True)
    sol = extract_cylindrical_solution(sp.And(x >= 0, x <= 1, sp.Eq(y**2, x**3)), [x, y])
    strat = stratify_intrinsic_solution(sol)
    singular = [s for s in strat.singular_strata if s.cell.sample[x] == 0 and s.cell.sample[y] == 0]
    assert len(singular) == 1
    assert singular[0].dimension == 0
    regular_curves = [s for s in strat.regular_strata if s.dimension == 1]
    assert len(regular_curves) == 2

    pieces = intrinsic_solution_integrals(sol, 1, dimension=1, evaluate=True)
    assert len(pieces) == 2
    total = sp.simplify(sum(piece.integral for piece in pieces))
    assert sp.simplify(total - (26 * sp.sqrt(13) - 16) / 27) == 0
