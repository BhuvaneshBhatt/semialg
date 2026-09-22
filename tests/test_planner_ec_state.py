import sympy as sp

from semialg.cad_algorithms.constants import PROJECTION_MCCALLUM
from semialg.planner.analyze import ProblemAnalysis
from semialg.planner.features import ProblemFeatures
from semialg.planner.select import select_strat_analysis


def test_equational_constraint_is_projection_state_not_backend_alias() -> None:
    x, y = sp.symbols("x y")
    features = ProblemFeatures(
        variables=(x, y),
        num_polynomials=2,
        num_atoms=2,
        equality_count=1,
        inequality_count=1,
        has_ecs=True,
        suggested_variable_order=(x, y),
    )
    selected = select_strat_analysis(ProblemAnalysis(features))
    assert selected.backend == PROJECTION_MCCALLUM
    assert selected.projection.operator == PROJECTION_MCCALLUM
    assert selected.projection.use_ecs is True
