import sympy as sp

from semialg._strategy_outcome import StrategyOutcome, StrategyStatus
from semialg.reconstruct.root_functions import root_of
from semialg.topology.sample_context import CylindricalSampleContext


def test_strategy_outcome_distinguishes_unknown_from_not_applicable():
    assert StrategyOutcome.not_applicable("shape").status is StrategyStatus.NOT_APPLICABLE
    assert StrategyOutcome.unknown("certificate").status is StrategyStatus.UNKNOWN
    assert StrategyOutcome.success(3).value == 3


def test_cylindrical_sample_context_closes_only_lower_prefix():
    x, y = sp.symbols("x y", real=True)
    ctx = CylindricalSampleContext((x, y), (sp.Integer(2), root_of(y**2 - x, y, 1)))
    assert ctx.assignment_before(y) == {x: 2}
    full = ctx.full_assignment()
    assert full[x] == 2
    assert sp.simplify(full[y] ** 2 - 2) == 0
