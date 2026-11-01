"""Operation-order invariants for exact graph decisions and cache resets."""

import sympy as sp
from hypothesis import settings
from hypothesis import strategies as st
from hypothesis.stateful import RuleBasedStateMachine, invariant, rule

from semialg.cache_control import clear_caches
from semialg.decision import equivalent
from semialg.function_verification import verify_real_function

x, y = sp.symbols("x y", real=True)


class GraphCacheMachine(RuleBasedStateMachine):
    def __init__(self):
        super().__init__()
        self.offset = 0

    @rule(offset=st.integers(-2, 2))
    def choose_graph(self, offset):
        self.offset = offset

    @rule()
    def reset(self):
        clear_caches(include_sympy=False, collect=False)

    @invariant()
    def exact_graph(self):
        graph = sp.And(sp.Eq(y * y, x + self.offset), y >= 0)
        result = verify_real_function(graph, input=x, output=y, domain=x >= -self.offset)
        assert result.total and result.single_valued
        assert equivalent(result.totality_condition, x >= -self.offset, (x,))


TestGraphCache = GraphCacheMachine.TestCase
TestGraphCache.settings = settings(max_examples=5, stateful_step_count=4, deadline=None)
