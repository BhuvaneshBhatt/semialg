import sympy as sp

from semialg.formula import parse_formula
from semialg.qe.complete import qe_by_complete_cad


def test_qe_reconstruction_keeps_original_free_polynomial_constraints():
    x, y = sp.symbols("x y", real=True)
    result = qe_by_complete_cad(
        [x, y],
        [("exists", y)],
        parse_formula(sp.And(sp.Eq(y**2, x), x <= 4)),
        free_variables=[x],
        return_result=True,
    )
    assert sp.simplify_logic(sp.Equivalent(result.formula, sp.And(x >= 0, x <= 4))) is sp.true
    assert not result.formula.has(sp.Function("root_of"))
