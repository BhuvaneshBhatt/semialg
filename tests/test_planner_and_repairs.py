import sympy as sp

from semialg.cad.reduced import decompose_reduced_safe
from semialg.formula import parse_formula
from semialg.planner import brown_variable_order, cand_var_orders, extract_problem_features
from semialg.simplify.implication import minimize_disj_by_impl
from semialg.simplify.result import simp_semialg_expr
from semialg.solve.reduce import reduce_text


def test_planner_and_01():
    x, y = sp.symbols("x y", real=True)
    expr = sp.Or(sp.And(x > 1, x > 0, y > 0), sp.And(x > 1, y > 0))
    assert simp_semialg_expr(expr) == sp.And(x > 1, y > 0)
    direct = minimize_disj_by_impl(expr)
    assert direct == sp.And(x > 1, y > 0)


def test_planner_and_02():
    x, y = sp.symbols("x y", real=True)
    polys = (y**2 - x, x**2 + y)
    formula = parse_formula((y**2 - x >= 0) & (x**2 + y >= 0))
    features = extract_problem_features(formula, variables=(y, x))
    orders = cand_var_orders(features, polys)
    assert orders
    assert brown_variable_order(polys, (y, x)) in {item.order for item in orders}
    assert orders[0].score <= orders[-1].score


def test_planner_and_03():
    solved = reduce_text("exists y. y^2 - x = 0", strategy="auto", return_result=True)
    assert solved.result == sp.Ge(sp.Symbol("x", real=True), 0)
    assert "strategy_selection" in solved.metadata


def test_planner_and_04():
    x, y = sp.symbols("x y", real=True)
    result = decompose_reduced_safe(
        [y**2 - x, y - 1],
        (x, y),
        backend="mccallum",
        equational_constraints=(y**2 - x,),
    )
    assert result.complete
    assert result.effective_backend in {
        "mccallum-reduced-certified",
        "mccallum-repaired-certified",
        "collins-complete",
    }
    if result.effective_backend == "mccallum-repaired-certified":
        assert result.reduced_projection.tower.metadata.get("delineating_repair") is True
