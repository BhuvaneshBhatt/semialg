import sympy as sp

from semialg import (
    classify_real_roots,
    find_instance,
    path_between,
    semialgebraic_maximize,
    semialgebraic_minimize,
)
from semialg.decomposition import component_instances, generic_cad
from semialg.formula import parse_formula
from semialg.geometry_queries import CADPathResult
from semialg.optimization_results import OptimizationResult
from semialg.qe import CompleteQEResult, qe_by_complete_cad
from semialg.root_classification import RootClassificationResult
from semialg.solve.find_instance import InstanceResult
from semialg.solve.zero_dimensional import ZeroDimensionalSolveResult, solve_zero_dimensional_system


def test_optimization_defaults_to_value_and_optimizer_points():
    x = sp.Symbol("x", real=True)
    assert semialgebraic_minimize((x - 2) ** 2, variables=[x]) == [0, [{x: 2}]]
    assert semialgebraic_maximize(x, [x >= 0, x <= 2], [x]) == [2, [{x: 2}]]
    assert isinstance(
        semialgebraic_minimize((x - 2) ** 2, variables=[x], return_result=True),
        OptimizationResult,
    )


def test_instance_defaults_to_instances():
    x = sp.Symbol("x", real=True)
    one = find_instance(sp.Eq(x, 3), [x])
    assert one == {x: 3}
    many = find_instance(sp.Or(sp.Eq(x, 1), sp.Eq(x, 2)), [x], count=2)
    assert isinstance(many, tuple)
    assert set(point[x] for point in many) == {1, 2}
    assert isinstance(find_instance(sp.Eq(x, 3), [x], return_result=True), InstanceResult)


def test_zero_dimensional_solver_defaults_to_exact_points():
    x = sp.Symbol("x", real=True)
    assert solve_zero_dimensional_system([x**2 - 1], variables=[x]) == ((-1,), (1,))
    assert isinstance(
        solve_zero_dimensional_system([x**2 - 1], variables=[x], return_result=True),
        ZeroDimensionalSolveResult,
    )


def test_component_and_generic_cad_defaults_are_mathematical_outputs():
    x = sp.Symbol("x", real=True)
    instances = component_instances(sp.Or(x < -1, x > 1), [x])
    assert isinstance(instances, tuple)
    assert len(instances) == 2
    assert generic_cad(x**2 < 1, [x]) == (x**2 - 1 < 0)


def test_complete_qe_defaults_to_formula():
    x = sp.Symbol("x", real=True)
    formula = qe_by_complete_cad((x,), (("exists", x),), parse_formula(sp.Eq(x**2, 1)))
    assert formula is sp.true
    structured = qe_by_complete_cad(
        (x,), (("exists", x),), parse_formula(sp.Eq(x**2, 1)), return_result=True
    )
    assert isinstance(structured, CompleteQEResult)
    assert structured.formula is sp.true


def test_structured_objects_that_are_the_mathematical_answer_remain_structured():
    x, a = sp.symbols("x a", real=True)
    classification = classify_real_roots(x**2 - a, x, parameters=[a])
    assert isinstance(classification, RootClassificationResult)

    path = path_between(sp.And(x >= 0, x <= 2), {x: 0}, {x: 2}, [x])
    assert isinstance(path, CADPathResult)
    assert path.connected is True
