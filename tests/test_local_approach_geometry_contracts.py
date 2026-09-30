import pytest
import sympy as sp

import semialg
from semialg.geometry_queries import semialgebraic_image


def test_correlated_image_rejects_wrong_target_arity_and_source_collision():
    x, y = sp.symbols("x y", real=True)
    with pytest.raises(ValueError, match="same length"):
        semialg.correlated_map_image((x, x**2), variables=(x,), image_variables=(y,))
    with pytest.raises(ValueError, match="distinct from source"):
        semialg.correlated_map_image((x,), variables=(x,), image_variables=(x,))


def test_semialgebraic_image_rejects_parameter_target_collision():
    x, p = sp.symbols("x p", real=True)
    with pytest.raises(ValueError, match="distinct from parameters"):
        semialgebraic_image(x + p, sp.true, (x,), parameters=(p,), image_variables=(p,))


def test_projection_rejects_duplicate_elimination_variables():
    x, y = sp.symbols("x y", real=True)
    with pytest.raises(ValueError, match="must be distinct"):
        semialg.semialgebraic_projection(sp.Eq(x + y, 0), (x, x), (x, y))


def test_structured_diagnostics_do_not_truth_test_symbolic_metadata():
    x = sp.Symbol("x", real=True)

    class Result:
        method = "symbolic"
        diagnostics = {"critical_values": sp.Gt(x, 0)}
        proof_trace = ()

    diagnostics = semialg.structured_proof_diagnostics(Result())
    assert diagnostics.used_critical_values is True
    assert "critical_value_geometry" in diagnostics.steps


def test_correlated_image_accepts_python_scalar_mapping():
    x, y = sp.symbols("x y", real=True)
    result = semialg.correlated_map_image(1, sp.Eq(x, 0), (x,), image_variables=(y,))
    assert sp.simplify(result.formula.subs(y, 1)) is sp.true
    assert sp.simplify(result.formula.subs(y, 2)) is sp.false


def test_structured_diagnostics_recognize_critical_value_set_metadata():
    x = sp.Symbol("x", real=True)

    class Result:
        method = "symbolic"
        diagnostics = {"critical_value_set": sp.Eq(x, 0)}
        proof_trace = ()

    diagnostics = semialg.structured_proof_diagnostics(Result())
    assert diagnostics.used_critical_values is True
    assert "critical_value_geometry" in diagnostics.steps
