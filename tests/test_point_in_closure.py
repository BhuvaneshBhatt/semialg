import sympy as sp

from semialg import PointInClosureResult, point_in_closure


def test_point_in_closure_detects_arbitrarily_near_punctured_region():
    x, y = sp.symbols("x y", real=True)
    punctured_disk = sp.And(x**2 + y**2 < 1, sp.Ne(x**2 + y**2, 0))

    assert point_in_closure(punctured_disk, (0, 0), (x, y)) is True
    assert point_in_closure(punctured_disk, (2, 0), (x, y)) is False


def test_point_in_closure_returns_structured_cad_semantic_evidence():
    x = sp.Symbol("x", real=True)
    result = point_in_closure(x > 0, {x: 0}, return_result=True)

    assert isinstance(result, PointInClosureResult)
    assert result.in_closure is True
    assert bool(result) is True
    assert result.variables == (x,)
    assert result.point == {x: 0}
    assert result.method == "cad-semantic-closure"
    assert bool(result.closure.subs(x, 0))


def test_point_in_closure_respects_explicit_syntactic_strategy():
    x = sp.Symbol("x", real=True)
    result = point_in_closure(x > 0, (0,), (x,), strategy="syntactic", return_result=True)

    assert result.in_closure is True
    assert result.method == "syntactic-closure"


def test_point_in_closure_result_can_be_constructed_as_public_evidence_type():
    x = sp.Symbol("x", real=True)
    result = PointInClosureResult(True, x > 0, {x: 0}, (x,), x >= 0)
    assert result.in_closure and bool(result)
