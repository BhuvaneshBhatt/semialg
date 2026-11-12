"""Check algorithm dispatch and exact real-answer validation in the benchmark."""

import importlib.util
from pathlib import Path
from types import SimpleNamespace

import pytest
import sympy as sp
from sympy.polys.polyconfig import query, setup

MODULE_PATH = Path(__file__).resolve().parents[1] / "scripts/benchmark_field_pipeline.py"
spec = importlib.util.spec_from_file_location("field_pipeline_benchmark", MODULE_PATH)
benchmark = importlib.util.module_from_spec(spec)
spec.loader.exec_module(benchmark)


@pytest.mark.parametrize("method", ["buchberger", "f5b"])
def test_method_dispatch(method):
    x, y = sp.symbols("x y")
    previous = query("groebner")
    with benchmark.sympy_method(method), benchmark.audit_sympy_methods() as calls:
        basis = sp.groebner((x * y - 1, y**2 - x), x, y)
        assert basis.reduce(x * y - 1)[1] == 0
    assert calls == {"buchberger": int(method == "buchberger"), "f5b": int(method == "f5b")}
    assert query("groebner") == previous


def test_method_restoration():
    previous = query("groebner")
    setup("groebner", "f5b")
    try:
        with pytest.raises(RuntimeError, match="interrupted"):
            with benchmark.sympy_method("buchberger"):
                raise RuntimeError("interrupted")
        assert query("groebner") == "f5b"
    finally:
        setup("groebner", previous)


def test_missing_exceptional_point():
    t = sp.Symbol("t", real=True)
    answer = SimpleNamespace(formula=sp.Ne(t, sp.Rational(13, 37)), diagnostics={})
    with pytest.raises(ArithmeticError, match="wrong parameter condition"):
        benchmark.verify_solvability(t, sp.true, answer)


def test_strict_endpoint():
    t = sp.Symbol("t", real=True)
    expected = sp.And(t > 0, t < 1)
    answer = SimpleNamespace(formula=sp.And(t >= 0, t < 1), diagnostics={})
    with pytest.raises(ArithmeticError, match="wrong parameter condition"):
        benchmark.verify_solvability(t, expected, answer)


def test_f5b_record():
    cases = benchmark.workloads("sympy", None)
    row = benchmark.run_case("rur_rational", cases["rur_rational"], "sympy", "f5b", None, 1)
    assert row["construction_method"] == row["sympy_method"] == "f5b"
    assert row["check"]["dimension"] == row["check"]["geometric_points"] == 4
    assert row["method_audit"]["sympy_kernel_calls"]["f5b"] > 0
    assert row["method_audit"]["sympy_kernel_calls"]["buchberger"] == 0
    assert len(row["timing"]["samples_seconds"]) == len(row["groebner_profiles"]) == 1


@pytest.mark.parametrize(
    "name,source,variables,expected",
    benchmark.real_solvability_cases(),
    ids=[case[0] for case in benchmark.real_solvability_cases()],
)
def test_matched_real_cases(name, source, variables, expected):
    t = sp.Symbol("t", real=True)
    with benchmark.sympy_method("f5b"):
        result, records, stages = benchmark.solvability_operation(
            source, variables, (t,), "sympy", None, matched=True
        )
    checked = benchmark.verify_solvability(t, expected, result)
    assert checked["equivalent_on_all_real_parameters"]
    assert result.diagnostics["groebner_preprocessing"]
    assert set(stages) == {"equality_preprocessing", "real_elimination"}
    assert {record["backend"] for record in records} == {"sympy"}
