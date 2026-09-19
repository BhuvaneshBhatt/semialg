import sympy as sp

from semialg.benchmarks import (
    BenchmarkDomain,
    all_cases,
    cad_cases,
    polynomial_data,
    qe_cases,
    run_groebner,
    topology_cases,
)


def test_catalog_metadata():
    cases = all_cases()
    assert len(cases) >= 10
    assert len({case.name for case in cases}) == len(cases)
    assert all(case.difficulty is not None for case in cases)
    assert all(case.variables for case in cases)
    assert {case.domain for case in cases} >= {
        BenchmarkDomain.CAD,
        BenchmarkDomain.QE,
        BenchmarkDomain.ROADMAP,
        BenchmarkDomain.COMPONENTS,
        BenchmarkDomain.GROEBNER,
    }


def test_published_cad_provenance():
    assert all(case.provenance is not None for case in cad_cases())
    assert any(case.name == "davenport_heintz" for case in cad_cases())


def test_qe_formulas_parse():
    assert all(case.parse() is not None for case in qe_cases())


def test_topology_expectations():
    expected = {case.name: case.expected.get("components") for case in topology_cases()}
    assert expected == {"circle": 1, "two_circles": 2, "isolated_plus_circle": 2}


def test_cyclic_groebner_conformance():
    case = next(case for case in all_cases() if case.name == "groebner_cyclic_3")
    variables, polynomials = polynomial_data(case)
    assert variables == sp.symbols("x1 x2 x3", real=True)
    assert len(polynomials) == 3
    result = run_groebner(case)
    assert result.correct is True
    assert result.metrics.counters["basis_polynomials"] > 0
    assert result.metrics.cpu_seconds >= 0
