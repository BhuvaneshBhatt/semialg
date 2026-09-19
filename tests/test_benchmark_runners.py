import pytest

from semialg.benchmarks import cad_cases, groebner_family, run_cad, run_groebner, tticad_cases


def test_cad_runner_counts_cells():
    case = next(c for c in cad_cases() if c.provenance.identifier == "CMXYExamples:1")
    result = run_cad(case)
    assert result.correct is True
    assert result.metrics.counters["cells"] > 0
    assert result.metrics.counters["projection_poly_count_by_level"]


@pytest.mark.parametrize("family,n", [("cyclic", 3), ("katsura", 3), ("eco", 3), ("noon", 3)])
def test_groebner_families(family, n):
    case = groebner_family(family, n)
    result = run_groebner(case)
    assert result.correct is True
    assert result.metrics.counters["basis_polynomials"] > 0


def test_tticad_reference_counts_are_configuration_specific():
    case = tticad_cases()[0]
    assert case.expected["reference_cell_counts"] == {
        "published_full_cad": 3707,
        "published_tticad": 269,
    }
