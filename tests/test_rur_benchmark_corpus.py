from semialg.benchmarks.rur_corpus import rur_benchmark_cases


def test_rur_corpus_has_unique_named_systems():
    cases = rur_benchmark_cases()
    assert len(cases) >= 5
    assert len({case.name for case in cases}) == len(cases)
    assert all(case.equations and case.variables for case in cases)
