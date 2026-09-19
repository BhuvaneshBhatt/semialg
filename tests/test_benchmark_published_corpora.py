from semialg.benchmarks.cases import cad_cases, tticad_cases


def test_wilson_bank_complete():
    cases = cad_cases()
    assert len(cases) == 68
    ids = {case.expected["source_identifier"] for case in cases}
    assert "CMXYExamples:1" in ids
    assert "BHExamples:10" in ids
    assert "OtherExamples:14" in ids
    assert "JoukowskyTransformation:3" in ids
    assert "TTICADExamples:18" in ids
    assert all(case.variables == case.expected["variable_order"] for case in cases)


def test_tticad_section82_complete():
    cases = tticad_cases()
    assert len(cases) == 29
    assert cases[0].provenance.identifier.startswith("Section82:1:")
    assert cases[-1].provenance.identifier.startswith("Section82:29:")
    assert cases[4].expected["formula_count"] == 2
    assert cases[14].expected["formula_count"] == 7
