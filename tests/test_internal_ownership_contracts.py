from pathlib import Path


def test_cylindrical_sample_specialization_has_one_owner():
    root = Path(__file__).parents[1] / "src" / "semialg"
    incidence = (root / "topology" / "incidence.py").read_text()
    assert "CylindricalSampleContext" in incidence
    assert (
        "def _closed_assignment_items" in incidence
    )  # retained cache API, not a second policy owner


def test_function_graph_documents_exact_rational_power_boundary():
    text = (Path(__file__).parents[1] / "src" / "semialg" / "function_graph.py").read_text()
    assert "principal rational power" in text
    assert "explicit real-root" in text
