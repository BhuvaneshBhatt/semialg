import sympy as sp

from semialg.cad_algorithms.cells import CylindricalSolution, StructuredCADCell
from semialg.exact_arithmetic import compare_extended_reals
from semialg.instances.real_fallbacks import CoordinateBounds
from semialg.normalization import normalize_parameters
from semialg.reasoning import SignProofResult, SimplifiedSystem
from semialg.solve.integer._linear_candidates import linear_equality_candidate


def test_shared_parameter_normalization_preserves_assumption_identity():
    a = sp.Symbol("a", positive=True)
    assert normalize_parameters(["a"], a + 1) == (a,)


def test_extended_real_comparison_covers_infinities_and_finite_values():
    assert compare_extended_reals(-sp.oo, sp.Integer(0)) == -1
    assert compare_extended_reals(sp.Integer(0), sp.oo) == -1
    assert compare_extended_reals(sp.oo, -sp.oo) == 1
    assert compare_extended_reals(sp.sqrt(2), sp.Rational(3, 2)) == -1


def test_public_dataclasses_report_defining_modules():
    assert StructuredCADCell.__module__ == "semialg.cad_algorithms.structured_cells"
    assert CylindricalSolution.__module__ == "semialg.cad_algorithms.cylindrical_solution"
    assert CoordinateBounds.__module__ == "semialg.instances.real_utils"
    assert SimplifiedSystem.__module__ == "semialg.reasoning_results"
    assert SignProofResult.__module__ == "semialg.reasoning_results"


def test_shared_integer_linear_candidate_matches_expected_elimination_data():
    x, y = sp.symbols("x y", integer=True)
    data = linear_equality_candidate(sp.expand(3 * x + 2 * y - 7), x)
    assert data is not None
    coefficient, numerator, denominator, replacement = data
    assert coefficient == 3
    assert numerator == 7 - 2 * y
    assert denominator == 3
    assert sp.simplify(replacement - (7 - 2 * y) / 3) == 0
