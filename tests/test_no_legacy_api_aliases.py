"""The 1.2 API intentionally carries no compatibility aliases."""

import semialg
import semialg.algebraic as algebraic
import semialg.reconstruct as reconstruct
import semialg.simplify as simplify
from semialg.planner import select


def test_only_canonical_names_are_exposed() -> None:
    for name in (
        "RootFunction",
        "compare_alg_numbers",
        "sort_algebraic_numbers",
        "canonicalize_qe_formula",
    ):
        assert name not in semialg.__all__
        assert not hasattr(semialg, name)

    assert not hasattr(reconstruct, "RootFunction")
    assert not hasattr(algebraic, "compare_alg_numbers")
    assert not hasattr(algebraic, "sort_algebraic_numbers")
    assert not hasattr(simplify, "canonicalize_qe_formula")
    assert not hasattr(select, "PROJECTION_EC")
