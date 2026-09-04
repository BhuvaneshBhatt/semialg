from __future__ import annotations

from types import SimpleNamespace

import sympy as sp

from semialg.cad_algorithms import cell_complex
from semialg.context import SemialgebraicContext


def test_semialgebraic_context_distinguishes_omitted_from_empty_variables() -> None:
    x = sp.Symbol("x", real=True)

    inferred = SemialgebraicContext(x >= 0)
    explicit_empty = SemialgebraicContext(x >= 0, ())

    assert inferred.variables == (x,)
    assert explicit_empty.variables == ()


def test_cell_complex_respects_explicit_empty_variables() -> None:
    x = sp.Symbol("x", real=True)
    region = SimpleNamespace(formula=sp.true, variables=(x,))

    complex_ = cell_complex.build_cad_cell_complex(region, variables=())

    assert complex_.variables == ()
    assert len(complex_.cells) == 1


def test_cell_complex_omitted_variables_still_uses_region_variables() -> None:
    x = sp.Symbol("x", real=True)
    region = SimpleNamespace(formula=sp.true, variables=(x,))

    complex_ = cell_complex.build_cad_cell_complex(region)

    assert complex_.variables == (x,)


def test_symbolic_region_distinguishes_omitted_from_empty_variables() -> None:
    from semialg.symbolic_regions import SemialgebraicRegion

    x = sp.Symbol("x", real=True)

    inferred = SemialgebraicRegion(x >= 0)
    explicit_empty = SemialgebraicRegion(x >= 0, ())

    assert inferred.variables == (x,)
    assert explicit_empty.variables == ()
