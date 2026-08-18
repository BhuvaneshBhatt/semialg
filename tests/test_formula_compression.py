from __future__ import annotations

import pytest
import sympy as sp

from semialg import cad
from semialg.reconstruct.merge import compressed_formula_from_cells, dnf_formula_from_cells
from semialg.reconstruct.nested import nested_formula_from_cells

pytestmark = pytest.mark.slow


def test_nested_reconstruction_reduces_repeated_prefixes_for_ball():
    x, y, z = sp.symbols("x y z", real=True)
    result = cad(x**2 + y**2 + z**2 < 1, [x, y, z], return_result=True)
    nested = nested_formula_from_cells(result.cells, result.variables, result.cad.cells_by_level)
    dnf = dnf_formula_from_cells(result.cells, result.variables, result.cad.cells_by_level)
    assert nested.formula != sp.false
    assert nested.stats.emitted_blocks <= max(1, len(result.cells))
    assert len(sp.sstr(nested.formula)) <= len(sp.sstr(dnf))


def test_cad_accepts_dnf_and_nested_formula_forms():
    x, y, z = sp.symbols("x y z", real=True)
    expr = x**2 + y**2 + z**2 < 1
    nested = cad(expr, [x, y, z], formula_form="nested", return_result=True)
    dnf = cad(expr, [x, y, z], formula_form="dnf", return_result=True)
    assert nested.formula != sp.false
    assert dnf.formula != sp.false
    assert nested.diagnostics["formula_form"] == "nested"
    assert dnf.diagnostics["formula_form"] == "dnf"


def test_compressed_formula_has_size_guard_fallback():
    x, y = sp.symbols("x y", real=True)
    result = cad((x**2 + y**2 < 1) | (x > 2), [x, y], return_result=True)
    compressed = compressed_formula_from_cells(
        result.cells,
        result.variables,
        result.cad.cells_by_level,
        max_terms=0,
    )
    assert compressed.formula != sp.false
    assert compressed.stats.fallback_used
