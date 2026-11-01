"""Projection and lifting replay checks mathematical data, not success flags."""

from dataclasses import replace

import pytest
import sympy as sp

from semialg.cad_algorithms import decomp_collins_complete, decomp_lazard_complete, verify_cad
from semialg.cad_algorithms.reduced import decomp_form_reduced_safe
from semialg.formula import parse_formula


@pytest.mark.parametrize("backend", [decomp_collins_complete, decomp_lazard_complete])
def test_ambient_geometry_replay(backend):
    x, y = sp.symbols("x y", real=True)
    source = backend((y**2 - x, x - 1), (x, y))
    result = verify_cad(source, polynomials=(y**2 - x, x - 1), return_result=True)
    assert result.valid, result.failures
    assert result.checked_cells == sum(map(len, source.cells_by_level.values()))
    assert result.checked_root_selectors > 0
    assert verify_cad(replace(source, complete=False), polynomials=(y**2 - x, x - 1))
    assert not verify_cad(source, polynomials=(y**2 + x, x - 1))


@pytest.mark.parametrize(
    "mutation",
    [
        "projection",
        "provenance",
        "rank",
        "selector",
        "interval",
        "bounds",
        "sign",
        "coverage",
        "expression",
        "context",
    ],
)
def test_geometric_mutations(mutation):
    from semialg.algebraic.intervals import RationalInterval
    from semialg.algebraic.samples import RationalSample

    x, y = sp.symbols("x y", real=True)
    source = decomp_lazard_complete((y**2 - x,), (x, y))
    if mutation in {"projection", "provenance"}:
        level = source.tower.level(1)
        if mutation == "projection":
            changed = replace(level, polynomials=())
        else:
            changed = replace(
                level, entries=(replace(level.entries[0], parents=("wrong-parent",)),)
            )
        tower = replace(source.tower, levels=(changed, source.tower.level(2)))
        mutated = replace(source, tower=tower)
    else:
        cells = list(source.cells_by_level[2])
        position = next(i for i, c in enumerate(cells) if c.is_section)
        cell = cells[position]
        if mutation == "rank":
            root = replace(cell.sample[-1], root_index=cell.sample[-1].root_index + 1)
            cells[position] = replace(cell, sample=(*cell.sample[:-1], root))
        elif mutation == "selector":
            cells[position] = replace(cell, section_polynomial=y**2 + x)
        elif mutation == "interval":
            root = replace(cell.sample[-1], interval=RationalInterval(sp.S(20), sp.S(21)))
            cells[position] = replace(cell, sample=(*cell.sample[:-1], root))
        elif mutation == "expression":
            root = replace(cell.sample[-1], root_expr=sp.S(42))
            cells[position] = replace(cell, sample=(*cell.sample[:-1], root))
        elif mutation == "context":
            context = replace(
                cell.sample[-1].fiber_context, parent_samples=(RationalSample(sp.S(42)),)
            )
            root = replace(cell.sample[-1], fiber_context=context)
            cells[position] = replace(cell, sample=(*cell.sample[:-1], root))
        elif mutation == "bounds":
            neighbor = position - 1
            cells[neighbor] = replace(cells[neighbor], interval=(None, None))
        elif mutation == "sign":
            signs = dict(cell.signs)
            signs[next(iter(signs))] = 1
            cells[position] = replace(cell, signs=signs)
        else:
            cells.pop(position)
        mutated = replace(source, cells_by_level={1: source.cells_by_level[1], 2: tuple(cells)})
    assert not verify_cad(mutated)


def test_valuation_and_residue_replay():
    x, y, t, z = sp.symbols("x y t z", real=True)
    source = decomp_lazard_complete((x * z + y,), (x, y, t, z))
    result = verify_cad(source, return_result=True)
    assert result.valid, result.failures
    assert result.checked_valuation_claims > 0
    records = source.tower.metadata["lazard_lifts"]
    for changed in [
        (),
        (replace(records[0], valuation=(9,) * len(records[0].valuation)), *records[1:]),
        (replace(records[0], residue=sp.Poly(z + 1, *records[0].variables)), *records[1:]),
    ]:
        tower = replace(
            source.tower, metadata={**dict(source.tower.metadata), "lazard_lifts": changed}
        )
        assert not verify_cad(replace(source, tower=tower))


def test_recursive_ec_projection_replay():
    x, y, z = sp.symbols("x y z", real=True)
    matrix = parse_formula(sp.And(sp.Eq(y - x, 0), sp.Eq(z**2 - y, 0)))
    source = decomp_form_reduced_safe(matrix, (x, y, z), backend="mccallum").cad
    result = verify_cad(source, formula=matrix, return_result=True)
    assert result.valid, result.failures
    assert result.invariant == "truth"
    assert not verify_cad(source)
    tower = replace(source.tower, metadata={**dict(source.tower.metadata), "complete": True})
    assert not verify_cad(replace(source, tower=tower), formula=matrix)


def test_local_delineating_family_replay():
    x, y, z, w = sp.symbols("x y z w", real=True)
    source = decomp_collins_complete((w**2 + y * z - x**2,), (x, y, z, w))
    result = verify_cad(source, return_result=True)
    assert result.valid, result.failures
    local = [c for c in source.cells_by_level[3] if c.index[:2] == (1, 1) and c.is_section]
    assert len(local) == 1
    assert sp.Poly(local[0].section_polynomial.subs({x: 0, y: 0}), z).monic().as_expr() == z


def test_lazard_field_mutations():
    from semialg.cad_algorithms.lifting.lazard import native_lazard_evaluate

    from ._certificate_mutation import assert_certificate_fields_reject_mutation

    x, y, t, z = sp.symbols("x y t z", real=True)
    source = decomp_lazard_complete((x * z + y,), (x, y, t, z))
    claim = source.tower.metadata["lazard_lifts"][0]
    assert_certificate_fields_reject_mutation(
        claim,
        lambda candidate: (
            candidate
            == native_lazard_evaluate(claim.polynomial, claim.variables, claim.parent_samples)
        ),
    )


def test_zero_coordinate_statement():
    source = decomp_lazard_complete((1,), ())
    assert verify_cad(source, polynomials=(2,), formula=sp.true)
    x = sp.Symbol("x", real=True)
    assert not verify_cad(source, polynomials=(x,))
    assert not verify_cad(source, formula=x > 0)
    changed = replace(source.cells[0], level=1)
    assert not verify_cad(replace(source, cells_by_level={0: (changed,)}))
