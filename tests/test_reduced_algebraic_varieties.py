import sympy as sp

from semialg import (
    certified_radicalization,
    irreducible_components,
    reduced_component_singular_loci,
    singular_locus,
)


def test_certified_radicalization_removes_nilpotent_multiplicity():
    x, y = sp.symbols("x y", real=True)
    reduced = certified_radicalization((x**2,), (x, y))
    assert reduced.equations == (x,)
    assert reduced.formula == sp.Eq(x, 0)


def test_singular_locus_uses_reduced_presentation_by_default():
    x, y = sp.symbols("x y", real=True)
    assert singular_locus((x**2,), (x, y)) is sp.false


def test_reducible_nonradical_variety_has_certified_minimal_primes():
    x, y = sp.symbols("x y", real=True)
    components = irreducible_components((x**2 * y**2,), (x, y))
    assert {component.equations for component in components} == {(x,), (y,)}
    assert all(component.dimension == 1 for component in components)


def test_component_crossing_is_not_component_singularity():
    x, y = sp.symbols("x y", real=True)
    loci = reduced_component_singular_loci((x**2 * y**2,), (x, y))
    assert len(loci) == 2
    assert all(item.formula is sp.false for item in loci)


def test_cusp_remains_singular_after_radicalization():
    x, y = sp.symbols("x y", real=True)
    locus = singular_locus((y**2 - x**3,), (x, y))
    assert bool(locus.subs({x: 0, y: 0}))
    assert not bool(locus.subs({x: 1, y: 1}))
