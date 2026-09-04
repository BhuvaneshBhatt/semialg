import pytest

hypothesis = pytest.importorskip("hypothesis")
import sympy as sp
from hypothesis import given, settings
from hypothesis import strategies as st

from semialg import equivalent, is_satisfiable, is_tautology

x = sp.Symbol("x", real=True)


@st.composite
def small_relations(draw):
    a = draw(st.integers(-3, 3).filter(lambda n: n != 0))
    b = draw(st.integers(-5, 5))
    c = draw(st.integers(-5, 5))
    degree = draw(st.integers(1, 2))
    poly = a * x**degree + b * x + c
    op = draw(st.sampled_from(["<", "<=", "==", ">=", ">"]))
    return {"<": poly < 0, "<=": poly <= 0, "==": sp.Eq(poly, 0), ">=": poly >= 0, ">": poly > 0}[
        op
    ]


@settings(max_examples=30, deadline=None)
@given(small_relations())
def test_boolean_identities_hold_for_small_semialgebraic_atoms(atom):
    assert equivalent(sp.And(atom, sp.true), atom, [x])
    assert equivalent(sp.Or(atom, sp.false), atom, [x])
    assert not is_satisfiable(sp.And(atom, sp.Not(atom)), [x])
    assert is_tautology(sp.Or(atom, sp.Not(atom)), [x])


@settings(max_examples=25, deadline=None)
@given(small_relations())
def test_variable_renaming_preserves_satisfiability(atom):
    y = sp.Symbol("y", real=True)
    renamed = atom.xreplace({x: y})
    assert is_satisfiable(atom, [x]) == is_satisfiable(renamed, [y])
