"""Specialization changes multiplicity while leaving distinct-root ranks intact."""

import pickle
from unittest.mock import PropertyMock, patch

import pytest
import sympy as sp

from semialg.algebraic import isolate_real_roots, root_multiplicity
from semialg.algebraic.roots import refine_root_once
from semialg.algebraic.samples import AlgebraicRoot, FiberRootContext, RationalSample
from semialg.cad_algorithms.decomposition import _roots_in_parent_tower
from semialg.errors import ExactComputationUnavailable


def _fiber(power):
    x, y = sp.symbols("x y", real=True)
    parent = isolate_real_roots(sp.Poly(x * x - 2, x))[1]
    context = FiberRootContext(sp.Poly(y**power - x * x + 2, x, y), (x, y), (parent,), y)
    return _roots_in_parent_tower(context)[0]


@pytest.mark.parametrize("power", [2, 3, 4])
def test_specialized_multiplicity(power):
    root = _fiber(power)
    assert root.multiplicity_hint is None
    assert root.root_index == 0 and root.as_expr() == 0
    assert root.multiplicity == root_multiplicity(root) == power
    assert pickle.loads(pickle.dumps(root)).multiplicity == power


def test_refinement_does_not_force_multiplicity():
    x, y = sp.symbols("x y", real=True)
    parent = isolate_real_roots(sp.Poly(x * x - 2, x))[1]
    context = FiberRootContext(sp.Poly(y * y - x, x, y), (x, y), (parent,), y)
    root = _roots_in_parent_tower(context)[0]
    with patch.object(
        AlgebraicRoot,
        "multiplicity",
        new_callable=PropertyMock,
        side_effect=AssertionError("eager derivative work"),
    ):
        refined = refine_root_once(root)
    assert refined.interval.left >= root.interval.left
    assert refined.interval.right <= root.interval.right
    assert root_multiplicity(refined) == 1


def test_nested_parent_multiplicity():
    x, y, z = sp.symbols("x y z", real=True)
    parent = isolate_real_roots(sp.Poly(x * x - 2, x))[1]
    second = _roots_in_parent_tower(
        FiberRootContext(sp.Poly(y * y - x, x, y), (x, y), (parent,), y)
    )[1]
    context = FiberRootContext(sp.Poly(z**3 - y * y + x, x, y, z), (x, y, z), (parent, second), z)
    root = _roots_in_parent_tower(context)[0]
    assert root.multiplicity == 3
    assert root.as_expr() == 0


def test_nullified_fiber_has_no_finite_multiplicity():
    x, y = sp.symbols("x y", real=True)
    from semialg.algebraic.intervals import RationalInterval

    context = FiberRootContext(sp.Poly(x * y, x, y), (x, y), (RationalSample(sp.S.Zero),), y)
    root = AlgebraicRoot(sp.Poly(y, y), RationalInterval(0, 0), 0, None, None, context)
    with pytest.raises(ExactComputationUnavailable, match="nullified"):
        root_multiplicity(root)
