"""Replay must not depend on construction caches or accept another statement."""

import pickle

import pytest
import sympy as sp

from semialg.algebraic.samples import RationalSample
from semialg.cache_control import clear_caches
from semialg.cache_utils import cache_disabled
from semialg.cad_algorithms import decomp_lazard_complete, verify_cad
from semialg.topology.branch_correspondence import (
    certify_branch_correspondence,
    replay_branch_correspondence,
)


def test_native_cad_context():
    x, y = sp.symbols("x y", real=True)
    polynomial = y * y - x
    original = decomp_lazard_complete((polynomial,), (x, y))
    _ = original._point_locator
    source = pickle.loads(pickle.dumps(original))
    assert "_point_locator" not in source.__dict__
    with pytest.raises(TypeError):
        source.cells_by_level[1] = ()
    with pytest.raises(TypeError):
        source.tower.metadata["complete"] = False
    clear_caches(collect=False)
    with cache_disabled():
        assert verify_cad(source, polynomials=(polynomial,))
        assert not verify_cad(source, polynomials=(y * y + x,))


def test_branch_context():
    x, y = sp.symbols("x y", real=True)
    proof = certify_branch_correspondence(
        ((y + 1) ** 2 + x) * (y - 2), y, 0, (x,), (RationalSample(sp.S.One),), (sp.S.Zero,)
    )
    restored = pickle.loads(pickle.dumps(proof))
    clear_caches(collect=False)
    with cache_disabled():
        assert replay_branch_correspondence(restored)
    assert restored.value == 2
    assert sp.expand(restored.polynomial.as_expr() - ((y + 1) ** 2 + x) * (y - 2)) == 0


@pytest.mark.parametrize("output", ["cells", "function", "tree"])
def test_public_cad_pickle(output):
    from semialg import cad, replay_certificate

    x, y = sp.symbols("x y", real=True)
    result = cad(sp.Eq(y * y, x), (x, y), output=output, return_result=True)
    restored = pickle.loads(pickle.dumps(result))
    assert replay_certificate(restored).verified is True
    clear_caches(collect=False)
    function = restored.as_function()
    assert function(1, 1) is True
    assert function(-1, 0) is False
    assert function(1, 0) is False
    assert verify_cad(restored.cad, polynomials=(y * y - x,))
    with pytest.raises(TypeError):
        restored.cell_set.cells_by_level[1] = ()
