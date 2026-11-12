"""Verified choices in complete GTZ and RUR pipelines, including cache failures."""

import pytest
import sympy as sp
from algroots import GroebnerLimitError, GroebnerVerificationError, NativeF4Limits
from algroots.native_f4 import compiled_available
from algroots.profiling import GroebnerProfiler

from semialg.algebraic.cache import CACHE
from semialg.algebraic.gtz_primary import (
    clear_gtz_caches,
    gtz_primary_decomposition,
    verify_gtz_primary_decomposition_certificate,
)
from semialg.algebraic.rational_univariate import compute_rur

x, y, t = sp.symbols("x y t")


@pytest.mark.parametrize("backend", ["native"])
@pytest.mark.parametrize("source", [(x * y,), ((x + y) ** 2, y * (x + y)), (x * (x - 1),)])
def test_gtz_backends(backend, source):
    if backend == "native" and not compiled_available():
        pytest.skip("native F4 is not installed")
    clear_gtz_caches()
    with GroebnerProfiler() as profile:
        result = gtz_primary_decomposition(source, (x, y), groebner_backend=backend)
    assert result.complete and result.irredundant
    assert verify_gtz_primary_decomposition_certificate(result.certificate)
    assert len(result.components) == 2
    assert all(row["completed"] for row in profile.records)
    assert any("(" in row["domain"] for row in profile.records)
    assert all(row["backend"] == backend for row in profile.records)


@pytest.mark.skipif(not compiled_available(), reason="native F4 is not installed")
def test_gtz_limit_after_cache():
    source = ((x + y) ** 2, y * (x + y))
    assert gtz_primary_decomposition(source, (x, y), groebner_backend="native").complete
    with pytest.raises(GroebnerLimitError):
        gtz_primary_decomposition(
            source, (x, y), groebner_backend="native", groebner_limits=NativeF4Limits(max_terms=1)
        )
    with pytest.raises(ValueError, match="backend"):
        gtz_primary_decomposition(source, (x, y), groebner_backend="missing")


@pytest.mark.skipif(not compiled_available(), reason="native F4 is not installed")
def test_gtz_verification_failure(monkeypatch):
    import algroots.groebner as module

    clear_gtz_caches()
    original = module.compute_groebner_basis

    def wrong(source, ring, **kwargs):
        if kwargs.get("backend") == "native":
            return (ring.one,)
        return original(source, ring, **kwargs)

    monkeypatch.setattr(module, "compute_groebner_basis", wrong)
    with pytest.raises(GroebnerVerificationError):
        gtz_primary_decomposition((x * y,), (x, y), groebner_backend="native")


@pytest.mark.parametrize("backend", ["sympy", "native"])
def test_rur_domain_cache(backend):
    if backend == "native" and not compiled_available():
        pytest.skip("native F4 is not installed")
    CACHE.rur.clear()
    first = compute_rur(
        [x * x - t, y - x], (x, y), domain=sp.QQ.frac_field(t), groebner_backend=backend
    )
    assert first.dimension == 2 and first.solution_count == 2
    with GroebnerProfiler() as profile:
        cached = compute_rur(
            [x * x - t, y - x], (x, y), domain=sp.QQ.frac_field(t), groebner_backend=backend
        )
    assert cached is first
    assert profile.records == []
    second = compute_rur(
        [x * x - t, y - x],
        (x, y),
        domain=sp.QQ.frac_field(t, sp.Symbol("u")),
        groebner_backend=backend,
    )
    assert second.defining_polynomial.domain != first.defining_polynomial.domain
    if backend == "native":
        with pytest.raises(GroebnerLimitError):
            compute_rur(
                [x * x - t, y - x],
                (x, y),
                domain=sp.QQ.frac_field(t),
                groebner_backend=backend,
                groebner_limits=NativeF4Limits(max_terms=1),
            )
