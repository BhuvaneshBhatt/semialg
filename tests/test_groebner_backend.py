import pytest
import sympy as sp
from algroots.native_f4 import compiled_available

from semialg.algebraic.groebner_utils import compute_groebner_basis


@pytest.mark.skipif(not compiled_available(), reason="native F4 extension unavailable")
def test_explicit_backend():
    x, y = sp.symbols("x y")
    system = (x**2 - y, y**2 - 1)
    result = compute_groebner_basis(iter(system), (x, y), backend="native")
    assert result == sp.groebner(system, x, y, domain=sp.QQ, order="grevlex")


def test_ambiguous_options():
    x = sp.Symbol("x")
    with pytest.raises(ValueError, match="domain"):
        compute_groebner_basis((x,), (x,), backend="native", domain=sp.QQ, modulus=101)


@pytest.mark.skipif(not compiled_available(), reason="native F4 extension unavailable")
def test_rur_cache_selection(monkeypatch):
    from semialg.algebraic.cache import clear_algebraic_caches
    from semialg.algebraic.rational_univariate import construction

    x, y = sp.symbols("x y")
    calls = []
    original = construction.compute_rational_univariate_representation

    def record(*args, **kwargs):
        calls.append(kwargs["groebner_backend"])
        return original(*args, **kwargs)

    monkeypatch.setattr(construction, "compute_rational_univariate_representation", record)
    clear_algebraic_caches()
    for backend in ["native", "sympy", "native"]:
        result = construction.compute_rur([x * x - y, y * y - 1], (x, y), groebner_backend=backend)
        assert result.solution_count == 4
    assert calls == ["native", "sympy"]
    with pytest.raises(ValueError, match="backend"):
        construction.compute_rur([x * x - y, y * y - 1], (x, y), groebner_backend="unknown")


@pytest.mark.parametrize("domain", [sp.QQ, sp.GF(2), sp.GF(101)])
@pytest.mark.parametrize("order", ["lex", "grlex", "grevlex"])
@pytest.mark.parametrize("backend", ["sympy", "native"])
def test_basis_parity(domain, order, backend):
    from algroots import groebner_backends, groebner_basis

    if not groebner_backends()[backend].available:
        pytest.skip(f"{backend} is not installed")
    x, y = sp.symbols("x y")
    equations = [x * x - y, y * y - 1]
    basis = compute_groebner_basis(equations, (x, y), domain=domain, order=order, backend=backend)
    assert basis == groebner_basis(equations, (x, y), domain=domain, order=order, backend="sympy")
    assert basis.reduce(x**4 - 1)[1] == 0


@pytest.mark.parametrize("backend", ["sympy", "native"])
@pytest.mark.parametrize("multiple", [False, True])
def test_rur_parity(backend, multiple):
    from algroots import GroebnerLimitError, NativeF4Limits, groebner_backends

    from semialg.algebraic.rational_univariate.construction import compute_rur

    if not groebner_backends()[backend].available:
        pytest.skip(f"{backend} is not installed")
    x, y = sp.symbols("x y")
    equations = [x * x if multiple else x * x - 1, y - x]
    result = compute_rur(equations, (x, y), groebner_backend=backend)
    reference = compute_rur(equations, (x, y), groebner_backend="sympy")
    assert result.solution_count == reference.solution_count == (1 if multiple else 2)
    assert result.defining_polynomial == reference.defining_polynomial
    if backend == "native":
        with pytest.raises(GroebnerLimitError):
            compute_rur(
                equations,
                (x, y),
                groebner_backend=backend,
                groebner_limits=NativeF4Limits(max_terms=1),
            )


@pytest.mark.parametrize("backend", ["native", "sympy"])
def test_primary_policy(backend):
    from algroots import GroebnerBackendUnavailable
    from algroots.profiling import GroebnerProfiler

    from semialg.algebraic_decomposition import primary_decomposition

    x, y = sp.symbols("x y")
    options = {} if backend == "native" else {"groebner_backend": "sympy"}
    with GroebnerProfiler() as profile:
        if backend == "native" and not compiled_available():
            with pytest.raises(GroebnerBackendUnavailable):
                primary_decomposition([(x + y) ** 2, y * (x + y)], (x, y), **options)
        else:
            result = primary_decomposition([(x + y) ** 2, y * (x + y)], (x, y), **options)
            assert result.complete
            assert sorted(c.dimension for c in result.components) == [0, 1]
    if backend == "native":
        assert profile.records
    assert {row["backend"] for row in profile.records} <= {backend}


def test_policy_reset():
    from algroots import NativeF4Limits

    from semialg.algebraic.groebner_utils import (
        _backend_policy,
        _selected_backend,
        _selected_limits,
    )

    @_backend_policy
    def probe(*, groebner_backend="native", groebner_limits=None):
        assert _selected_backend(None) == groebner_backend
        assert _selected_limits(None) == groebner_limits
        raise RuntimeError("interrupt")

    with pytest.raises(RuntimeError, match="interrupt"):
        probe(groebner_backend="sympy")
    assert _selected_backend(None) == "native"
    assert _selected_limits(None) is None
    with pytest.raises(ValueError, match="limits"):
        probe(groebner_backend="sympy", groebner_limits=NativeF4Limits())
