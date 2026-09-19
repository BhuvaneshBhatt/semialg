from __future__ import annotations

import random

import sympy as sp

import semialg
from semialg.algebraic.groebner_utils import compute_groebner_basis
from semialg.algebraic_decomposition import (
    associated_primes,
    certify_triangular_primality,
    primary_decomposition,
    verify_triangular_primality_certificate,
)
from semialg.algebraic_function_fields import (
    MonogenicFunctionField,
    RationalFunctionField,
    maybe_compress_primitive_element,
)
from semialg.cache_control import cache_report, clear_caches


def _basis_exprs(gb) -> tuple[sp.Expr, ...]:
    return tuple(sp.expand(poly.as_expr()) for poly in gb.polys)


def test_unified_cache_controls_clear_every_registered_family():
    x, y = sp.symbols("x y")
    clear_caches(include_sympy=True)

    # Populate GTZ/QQ canonical-basis and several exact algebraic caches.
    result = primary_decomposition(((x + y) ** 2, y * (x + y)), (x, y))
    assert result.complete
    assert semialg.replay_certificate(result).verified

    before = cache_report()
    assert before["gtz"]["qq"]["currsize"] >= 1
    for family in ("algebraic", "cad"):
        sizes = before[family]["sizes"]
        limits = before[family]["limits"]
        assert all(sizes[name] <= limits[name] for name in sizes)

    clear_caches(include_sympy=True)
    after = cache_report()
    assert all(info["currsize"] == 0 for info in after["gtz"].values())
    assert all(size == 0 for size in after["algebraic"]["sizes"].values())
    assert all(size == 0 for size in after["cad"]["sizes"].values())
    assert all(info["currsize"] == 0 for info in after["functools"].values())

    # Cache teardown cannot invalidate an existing exact certificate.
    assert semialg.replay_certificate(result).verified


def test_public_primitive_compression_heuristic_has_direct_contract_call():
    a = sp.Symbol("a")
    base = RationalFunctionField(tuple())
    field = MonogenicFunctionField(base, a, (sp.Integer(-2), sp.Integer(0), sp.Integer(1)))
    assert maybe_compress_primitive_element(field) is None


def test_public_triangular_primality_verifier_has_direct_contract_call():
    x, y = sp.symbols("x y")
    certificate = certify_triangular_primality((y - x**2,), (x, y))
    assert certificate.complete and certificate.prime
    assert verify_triangular_primality_certificate(certificate)


def test_modular_and_direct_groebner_differential_seeded_corpus():
    rng = random.Random(20260906)
    x, y = sp.symbols("x y")
    for _ in range(12):
        a, b, c, d = [rng.randint(-3, 3) for _ in range(4)]
        f = x**2 + a * x * y + b * y + c
        g = y**2 + d * x + a
        direct = compute_groebner_basis(
            (f, g), (x, y), order="grevlex", domain=sp.QQ, modular=False
        )
        modular = compute_groebner_basis(
            (f, g), (x, y), order="grevlex", domain=sp.QQ, modular=True
        )
        assert _basis_exprs(modular) == _basis_exprs(direct)


def test_generated_primary_decomposition_reconstruction_properties():
    rng = random.Random(6151)
    x, y = sp.symbols("x y")
    for _ in range(8):
        a = rng.choice((-2, -1, 1, 2))
        b = rng.choice((-2, -1, 0, 1, 2))
        c = a + rng.choice((1, 2, 3))
        d = b + rng.choice((1, 2))
        p = x + a * y + b
        q = x + c * y + d
        e1 = rng.choice((1, 2, 3))
        e2 = rng.choice((1, 2))
        source = sp.expand(p**e1 * q**e2)

        result = primary_decomposition((source,), (x, y))
        assert result.complete
        replay = semialg.replay_certificate(result)
        assert replay.verified
        primes = associated_primes((source,), (x, y))
        assert primes.complete
        assert len(primes.primes) == 2

        # Multiplication by a nonzero rational unit cannot change the ideal.
        scaled = primary_decomposition((sp.Rational(7, 3) * source,), (x, y))
        assert scaled.complete
        assert semialg.replay_certificate(scaled).verified
        scaled_primes = associated_primes((sp.Rational(7, 3) * source,), (x, y))
        assert scaled_primes.complete
        assert {tuple(map(sp.expand, prime)) for prime in scaled_primes.primes} == {
            tuple(map(sp.expand, prime)) for prime in primes.primes
        }


def test_repeated_workload_cache_bounds_and_teardown():
    x, y = sp.symbols("x y")
    clear_caches(include_sympy=True)
    for k in range(24):
        f = (x + (k % 5) * y + 1) * (x - y + (k % 3))
        result = primary_decomposition((sp.expand(f),), (x, y))
        assert result.complete
    report = cache_report()
    for family in ("algebraic", "cad"):
        assert all(
            report[family]["sizes"][name] <= report[family]["limits"][name]
            for name in report[family]["sizes"]
        )
    for info in report["gtz"].values():
        assert info["currsize"] <= info["maxsize"]
    for info in report["functools"].values():
        assert info["currsize"] <= info["maxsize"]
    clear_caches(include_sympy=True)


def test_functools_cache_registry_drives_lifecycle_and_reporting():
    from semialg.cache_control import FUNCTOOLS_CACHE_REGISTRY

    namespaces = tuple(spec.namespace for spec in FUNCTOOLS_CACHE_REGISTRY)
    sources = tuple(spec.source for spec in FUNCTOOLS_CACHE_REGISTRY)
    assert len(namespaces) == len(set(namespaces))
    assert len(sources) == len(set(sources))

    report = cache_report()["functools"]
    assert set(report) == set(namespaces)

    clear_caches(collect=False)
    cleared = cache_report()["functools"]
    assert all(info["currsize"] == 0 for info in cleared.values())


def test_every_lru_cache_is_registered_with_unified_teardown():
    """A new functools LRU must be added to the package-wide lifecycle hook."""
    import ast
    from pathlib import Path

    from semialg.cache_control import REGISTERED_LRU_CACHES

    root = Path(__file__).resolve().parents[1] / "src" / "semialg"
    found: set[tuple[str, str]] = set()
    for path in root.rglob("*.py"):
        tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
        for node in ast.walk(tree):
            if not isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)):
                continue
            for decorator in node.decorator_list:
                target = decorator.func if isinstance(decorator, ast.Call) else decorator
                if isinstance(target, ast.Name) and target.id == "lru_cache":
                    found.add((path.relative_to(root).as_posix(), node.name))
    assert found == set(REGISTERED_LRU_CACHES)


def test_every_bounded_lru_is_registered_with_unified_teardown():
    import ast
    from pathlib import Path

    from semialg.cache_control import BOUNDED_CACHE_REGISTRY

    root = Path(__file__).resolve().parents[1] / "src" / "semialg"
    found: set[str] = set()
    for path in root.rglob("*.py"):
        tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
        for node in ast.walk(tree):
            if not isinstance(node, ast.Call):
                continue
            func = node.func
            is_constructor = isinstance(func, ast.Name) and func.id == "BoundedLRU"
            is_policy_constructor = (
                isinstance(func, ast.Attribute)
                and isinstance(func.value, ast.Name)
                and func.value.id == "BoundedLRU"
                and func.attr in {"immutable", "heavyweight", "operation_local"}
            )
            if not (is_constructor or is_policy_constructor):
                continue
            if len(node.args) >= 2 and isinstance(node.args[1], ast.Constant):
                if isinstance(node.args[1].value, str):
                    found.add(node.args[1].value)
    assert found == set(BOUNDED_CACHE_REGISTRY)
