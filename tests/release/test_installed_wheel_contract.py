"""Release contract executed against an installed wheel outside the source tree."""

from __future__ import annotations

import importlib.util
import os
from importlib import metadata
from pathlib import Path

import pytest
import sympy as sp

from semialg.cache_control import clear_caches

if os.environ.get("SEMIALG_INSTALLED_WHEEL_TEST") != "1":
    pytest.skip(
        "installed-wheel contract runs through scripts/run_installed_wheel_contract.py",
        allow_module_level=True,
    )

import semialg
from semialg.algebraic import (
    certify_polynomial_root_interval,
    verify_polynomial_root_interval_certificate,
)


def test_import_resolves_to_installed_package_not_source_checkout():
    module_path = Path(semialg.__file__).resolve()
    source_root = os.environ.get("SEMIALG_SOURCE_ROOT")
    if source_root:
        assert not module_path.is_relative_to((Path(source_root) / "src" / "semialg").resolve())
    assert "site-packages" in str(module_path) or os.environ.get("SEMIALG_TARGET_INSTALL") in str(
        module_path
    )


def test_distribution_metadata_and_public_surface_are_present():
    assert metadata.version("semialg") == semialg.__version__
    required = {
        "cad",
        "function_range",
        "is_satisfiable",
        "region_union",
        "replay_certificate",
    }
    assert required <= set(semialg.__all__)
    assert all(hasattr(semialg, name) for name in required)


def test_release_does_not_ship_repository_test_or_development_packages():
    package_root = Path(semialg.__file__).resolve().parent
    assert not (package_root / "tests").exists()
    assert importlib.util.find_spec("tests") is None


def test_optional_transcendental_subpackage_is_shipped():
    assert importlib.util.find_spec("semialg.solve.transcendental") is not None


def test_certificate_replay_survives_cache_reset_in_installed_wheel():
    x = sp.Symbol("x")
    certificate = certify_polynomial_root_interval(x**2 - 2, 0, 2, var=x)
    clear_caches(include_sympy=False, collect=False)
    assert verify_polynomial_root_interval_certificate(
        certificate,
        polynomial=x**2 - 2,
        var=x,
        left=0,
        right=2,
        include_left=True,
        include_right=True,
    )


def test_representative_exact_public_computations():
    x = sp.Symbol("x", real=True)
    assert semialg.is_satisfiable((x >= -1) & (x <= 1), (x,))
    assert semialg.is_equal(x**2 <= 1, (x >= -1) & (x <= 1), (x,))
    result = semialg.polynomial_nonnegative(x**2 + 1, (x,), sos_backend="none")
    assert result is True


def test_pruned_expert_api_is_not_reexported_from_installed_package_root():
    from importlib import import_module

    from semialg._api_policy import EXPERT_EXPORTS, PRUNED_ROOT_EXPERTS

    assert len(PRUNED_ROOT_EXPERTS) == 88
    for name in PRUNED_ROOT_EXPERTS:
        assert name not in semialg.__all__
        assert not hasattr(semialg, name)
        module = import_module(EXPERT_EXPORTS[name], "semialg")
        assert getattr(module, name) is not None
