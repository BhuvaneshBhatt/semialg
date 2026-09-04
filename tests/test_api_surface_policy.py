"""Contracts for the tiered package API surface."""

from __future__ import annotations

import semialg
from semialg._api_policy import (
    ALL_API_NAMES,
    EXPERT_EXPORTS,
    INTERNAL_EXPORTS,
    PRIMARY_EXPORTS,
)


def test_api_policy_partitions_the_audited_surface():
    primary = set(PRIMARY_EXPORTS)
    expert = set(EXPERT_EXPORTS)
    internal = set(INTERNAL_EXPORTS)

    assert primary.isdisjoint(expert)
    assert primary.isdisjoint(internal)
    assert expert.isdisjoint(internal)
    assert ALL_API_NAMES == primary | expert | internal


def test_only_primary_names_are_exported_at_package_root():
    assert set(semialg.__all__) == {"__version__", *PRIMARY_EXPORTS}
    assert not (set(EXPERT_EXPORTS) & set(semialg.__all__))
    assert not (set(INTERNAL_EXPORTS) & set(semialg.__all__))


def test_dir_exposes_only_primary_names_and_module_metadata():
    visible = set(dir(semialg))
    assert set(PRIMARY_EXPORTS) <= visible
    assert not (set(EXPERT_EXPORTS) & visible)
    assert not (set(INTERNAL_EXPORTS) & visible)


def test_cad_function_survives_algorithm_package_import():
    """Importing CAD internals must not replace the root ``cad`` function."""
    import semialg.cad_algorithms as cad_algorithms
    from semialg import cad as root_cad

    assert callable(root_cad)
    assert semialg.cad is root_cad
    assert cad_algorithms.__name__ == "semialg.cad_algorithms"
    assert semialg.cad is root_cad
