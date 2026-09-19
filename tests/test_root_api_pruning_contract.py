"""Contracts for expert APIs intentionally pruned from the package root."""

from __future__ import annotations

from importlib import import_module

import pytest

import semialg
from semialg._api_policy import EXPERT_EXPORTS, PRIMARY_EXPORTS, PRUNED_ROOT_EXPERTS


def test_pruned_expert_names_are_absent_from_root_and_registered_as_expert():
    assert len(PRUNED_ROOT_EXPERTS) == 88
    assert PRUNED_ROOT_EXPERTS.isdisjoint(PRIMARY_EXPORTS)
    assert PRUNED_ROOT_EXPERTS <= set(EXPERT_EXPORTS)
    assert len(PRIMARY_EXPORTS) >= len(PRUNED_ROOT_EXPERTS)

    for name in PRUNED_ROOT_EXPERTS:
        assert name not in semialg.__all__
        with pytest.raises(AttributeError):
            getattr(semialg, name)


def test_pruned_expert_names_remain_public_from_their_defining_submodules():
    for name in PRUNED_ROOT_EXPERTS:
        module_name = EXPERT_EXPORTS[name]
        module = import_module(module_name, "semialg")
        assert getattr(module, name) is not None
