from __future__ import annotations

import inspect

import semialg
from semialg._api_policy import (
    DETAILED_WRAPPER_APIS,
    EVERYDAY_EXPORTS,
    PRIMARY_EXPORTS,
    STRUCTURED_RESULT_APIS,
)


def test_everyday_exports_are_a_small_supported_subset_of_root_api():
    assert EVERYDAY_EXPORTS <= set(PRIMARY_EXPORTS)
    assert len(EVERYDAY_EXPORTS) < len(PRIMARY_EXPORTS) // 2


def test_sole_result_returning_root_functions_are_explicitly_classified():
    candidates = set()
    for name in PRIMARY_EXPORTS:
        obj = getattr(semialg, name)
        if not inspect.isfunction(obj):
            continue
        signature = inspect.signature(obj)
        if "return_result" in signature.parameters:
            continue
        if "Result" in str(signature.return_annotation):
            candidates.add(name)

    classified = STRUCTURED_RESULT_APIS | DETAILED_WRAPPER_APIS
    assert candidates <= classified


def test_result_policy_groups_are_root_names_and_disjoint():
    assert STRUCTURED_RESULT_APIS <= set(PRIMARY_EXPORTS)
    assert DETAILED_WRAPPER_APIS <= set(PRIMARY_EXPORTS)
    assert STRUCTURED_RESULT_APIS.isdisjoint(DETAILED_WRAPPER_APIS)
