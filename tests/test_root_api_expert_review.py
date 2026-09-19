"""Contracts for specialist workflows pruned from the package root."""

from importlib import import_module

import pytest

import semialg
from semialg._api_policy import EXPERT_EXPORTS, PRIMARY_EXPORTS, ROOT_EXPERT_CANDIDATES

SPECIALIST_WORKFLOWS = {
    "thom_encoding": ".algebraic",
    "thom_encodings": ".algebraic",
    "roadmap": ".roadmaps",
    "dimension_strata": ".topology.semialgebraic",
    "component_decomposition": ".topology.semialgebraic",
    "hardt_trivialization": ".topology.semialgebraic",
    "triangulate_region": ".topology.semialgebraic",
    "simplicial_betti_numbers": ".topology.semialgebraic",
    "triangulation_betti_numbers": ".topology.semialgebraic",
    "root_count_conditions": ".parameters",
    "bounded_parametric_cover": ".parametric_geometry",
    "intrinsic_parametric_cover": ".parametric_geometry",
    "parametric_map_degree": ".map_degree",
}


def test_specialist_review_is_complete_and_namespace_only():
    assert ROOT_EXPERT_CANDIDATES == {}
    for name, namespace in SPECIALIST_WORKFLOWS.items():
        assert name not in PRIMARY_EXPORTS
        assert EXPERT_EXPORTS[name] == namespace
        assert name not in semialg.__all__
        with pytest.raises(AttributeError):
            getattr(semialg, name)
        assert getattr(import_module(namespace, "semialg"), name) is not None
