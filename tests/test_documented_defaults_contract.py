from __future__ import annotations

import inspect
import re
from pathlib import Path

import semialg
from semialg.parameters import root_count_conditions

ROOT = Path(__file__).resolve().parents[1]

EXPECTED_DEFAULTS = {
    "cad": {
        "output": "formula",
        "strategy": "auto",
        "domain": "reals",
        "diagnostics": True,
        "strict": False,
        "return_result": False,
        "formula_form": "nested",
        "max_formula_terms": 512,
        "max_preprocess_aux_vars": 3,
    },
    "function_range": {
        "constraints": None,
        "variables": None,
        "domain": "reals",
        "method": "qe",
        "return_result": False,
        "parameters": None,
        "return_stratified": False,
        "eliminate_quantifiers": False,
    },
    "root_count_conditions": {
        "parameters": None,
        "return_result": False,
        "return_stratified": False,
    },
    "classify_real_roots": {"parameters": None},
    "semialgebraic_minimize": {
        "constraints": None,
        "variables": None,
        "domain": "reals",
        "return_result": False,
        "certification": "auto",
        "range_cost_limit": 2500,
        "recursion_limit": 4,
        "max_boolean_branches": 32,
        "parameters": None,
        "return_stratified": False,
        "eliminate_quantifiers": False,
    },
    "semialgebraic_maximize": {
        "constraints": None,
        "variables": None,
        "domain": "reals",
        "return_result": False,
        "certification": "auto",
        "range_cost_limit": 2500,
        "recursion_limit": 4,
        "max_boolean_branches": 32,
        "parameters": None,
        "return_stratified": False,
        "eliminate_quantifiers": False,
    },
    "integrate_over_region": {
        "bounds": None,
        "method": "symbolic",
        "precision": 50,
        "measure_dimension": "ambient",
        "return_result": False,
        "parameters": None,
        "return_stratified": False,
    },
}


def test_documented_high_level_defaults_match_python_signatures() -> None:
    for name, expected in EXPECTED_DEFAULTS.items():
        function = (
            root_count_conditions if name == "root_count_conditions" else getattr(semialg, name)
        )
        parameters = inspect.signature(function).parameters
        actual = {parameter: parameters[parameter].default for parameter in expected}
        assert actual == expected, name


def test_default_contract_page_contains_every_guarded_default() -> None:
    text = (ROOT / "docs" / "reference" / "defaults_contract.md").read_text(encoding="utf-8")
    normalized = text.replace('"', "'")
    for name, defaults in EXPECTED_DEFAULTS.items():
        documented_name = (
            "semialg.parameters.root_count_conditions" if name == "root_count_conditions" else name
        )
        assert f"`{documented_name}`" in text
        for parameter, value in defaults.items():
            token = f"`{parameter}={value!r}`".replace('"', "'")
            assert re.search(re.escape(token), normalized), f"{name}: missing {token}"
