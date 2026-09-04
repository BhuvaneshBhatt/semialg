from __future__ import annotations

import inspect

import semialg

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
        parameters = inspect.signature(getattr(semialg, name)).parameters
        actual = {parameter: parameters[parameter].default for parameter in expected}
        assert actual == expected, name
