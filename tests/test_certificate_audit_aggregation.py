"""Audit aggregation must expose inactive branches and schema-only evidence."""

import importlib.util
from pathlib import Path

_spec = importlib.util.spec_from_file_location(
    "audit_checker", Path(__file__).parents[1] / "tools/check_certificate_audit.py"
)
_checker = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_checker)


def _record(rows):
    return {"certificate": "Proof", "fields": rows}


def test_inactive_branch():
    result = _checker.summarize(
        [
            _record(
                [
                    {
                        "definition": "Proof.branch",
                        "path": "branch",
                        "status": "inactive",
                        "role": "conditional",
                    }
                ]
            )
        ]
    )
    assert result["omissions"] == [
        "inactive proof branch: Proof.branch",
        "no well-formed mathematical replay: Proof",
    ]


def test_schema_is_not_replay():
    rows = [
        {
            "definition": "Proof.value",
            "path": "value",
            "status": "active",
            "role": "proof",
            "nonempty": True,
        },
        {"path": "value", "status": "schema-rejected"},
    ]
    result = _checker.summarize([_record(rows)])
    assert result["families"]["Proof"]["schema"] == 1
    assert result["families"]["Proof"]["well_formed_replay"] == 0
    assert result["omissions"] == ["no well-formed mathematical replay: Proof"]


def test_cross_fixture_activation():
    inactive = {
        "definition": "Proof.branch",
        "path": "branch",
        "status": "inactive",
        "role": "conditional",
    }
    active = dict(inactive, status="active", nonempty=True)
    rejected = {"path": "branch.value", "status": "replay-rejected", "well_formed": True}
    result = _checker.summarize([_record([inactive]), _record([active, rejected])])
    assert result["omissions"] == []
    assert result["fields"]["Proof.branch"]["active"] == 1
    assert result["fields"]["Proof.branch"]["inactive"] == 1


def test_missing_family_is_visible():
    required = {"families": ["DeletedProof"], "proof_fields": ["DeletedProof.value"]}
    result = _checker.summarize([], required)
    assert "missing certificate family: DeletedProof" in result["omissions"]
    assert "missing proof field: DeletedProof.value" in result["omissions"]
