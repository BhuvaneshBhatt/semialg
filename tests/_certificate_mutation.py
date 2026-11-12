"""Explicit certificate-field policies and observable mutation coverage.

No field is skipped by its name. New dataclass fields require a policy entry.
Inactive optional branches and non-proof metadata appear in the audit report.
Unexpected verifier exceptions fail the test instead of counting as rejection.
"""

from __future__ import annotations

import json
import os
import re
from dataclasses import fields, is_dataclass, replace
from pathlib import Path

import sympy as sp

POLICIES = json.loads(
    (Path(__file__).with_name("certificate_fields.json")).read_text(encoding="utf-8")
)
AUDIT = []


def _policy(value):
    cls = type(value)
    name = cls.__module__ + "." + cls.__qualname__
    assert name in POLICIES, f"unclassified certificate dataclass: {name}"
    policy = POLICIES[name]
    actual = {field.name for field in fields(value)}
    assert set(policy) == actual, (
        f"certificate policy mismatch for {name}: missing={actual - set(policy)}, stale={set(policy) - actual}"
    )
    for field, entry in policy.items():
        assert entry["role"] in {
            "proof",
            "conditional",
            "provenance",
            "derived",
            "diagnostic",
            "context",
        }, (name, field)
        if entry["role"] != "proof":
            assert entry.get("reason"), (name, field)
    return policy


def _atom(value):
    if isinstance(value, bool):
        return not value, False
    if isinstance(value, (int, sp.Integer)):
        return value + 1, False
    if isinstance(value, sp.Symbol):
        return sp.Dummy("altered_" + value.name), True
    if isinstance(value, sp.Poly):
        return sp.Poly(value.as_expr() + 1, *value.gens, domain=value.domain), False
    if isinstance(value, sp.MatrixBase):
        matrix = sp.MutableDenseMatrix(value)
        if matrix.rows and matrix.cols:
            matrix[0, 0] += 1
        else:
            matrix = sp.Matrix([[1]])
        return sp.ImmutableMatrix(matrix), False
    if isinstance(value, sp.Basic):
        return sp.expand(value + 1), False
    if isinstance(value, str):
        return value + "_altered", False
    if value is None:
        return sp.S.One, True
    raise TypeError(f"no registered mutator for {type(value).__name__}")


def _mutations(value, path, report, root_name):
    """Traverse proof fields while preserving enclosing dataclass and container shapes."""
    if is_dataclass(value) and not isinstance(value, type):
        for field, declared in _policy(value).items():
            entry = declared
            field_path = path + "." + field if path else field
            normalized = re.sub(r"\[\d+\]", "[]", field_path)
            entry = POLICIES.get("$overrides", {}).get(root_name, {}).get(normalized, entry)
            current = getattr(value, field)
            definition = type(value).__module__ + "." + type(value).__qualname__ + "." + field
            child = path + "." + field if path else field
            if entry["role"] not in {"proof", "conditional"}:
                report.append({"path": child, "status": entry["role"], "reason": entry["reason"]})
                continue
            if entry["role"] == "conditional" and current is None:
                report.append(
                    {
                        "path": child,
                        "definition": definition,
                        "status": "inactive",
                        "role": entry["role"],
                        "reason": entry["reason"],
                    }
                )
                continue
            report.append(
                {
                    "path": child,
                    "definition": definition,
                    "status": "active",
                    "role": entry["role"],
                    "nonempty": current is not None and current != () and current != frozenset(),
                }
            )
            for subpath, changed, malformed in _mutations(current, child, report, root_name):
                try:
                    mutant = replace(value, **{field: changed})
                except (TypeError, ValueError) as exc:
                    report.append(
                        {
                            "path": subpath,
                            "status": "constructor-rejected",
                            "exception": type(exc).__name__,
                        }
                    )
                    continue
                yield subpath, mutant, malformed
        return
    if isinstance(value, (tuple, frozenset)):
        sequence = tuple(value)
        if not sequence:
            yield (
                path + "[empty]",
                (sp.S.One,) if isinstance(value, tuple) else frozenset({sp.S.One}),
                True,
            )
        for i, item in enumerate(sequence):
            for child, changed, malformed in _mutations(item, f"{path}[{i}]", report, root_name):
                sequence_changed = (*sequence[:i], changed, *sequence[i + 1 :])
                yield (
                    child,
                    sequence_changed if isinstance(value, tuple) else frozenset(sequence_changed),
                    malformed,
                )
        return
    try:
        changed, malformed = _atom(value)
    except TypeError as exc:
        raise AssertionError(f"uncovered proof field {path}: {exc}") from exc
    yield path, changed, malformed


def field_mutations(certificate, *, ignored=frozenset()):
    """Yield field mutations; exclusions must already be explicit in the class policy."""
    assert is_dataclass(certificate)
    policy = _policy(certificate)
    assert all(policy[name]["role"] not in {"proof", "conditional"} for name in ignored), ignored
    report = []
    yield from _mutations(
        certificate, "", report, type(certificate).__module__ + "." + type(certificate).__qualname__
    )


def assert_certificate_fields_reject_mutation(certificate, verifier, *, ignored=frozenset()):
    """Check registered proof mutations and record every classified or inactive field."""
    assert verifier(certificate), f"baseline {type(certificate).__name__} does not verify"
    policy = _policy(certificate)
    assert all(policy[name]["role"] not in {"proof", "conditional"} for name in ignored), ignored
    report = []
    for path, mutant, malformed in _mutations(
        certificate, "", report, type(certificate).__module__ + "." + type(certificate).__qualname__
    ):
        try:
            accepted = bool(verifier(mutant))
        except (TypeError, ValueError, sp.polys.polyerrors.BasePolynomialError) as exc:
            if not malformed:
                raise
            report.append(
                {"path": path, "status": "schema-rejected", "exception": type(exc).__name__}
            )
            continue
        assert not accepted, f"altered proof accepted: {type(certificate).__name__}.{path}"
        report.append({"path": path, "status": "replay-rejected", "well_formed": not malformed})
    assert any(row["status"].endswith("rejected") for row in report)
    record = {
        "certificate": type(certificate).__module__ + "." + type(certificate).__qualname__,
        "fields": report,
    }
    AUDIT.append(record)
    destination = os.environ.get("CERTIFICATE_AUDIT_REPORT")
    if destination:
        with Path(destination).open("a", encoding="utf-8") as stream:
            stream.write(json.dumps(record) + "\n")
    return record


def assert_field_replacements(certificate, verifier, replacements):
    """Audit explicit whole-field mathematical mutations of structured payloads.

    Native CAD payloads contain graphs and read-only mappings. Their specialized
    geometric replay tests audit internals; this helper binds every enclosing
    statement field without reducing graphs to invalid scalar placeholders.
    """
    policy = _policy(certificate)
    required = {name for name, entry in policy.items() if entry["role"] in {"proof", "conditional"}}
    assert set(replacements) == required
    assert verifier(certificate)
    family = type(certificate).__module__ + "." + type(certificate).__qualname__
    report = []
    for name, entry in policy.items():
        if name not in replacements:
            report.append({"path": name, "status": entry["role"], "reason": entry["reason"]})
            continue
        report.append(
            {
                "path": name,
                "definition": family + "." + name,
                "status": "active",
                "role": entry["role"],
                "nonempty": True,
            }
        )
        mutant = replace(certificate, **{name: replacements[name]})
        assert not verifier(mutant), f"altered statement accepted: {family}.{name}"
        report.append({"path": name, "status": "replay-rejected", "well_formed": True})
    record = {"certificate": family, "fields": report}
    AUDIT.append(record)
    destination = os.environ.get("CERTIFICATE_AUDIT_REPORT")
    if destination:
        with Path(destination).open("a", encoding="utf-8") as stream:
            stream.write(json.dumps(record) + "\n")
    return record
