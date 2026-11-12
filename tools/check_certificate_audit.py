"""Aggregate mutation evidence without confusing constructor rejection with replay."""

from __future__ import annotations

import argparse
import json
from pathlib import Path


def summarize(records, required=None):
    """Collect field activation and rejection evidence across all certificate fixtures."""
    definitions = {}
    families = {}
    for record in records:
        root = record["certificate"]
        family = families.setdefault(
            root,
            {"fixtures": 0, "replay": 0, "well_formed_replay": 0, "schema": 0, "constructor": 0},
        )
        family["fixtures"] += 1
        rows = record["fields"]
        for row in rows:
            status = row["status"]
            if status in {"active", "inactive"}:
                field = definitions.setdefault(
                    row["definition"],
                    {"role": row["role"], "active": 0, "inactive": 0, "nonempty": 0, "rejected": 0},
                )
                field[status] += 1
                field["nonempty"] += int(row.get("nonempty", False))
                prefix = row["path"]
                field["rejected"] += int(
                    any(
                        r["status"].endswith("rejected")
                        and (
                            r["path"] == prefix
                            or r["path"].startswith(prefix + ".")
                            or r["path"].startswith(prefix + "[")
                        )
                        for r in rows
                    )
                )
            elif status == "replay-rejected":
                family["replay"] += 1
                family["well_formed_replay"] += int(row.get("well_formed", False))
            elif status == "schema-rejected":
                family["schema"] += 1
            elif status == "constructor-rejected":
                family["constructor"] += 1
    omissions = []
    for name, field in definitions.items():
        if not field["active"] or (field["role"] == "conditional" and not field["nonempty"]):
            omissions.append(f"inactive proof branch: {name}")
        elif not field["rejected"]:
            omissions.append(f"unmutated proof field: {name}")
    for name, family in families.items():
        if not family["well_formed_replay"]:
            omissions.append(f"no well-formed mathematical replay: {name}")
    if required is not None:
        omissions.extend(
            f"missing certificate family: {name}"
            for name in set(required["families"]) - families.keys()
        )
        omissions.extend(
            f"missing proof field: {name}"
            for name in set(required["proof_fields"]) - definitions.keys()
        )
    if not families or not definitions:
        omissions.append("no certificate mutation evidence")
    return {"families": families, "fields": definitions, "omissions": sorted(omissions)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("report", type=Path)
    parser.add_argument("--output", type=Path, default=Path("certificate-summary.json"))
    parser.add_argument(
        "--inventory",
        type=Path,
        default=Path(__file__).resolve().parents[1] / "tests/certificate_audit_inventory.json",
    )
    args = parser.parse_args()
    result = summarize(
        (json.loads(line) for line in args.report.read_text(encoding="utf-8").splitlines() if line),
        json.loads(args.inventory.read_text(encoding="utf-8")),
    )
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    for omission in result["omissions"]:
        print(omission)
    print(f"{len(result['families'])} families, {len(result['fields'])} proof fields")
    return int(bool(result["omissions"]))


if __name__ == "__main__":
    raise SystemExit(main())
