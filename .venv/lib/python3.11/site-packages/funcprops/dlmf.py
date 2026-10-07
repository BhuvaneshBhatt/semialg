"""Access to the generated DLMF property-table audit manifest."""

from __future__ import annotations

import json
from functools import lru_cache
from importlib.resources import files

from .real_analysis.special_functions import DLMFFact, dlmf_facts


@lru_cache(maxsize=1)
def _manifest():
    path = files("funcprops").joinpath("data/dlmf_audit.json")
    return json.loads(path.read_text(encoding="utf-8"))


def dlmf_audit_report(family: str | None = None):
    """Return the generated DLMF audit summary or one family's audit row."""
    data = _manifest()
    if family is None:
        return data
    for row in data["families"]:
        if row["family_id"] == family:
            return row
    return None


__all__ = ["DLMFFact", "dlmf_audit_report", "dlmf_facts"]
