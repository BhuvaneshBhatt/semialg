"""Canonical CAD projection mode names used by native backends."""

PROJECTION_COLLINS = "collins"
PROJECTION_MCCALLUM = "mccallum"
PROJECTION_LAZARD = "lazard"
PROJECTION_TTICAD = "tticad"
VALID_PROJECTION_MODES = frozenset(
    {
        PROJECTION_COLLINS,
        PROJECTION_MCCALLUM,
        PROJECTION_LAZARD,
        PROJECTION_TTICAD,
    }
)

__all__ = [
    "PROJECTION_COLLINS",
    "PROJECTION_MCCALLUM",
    "PROJECTION_LAZARD",
    "PROJECTION_TTICAD",
    "VALID_PROJECTION_MODES",
]
