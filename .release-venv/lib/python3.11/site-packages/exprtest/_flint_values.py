"""Direct classification of python-flint scalar values."""

from __future__ import annotations

try:
    import flint
except ImportError:  # exact-only source use without the optional runtime backend
    flint = None

from ._result import Verdict, ZeroClassification


def classify_flint_value(value):
    """Classify supported FLINT scalars without converting through SymPy.

    Returns ``None`` when ``value`` is not a supported FLINT scalar.
    """
    if flint is None:
        return None
    if isinstance(value, flint.fmpz):
        zero = value == 0
        return ZeroClassification(
            Verdict.ZERO_PROVEN if zero else Verdict.NONZERO_PROVEN,
            "flint-exact",
            detail="exact FLINT integer compared directly with zero",
            evidence="exact-arithmetic",
        )
    if isinstance(value, flint.fmpq):
        zero = value == 0
        return ZeroClassification(
            Verdict.ZERO_PROVEN if zero else Verdict.NONZERO_PROVEN,
            "flint-exact",
            detail="exact FLINT rational compared directly with zero",
            evidence="exact-arithmetic",
        )

    qqbar_type = getattr(flint, "qqbar", None)
    if qqbar_type is not None and isinstance(value, qqbar_type):
        zero = value == 0
        return ZeroClassification(
            Verdict.ZERO_PROVEN if zero else Verdict.NONZERO_PROVEN,
            "flint-exact",
            detail="exact FLINT algebraic number compared directly with zero",
            evidence="exact-algebraic",
        )

    arb_type = getattr(flint, "arb", None)
    acb_type = getattr(flint, "acb", None)
    if (arb_type is not None and isinstance(value, arb_type)) or (
        acb_type is not None and isinstance(value, acb_type)
    ):
        if not value.is_finite:
            return ZeroClassification(
                Verdict.UNKNOWN,
                "flint-ball",
                detail="input FLINT ball is non-finite",
                evidence="rigorous-enclosure",
            )
        if 0 not in value:
            return ZeroClassification(
                Verdict.NONZERO_PROVEN,
                "flint-ball",
                detail=f"input FLINT ball {value} rigorously excludes zero",
                evidence="rigorous-enclosure",
            )
        if value.rad() == 0:
            return ZeroClassification(
                Verdict.ZERO_PROVEN,
                "flint-ball",
                detail="input FLINT ball has zero radius and contains zero",
                evidence="exact-enclosure",
            )
        return ZeroClassification(
            Verdict.UNKNOWN,
            "flint-ball",
            detail="input FLINT ball contains both zero and nonzero values",
            evidence="rigorous-enclosure",
            enclosure_history=(str(value),),
        )
    return None
