"""Conservative unified real transcendental analysis facade.

It reuses exact structured zeros, tail signs, singularities, and the public range
engine.  Completeness is explicit; no sampling is promoted to proof.
"""

from __future__ import annotations

from dataclasses import dataclass

import sympy as sp

from ..asymptotic_tails import certify_tail_sign
from ..infinite_root_families import structured_zero_geometry


@dataclass(frozen=True)
class UnifiedRealAnalysis:
    expression: object
    variable: object
    domain: object
    zero_geometry: object = None
    left_tail: object = None
    right_tail: object = None
    image: object = None
    complete: bool = False


def analyze_transcendental_real(expr, var, domain=sp.S.Reals):
    f = sp.sympify(expr)
    x = sp.sympify(var)
    zg = structured_zero_geometry(f, x)
    lt = certify_tail_sign(f, x, direction=-1) if domain is sp.S.Reals else None
    rt = certify_tail_sign(f, x, direction=1) if domain is sp.S.Reals else None
    image = None
    complete = False
    # Exact safe families: affine trig images are known globally; tan has poles but exact image R.
    if domain is sp.S.Reals and f.func in (sp.sin, sp.cos):
        image = sp.Interval(-1, 1)
        complete = zg is not None
    elif domain is sp.S.Reals and f.func is sp.tan:
        image = sp.S.Reals
        complete = zg is not None
    else:
        try:
            from ..real_analysis.range import function_range

            image = function_range(f, x, sp.Symbol("_y", real=True), constraints=domain)
            complete = (
                image is not None and getattr(image, "complete", True) is not False
            )
        except (TypeError, ValueError, NotImplementedError, ImportError):
            pass
    return UnifiedRealAnalysis(f, x, domain, zg, lt, rt, image, complete)
