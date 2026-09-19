from __future__ import annotations

from .eliminate import (
    can_use_quadratic_vs,
    eliminate_exists_quadratic_variable,
    eliminate_quadratic_variable,
    try_quadratic_virtual_substitution_qe,
)
from .results import (
    QuadraticVirtualSubstitutionResult,
    VirtualSubstitutionError,
    VirtualSubstitutionQEResult,
    VirtualSubstitutionWitnessResult,
)
from .substitution import (
    substitute_infinity,
    substitute_perturbed_quadratic_root,
    substitute_quadratic_root,
)
from .witness import (
    reconstruct_vs_value,
    try_quadratic_virtual_substitution_witness,
)

__all__ = [
    "QuadraticVirtualSubstitutionResult",
    "VirtualSubstitutionError",
    "VirtualSubstitutionQEResult",
    "VirtualSubstitutionWitnessResult",
    "can_use_quadratic_vs",
    "eliminate_exists_quadratic_variable",
    "eliminate_quadratic_variable",
    "reconstruct_vs_value",
    "substitute_infinity",
    "substitute_perturbed_quadratic_root",
    "substitute_quadratic_root",
    "try_quadratic_virtual_substitution_qe",
    "try_quadratic_virtual_substitution_witness",
]
