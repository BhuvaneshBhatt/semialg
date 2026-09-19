from __future__ import annotations

from .blocks import QuantifierBlock, blocks_to_quantifiers, norm_quant_blocks, quantifiers_to_blocks
from .complete import (
    CellUnion,
    CompleteQEResult,
    QEDiagnostics,
    cells_to_formula,
    evaluate_formula_on_cell,
    qe_by_complete_cad,
    qe_from_cad,
)
from .prenex import eval_quantifier_free, qe_blocks, qe_parsed, qe_prenex, qe_prenex_suffix, qe_text

__all__ = [
    "CellUnion",
    "CompleteQEResult",
    "QEDiagnostics",
    "QuantifierBlock",
    "blocks_to_quantifiers",
    "cells_to_formula",
    "evaluate_formula_on_cell",
    "norm_quant_blocks",
    "qe_by_complete_cad",
    "qe_from_cad",
    "quantifiers_to_blocks",
    "eval_quantifier_free",
    "qe_blocks",
    "qe_parsed",
    "qe_prenex",
    "qe_prenex_suffix",
    "qe_text",
    "QuadraticVirtualSubstitutionResult",
    "VirtualSubstitutionError",
    "VirtualSubstitutionQEResult",
    "VirtualSubstitutionWitnessResult",
    "can_use_quadratic_vs",
    "eliminate_exists_quadratic_variable",
    "eliminate_quadratic_variable",
    "try_quadratic_virtual_substitution_qe",
    "try_quadratic_virtual_substitution_witness",
    "reconstruct_vs_value",
]

from .virtual_substitution import (
    QuadraticVirtualSubstitutionResult,
    VirtualSubstitutionError,
    VirtualSubstitutionQEResult,
    VirtualSubstitutionWitnessResult,
    can_use_quadratic_vs,
    eliminate_exists_quadratic_variable,
    eliminate_quadratic_variable,
    reconstruct_vs_value,
    try_quadratic_virtual_substitution_qe,
    try_quadratic_virtual_substitution_witness,
)
