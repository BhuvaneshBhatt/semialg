"""Experimental transcendental preprocessing and solving helpers.

These routines implement sound fragments and heuristics outside the certified
real-polynomial CAD/QE core.
"""

from .cleanup import (
    CleanupResult,
    finite_points_form,
    recon_solved_points,
    recon_univar_intv_form,
    remove_redundant_disjunc,
)
from .engine import (
    TranscendentalSolveTrace,
    TransReductionResult,
    reduce_trans_problem,
)
from .families import (
    TransFamDetection,
    TransFamHandler,
    classify_trans_fams,
    default_trans_handlers,
)
from .periodic import (
    PeriodicBoundingResult,
    compute_periodic_window,
    detect_real_period,
    find_periodic_variables,
    periodic_intv_form,
    recon_periodic_domain,
    recon_periodic_represent,
)
from .preprocess import (
    FamilyReplacementStep,
    QuantifierDispatchPlan,
    TransPrepResult,
    build_quantifier_plan,
    prep_trans_problem,
    replace_function_auxilia,
    simp_piecewise_subexprs,
)
from .quantifier_elimination import (
    QuantElimResult,
    eliminate_lead_block,
)
from .roots import (
    CertifiedIntervalRoot,
    RootIsolationResult,
    SampledTruthDecomp,
    decomp_univar_inequality,
    evaluate_form_points,
    isolate_univar_roots,
)
from .state import (
    QuantifierBlock,
    TransProblemState,
    build_trans_state,
    norm_quant_blocks,
)
from .system_roots import (
    CertifiedPoint,
    CompletenessCertificate,
    SearchBox,
    SystemRootFallbackResult,
    orchestrate_trans_search,
    solve_bounded_trans_sys,
)

EXPERIMENTAL = True

__all__ = [
    "QuantifierBlock",
    "TransProblemState",
    "norm_quant_blocks",
    "build_trans_state",
    "FamilyReplacementStep",
    "QuantifierDispatchPlan",
    "TransPrepResult",
    "simp_piecewise_subexprs",
    "build_quantifier_plan",
    "replace_function_auxilia",
    "prep_trans_problem",
    "CertifiedIntervalRoot",
    "RootIsolationResult",
    "SampledTruthDecomp",
    "isolate_univar_roots",
    "evaluate_form_points",
    "decomp_univar_inequality",
    "PeriodicBoundingResult",
    "detect_real_period",
    "compute_periodic_window",
    "periodic_intv_form",
    "recon_periodic_represent",
    "recon_periodic_domain",
    "find_periodic_variables",
    "TransFamDetection",
    "TransFamHandler",
    "classify_trans_fams",
    "default_trans_handlers",
    "QuantElimResult",
    "eliminate_lead_block",
    "CertifiedPoint",
    "CompletenessCertificate",
    "SystemRootFallbackResult",
    "SearchBox",
    "orchestrate_trans_search",
    "solve_bounded_trans_sys",
    "CleanupResult",
    "finite_points_form",
    "remove_redundant_disjunc",
    "recon_solved_points",
    "recon_univar_intv_form",
    "TransReductionResult",
    "TranscendentalSolveTrace",
    "reduce_trans_problem",
]
