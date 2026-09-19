"""Tunable constants shared across the cascade stages."""

# --- Numeric (Arb ball arithmetic) ----------------------------------------
INITIAL_PREC_BITS = 80  # starting Arb precision
MAX_PREC_BITS = 20_000  # give up escalating beyond this
PREC_ESCALATION_FACTOR = 4  # multiply precision by this on retry
GUARD_BITS = 30  # extra bits of separation required to trust a verdict

# How many *consecutive* escalations must show the ball radius shrinking by
# at least half the precision increase (in bits) before we trust "still
# contains 0" as evidence of an exact zero, rather than just numerical noise
# that happens not to have resolved yet.
REQUIRED_SHRINKS = 2
MIN_SHRINK_FRACTION = 0.5

# --- Nonzero witness probing -----------------------------------------------
NUM_WITNESS_POINTS = 3  # independent random points tried before giving up
WITNESS_PREC_BITS = 100

# --- Schwartz-Zippel / finite-field identity testing -----------------------
SZ_TARGET_FALSE_POSITIVE = 1e-9  # target Schwartz-Zippel error probability
SZ_FIELD_MIN_BITS = 61  # size of the prime field used for SZ testing
SZ_MAX_TRIALS = 8  # maximum independent trials in adaptive planning
SZ_DENOM_RETRIES = 32

# --- Symbolic fallback ------------------------------------------------------
SIMPLIFY_TIMEOUT_SECONDS = 3.0
FUNCTION_EXPAND_TIMEOUT = 1.5
COMPLEXITY_THRESHOLD = 800.0  # gate potentially expensive symbolic expansion

# --- Caches ------------------------------------------------------------------
CANONICALIZE_CACHE_SIZE = 4096
RESULT_CACHE_SIZE = 2048
EXACT_MINPOLY_CACHE_SIZE = 1024
EXACT_CYCLO_CACHE_SIZE = 256
EXACT_ROOT_CACHE_SIZE = 512
GENERATOR_CACHE_SIZE = 1024
EXACT_BOUND_CACHE_SIZE = 512

# --- Exact algebraic-number reduction ---------------------------------------
# Algebraic modeling is intentionally budgeted: primitive-element degree can
# grow as a product of generator degrees, even for expressions that print
# compactly.
ALGEBRAIC_MAX_GENERATORS = 16
ALG_MODEL_MAX_DEGREE = 4096
ALG_COMMON_MAX_DEGREE = 512

ALG_REDUCE_MAX_DEGREE = 1024
ALG_TOWER_MAX_GENS = 24

# Algebraic reduction growth controls.
ALG_RESULTANT_MAX_DEGREE = 128
ALG_RESULTANT_MAX_OPS = 2500

# Exact transcendental normalization.
LOG_NORMALIZE_PASSES = 4

# Exact closed-constant reductions.
ALG_CANON_MAX_DEGREE = 1024
TRIG_FIELD_MAX_DEGREE = 256

# Exact cyclotomic arithmetic.
CYCLOTOMIC_MAX_ORDER = 4096
CYCLOTOMIC_MAX_DEGREE = 1024

# --- Fast exact-oracle budgets ---------------------------------------------
# These are checked before operations such as minimal-polynomial construction,
# primitive-element conversion, resultants, and large cyclotomic reductions.
EXACT_MAX_OPS = 160
EXACT_MAX_NODES = 320
EXACT_MAX_DEPTH = 28
EXACT_MAX_GENERATORS = 8
EXACT_MAX_DEGREE_PRODUCT = 256
MAX_CYCLOTOMIC_ORDER = 1024
EXACT_MAX_LOG_TERMS = 12
EXACT_STAGE_TIMEOUT = 0.20
EXACT_MAX_POW_BITS = 12
EXACT_MAX_INT_BITS = 512
EXACT_MAX_POLY_TERMS = 512
EXACT_MAX_RESULTANT_VARS = 5
EXACT_MAX_PSLQ_TERMS = 8
EXACT_PSLQ_PREC_DIGITS = 90
EXACT_PSLQ_MAX_COEFF = 1000000
EXACT_MAX_SPECIAL_DEPTH = 20
EXACT_MAX_FACTOR_BITS = 96

# --- Fast-stage profiling / applicability caches ---------------------------
NEGATIVE_CACHE_SIZE = 2048

# Tower arithmetic gets a smaller gate than generic algebraic work because it
# is intended to stay on the fast oracle path.
EXACT_TOWER_MAX_OPS = 100
TOWER_MAX_GENS = 6
TOWER_SIGN_TIMEOUT = 0.08

# Certified dyadic root-radius refinement.  More steps give tighter algebraic
# separation bounds but increase integer arithmetic; keep this deliberately
# small on the fast path.
ALG_GAP_DYADIC_STEPS = 10
EXACT_SPECIAL_MAX_INDEX = 256

# --- Shared fast structural metadata ---------------------------------------
EXACT_FEATURE_CACHE_SIZE = 4096
EXACT_NONZERO_CACHE_SIZE = 4096
NORMAL_FORM_CACHE_SIZE = 2048
QUICK_KIND_MAX_NODES = 48

# Adaptive rigorous numerics. A certified separation hint is converted into
# a starting precision but is capped so a bad estimate cannot make a cheap
# evaluation start at an excessive precision.
ARB_HINT_MAX_START_BITS = 2048
ARB_HINT_MIN_GUARD_BITS = 24

# --- Fast theorem and exact-rewrite limits ----------------------------------
EXP_INDEP_MAX_TERMS = 6
EXP_INDEP_GAP_MAX_OPS = 48
EXP_INDEP_GAP_TIMEOUT = 0.03
SQRT_SUM_MAX_TERMS = 12
SQRT_SUM_RAD_BITS = 32
RADICAL_DENEST_MAX_OPS = 24
EXACT_REWRITE_CACHE_SIZE = 2048
EXACT_REWRITE_MAX_OPS = 80
CALL_MEMO_MIN_NODES = 24
TOWER_SPARSE_MAX_TERMS = 384
ARB_MIN_STEP_BITS = 48
ARB_MAX_GROWTH_FACTOR = 4
ARB_TARGET_GUARD_BITS = 32

# --- Small identity normalizers ---------------------------------------------
# These limits keep polynomial/rational canonicalization and exact elementary
# rewrites firmly on the latency-sensitive path.
IDENTITY_RAT_MAX_OPS = 64
IDENTITY_RAT_MAX_NODES = 128
IDENTITY_RAT_MAX_VARS = 6
IDENTITY_RAT_MAX_ADD_TERMS = 32
IDENTITY_RAT_MAX_POW_EXP = 16
IDENTITY_RAT_MAX_RESULT_OPS = 96
IDENTITY_RAT_MAX_RESULT_NODES = 192
IDENTITY_ELEM_MAX_OPS = 80
IDENTITY_ELEM_MAX_NODES = 160
