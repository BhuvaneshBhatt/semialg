"""DLMF-audited real-property facts for parameterized special functions.

The rules here are intentionally conservative.  Each fact has provenance:
``direct`` means DLMF states the property under essentially the same
hypotheses; ``derived`` means funcprops combines DLMF identities, endpoint
limits, or elementary differentiation without weakening the hypotheses.
"""

from __future__ import annotations

from dataclasses import dataclass

import sympy as sp
from sympy import S

from ..exactness import certify_identity_equal
from ..exceptions import SYMBOLIC_QUERY_ERRORS
from ..stratification import semantic_condition


@dataclass(frozen=True)
class DLMFFact:
    family: str
    property: str
    section: str
    url: str
    status: str = "direct"  # direct | derived
    note: str = ""


DLMF_FACTS = (
    DLMFFact(
        "erf",
        "odd/strictly increasing/range",
        "7.2,7.10",
        "https://dlmf.nist.gov/7",
        note="real x; derivative 2 exp(-x^2)/sqrt(pi)>0",
    ),
    DLMFFact(
        "erfc",
        "positive/strictly decreasing/range",
        "7.2,7.10",
        "https://dlmf.nist.gov/7",
        note="real x; derivative -2 exp(-x^2)/sqrt(pi)<0",
    ),
    DLMFFact(
        "polygamma",
        "sign/monotonicity",
        "5.15",
        "https://dlmf.nist.gov/5.15",
        note="positive real argument with order prerequisites",
    ),
    DLMFFact(
        "airyai",
        "real zeros/oscillation",
        "9.9",
        "https://dlmf.nist.gov/9.9",
        note="real zeros are negative and simple",
    ),
    DLMFFact(
        "airybi",
        "real zeros/oscillation",
        "9.9",
        "https://dlmf.nist.gov/9.9",
        note="real zeros are negative and simple",
    ),
    DLMFFact(
        "besseli",
        "positive/increasing",
        "10.37",
        "https://dlmf.nist.gov/10.37",
        note="fixed nu>=0, x>0",
    ),
    DLMFFact(
        "besselk",
        "positive/decreasing",
        "10.37",
        "https://dlmf.nist.gov/10.37",
        note="fixed nu>=0, x>0",
    ),
    DLMFFact(
        "besselk",
        "order symmetry",
        "10.27.3",
        "https://dlmf.nist.gov/10.27.E3",
        note="K_-nu(z)=K_nu(z)",
    ),
    DLMFFact(
        "lowergamma",
        "integral/domain",
        "8.2.1",
        "https://dlmf.nist.gov/8.2.E1",
        note="Re(a)>0",
    ),
    DLMFFact("uppergamma", "integral/domain", "8.2.2", "https://dlmf.nist.gov/8.2.E2"),
    DLMFFact(
        "beta",
        "Euler integral/domain",
        "5.12.1",
        "https://dlmf.nist.gov/5.12.E1",
        note="Re(a)>0, Re(b)>0",
    ),
    DLMFFact(
        "betainc",
        "integral/domain/range",
        "8.17.1",
        "https://dlmf.nist.gov/8.17.E1",
        note="a>0, b>0, 0<=x<=1",
    ),
    DLMFFact(
        "betainc_regularized",
        "normalization/domain/range",
        "8.17.2",
        "https://dlmf.nist.gov/8.17.E2",
        note="a>0, b>0, 0<=x<=1",
    ),
    DLMFFact("gamma", "meromorphic/no zeros/poles", "5.2", "https://dlmf.nist.gov/5.2"),
    DLMFFact("digamma", "meromorphic/poles", "5.2.2", "https://dlmf.nist.gov/5.2.E2"),
    DLMFFact(
        "hyper",
        "series/termination/principal branch",
        "16.2.1",
        "https://dlmf.nist.gov/16.2.E1",
        note="p<=q entire; p=q+1 principal branch cut [1,oo) unless terminating",
    ),
    DLMFFact(
        "assoc_legendre",
        "Ferrers zero distribution",
        "14.16",
        "https://dlmf.nist.gov/14.16",
        note="integer-degree/order specialization gives n-|m| interior zeros on (-1,1)",
    ),
    DLMFFact(
        "assoc_legendre",
        "principal complex branch",
        "14.21",
        "https://dlmf.nist.gov/14.21",
        note="principal branches use the standard cut to 1, with polynomial/degenerate exceptions handled separately",
    ),
    DLMFFact(
        "elliptic_f",
        "incomplete integral branch radicand",
        "19.2.4",
        "https://dlmf.nist.gov/19.2.E4",
        note="principal sqrt of 1-k^2 sin^2(theta)",
    ),
    DLMFFact(
        "elliptic_e",
        "incomplete integral branch radicand",
        "19.2.5",
        "https://dlmf.nist.gov/19.2.E5",
        note="principal sqrt of 1-k^2 sin^2(theta)",
    ),
    DLMFFact(
        "elliptic_pi",
        "incomplete third-kind path singularity",
        "19.2.7",
        "https://dlmf.nist.gov/19.2.E7",
        note="parameter singularity requires path-aware treatment; endpoint-only tests are insufficient",
    ),
    DLMFFact(
        "jacobi",
        "zeros/interlacing",
        "18.16.1",
        "https://dlmf.nist.gov/18.16.E1",
        note="n simple zeros in (-1,1) for alpha,beta>-1; consecutive degrees interlace by 18.2(vi)",
    ),
    DLMFFact(
        "assoc_laguerre",
        "zeros/interlacing",
        "18.16.9",
        "https://dlmf.nist.gov/18.16.E9",
        note="n simple positive zeros for alpha>-1",
    ),
    DLMFFact(
        "chebyshevt",
        "trigonometric representation/zeros",
        "18.5.1",
        "https://dlmf.nist.gov/18.5.E1",
    ),
    DLMFFact(
        "chebyshevu",
        "trigonometric representation/zeros",
        "18.5.2",
        "https://dlmf.nist.gov/18.5.E2",
    ),
    DLMFFact(
        "jn",
        "positive zeros",
        "10.58.1",
        "https://dlmf.nist.gov/10.58.E1",
        note="spherical Bessel zeros inherit ordinary Bessel zero structure",
    ),
    DLMFFact(
        "yn",
        "positive zeros",
        "10.58.1",
        "https://dlmf.nist.gov/10.58.E1",
        note="spherical Bessel zeros inherit ordinary Bessel zero structure",
    ),
    DLMFFact(
        "multigamma",
        "domain/product",
        "35.3.1,35.3.4",
        "https://dlmf.nist.gov/35.3.E1",
        note="Re(a)>(p-1)/2; product formula 35.3.E4",
    ),
)


def dlmf_facts(family: str | None = None):
    """Return executable and formula-audit DLMF facts.

    The generated audit contributes one provenance record per exact formula or
    narrowly scoped DLMF section.  Helper/wrapper families are marked
    ``derived``.
    """
    facts = list(DLMF_FACTS)
    try:
        import json
        from importlib.resources import files

        data = json.loads(
            files("funcprops")
            .joinpath("data/dlmf_audit.json")
            .read_text(encoding="utf-8")
        )
        seen = {(f.family, f.section) for f in facts}
        for row in data.get("families", ()):
            if row.get("status") not in {
                "formula-verified-core",
                "formula-verified-derived",
            }:
                continue
            status = (
                "derived" if row.get("verification_kind") == "derived" else "direct"
            )
            for ref in row.get("formula_refs", ()):
                key = (row["family_id"], ref)
                if key in seen:
                    continue
                facts.append(
                    DLMFFact(
                        row["family_id"],
                        row.get("formula_verified_properties", "formula audit"),
                        ref,
                        f"https://dlmf.nist.gov/{ref}",
                        status=status,
                        note="formula-level catalogue audit",
                    )
                )
                seen.add(key)
    except SYMBOLIC_QUERY_ERRORS:
        pass
    return tuple(f for f in facts if family is None or f.family == family)


def _real_guard(a):
    return S.true if a.is_real is True else sp.Contains(a, S.Reals)


def _besseli_successive_ratio_order(expr, var):
    """Return nu for I_{nu+1}(x)/I_nu(x), else None."""
    expr = sp.sympify(expr)
    var = sp.sympify(var)
    num, den = sp.fraction(expr)
    if (
        num.func is sp.besseli
        and den.func is sp.besseli
        and len(num.args) == 2
        and len(den.args) == 2
        and num.args[1] == var
        and den.args[1] == var
        and certify_identity_equal(num.args[0], den.args[0] + 1) is True
    ):
        return den.args[0]
    return None


def _positive_axis(domain, var):
    try:
        from .common import condition_to_set

        ds = condition_to_set(domain, var)
        if (
            isinstance(ds, sp.Interval)
            and ds.start == 0
            and ds.left_open
            and ds.end is S.Infinity
        ):
            return True
    except SYMBOLIC_QUERY_ERRORS:
        pass
    try:
        return sp.simplify_logic(sp.Implies(domain, var > 0)) is S.true
    except SYMBOLIC_QUERY_ERRORS:
        return False


def parameter_guard(expr, var, prop):
    expr = sp.sympify(expr)
    var = sp.sympify(var)
    if expr.func is sp.besselk and len(expr.args) == 2 and expr.args[1] == var:
        nu = expr.args[0]
        if prop in {
            "Positive",
            "NonNegative",
            "Range",
            "GlobalMonotonicity",
            "GlobalConvexity",
        }:
            # 10.37 gives nu>=0; 10.27.3 extends K to all real order.
            return _real_guard(nu)
    if expr.func is sp.besseli and len(expr.args) == 2 and expr.args[1] == var:
        nu = expr.args[0]
        if prop in {"Positive", "NonNegative", "Range", "GlobalMonotonicity"}:
            return sp.And(_real_guard(nu), nu >= 0)
    return S.true


def positive_axis_domain_strata(expr, var):
    """Sufficient exact real-domain strata on x>0.

    These are not used to claim maximal domains outside the stated guards;
    they provide parameter-conditioned branches that can be composed safely.
    """
    expr = sp.sympify(expr)
    var = sp.sympify(var)
    if (
        expr.func in (sp.lowergamma, sp.uppergamma)
        and len(expr.args) == 2
        and expr.args[1] == var
    ):
        a = expr.args[0]
        return ((var > 0, a > 0, "DLMF 8.2.1-8.2.2"),)
    from ..statistical_special_functions import RegularizedGammaP, RegularizedGammaQ

    if (
        expr.func in (RegularizedGammaP, RegularizedGammaQ)
        and len(expr.args) == 2
        and expr.args[1] == var
    ):
        a = expr.args[0]
        return ((var >= 0, a > 0, "regularized incomplete gamma integral"),)
    if expr.func is sp.beta and len(expr.args) == 2 and expr.args[1] == var:
        a = expr.args[0]
        return ((var > 0, a > 0, "DLMF 5.12.1 + Gamma quotient"),)
    if expr.func is sp.besseli and len(expr.args) == 2 and expr.args[1] == var:
        nu = expr.args[0]
        return ((var > 0, sp.And(_real_guard(nu), nu >= 0), "DLMF 10.37"),)
    if expr.func is sp.besselk and len(expr.args) == 2 and expr.args[1] == var:
        nu = expr.args[0]
        return ((var > 0, _real_guard(nu), "DLMF 10.37 + 10.27.3"),)
    if (
        expr.func in (sp.betainc, sp.betainc_regularized)
        and len(expr.args) == 4
        and expr.args[2] == 0
        and expr.args[3] == var
    ):
        a, b = expr.args[:2]
        return (
            (sp.And(var >= 0, var <= 1), sp.And(a > 0, b > 0), "DLMF 8.17.1-8.17.2"),
        )
    return ()


def positive_axis_range_strata(expr, var, y):
    """Exact ranges on x>0, including parameter-dependent endpoints."""
    expr = sp.sympify(expr)
    var = sp.sympify(var)
    y = sp.sympify(y)
    nu_ratio = _besseli_successive_ratio_order(expr, var)
    if nu_ratio is not None:
        return (
            (
                sp.And(y > 0, y < 1),
                sp.And(_real_guard(nu_ratio), nu_ratio >= -sp.Rational(1, 2)),
                "classical modified-Bessel quotient theorem",
            ),
        )
    if expr.func is sp.besselk and len(expr.args) == 2 and expr.args[1] == var:
        nu = expr.args[0]
        return ((y > 0, _real_guard(nu), "DLMF 10.37 + 10.27.3 + endpoint limits"),)
    if expr.func is sp.besseli and len(expr.args) == 2 and expr.args[1] == var:
        nu = expr.args[0]
        rg = _real_guard(nu)
        # I_0(x) -> 1 as x->0+, while I_nu(x)->0 for nu>0; all diverge to
        # +oo at +oo.  The endpoint is not attained because x>0.
        return (
            (y > 1, sp.And(rg, sp.Eq(nu, 0)), "DLMF 10.37 + 10.30 limiting forms"),
            (y > 0, sp.And(rg, nu > 0), "DLMF 10.37 + 10.30 limiting forms"),
        )
    from ..statistical_special_functions import RegularizedGammaP, RegularizedGammaQ

    if (
        expr.func in (RegularizedGammaP, RegularizedGammaQ)
        and len(expr.args) == 2
        and expr.args[1] == var
    ):
        a = expr.args[0]
        return (
            (
                sp.And(y > 0, y < 1),
                a > 0,
                "regularized incomplete gamma endpoint limits",
            ),
        )
    if expr.func is sp.beta and len(expr.args) == 2 and expr.args[1] == var:
        a = expr.args[0]
        return ((y > 0, a > 0, "DLMF 5.12.1 + Gamma quotient limits"),)
    if expr.func is sp.lowergamma and len(expr.args) == 2 and expr.args[1] == var:
        a = expr.args[0]
        return (
            (sp.And(y > 0, y < sp.gamma(a)), a > 0, "DLMF 8.2.1 + endpoint limits"),
        )
    if expr.func is sp.uppergamma and len(expr.args) == 2 and expr.args[1] == var:
        a = expr.args[0]
        return (
            (sp.And(y > 0, y < sp.gamma(a)), a > 0, "DLMF 8.2.2 + endpoint limits"),
        )
    return ()


def special_function_range_strata(expr, var, y, domain):
    """Exact DLMF-backed ranges with domain-sensitive endpoint attainment."""
    from .common import condition_to_set

    expr = sp.sympify(expr)
    var = sp.sympify(var)
    y = sp.sympify(y)
    ds = condition_to_set(domain, var)
    pos_open = (
        isinstance(ds, sp.Interval)
        and ds.start == 0
        and ds.left_open
        and ds.end is S.Infinity
    )
    pos_closed = (
        isinstance(ds, sp.Interval)
        and ds.start == 0
        and not ds.left_open
        and ds.end is S.Infinity
    )
    unit_open = (
        isinstance(ds, sp.Interval)
        and ds.start == 0
        and ds.end == 1
        and ds.left_open
        and ds.right_open
    )
    unit_closed = (
        isinstance(ds, sp.Interval)
        and ds.start == 0
        and ds.end == 1
        and not ds.left_open
        and not ds.right_open
    )

    from ..statistical_special_functions import RegularizedGammaP, RegularizedGammaQ

    if (
        expr.func in (RegularizedGammaP, RegularizedGammaQ)
        and len(expr.args) == 2
        and expr.args[1] == var
        and (pos_open or pos_closed)
    ):
        a = expr.args[0]
        if pos_open:
            cond = sp.And(y > 0, y < 1)
        elif expr.func is RegularizedGammaP:
            cond = sp.And(y >= 0, y < 1)
        else:
            cond = sp.And(y > 0, y <= 1)
        return ((cond, a > 0, "regularized incomplete gamma endpoint limits"),)

    if (
        expr.func is sp.marcumq
        and len(expr.args) == 3
        and expr.args[2] == var
        and (pos_open or pos_closed)
    ):
        m, a, _ = expr.args
        cond = sp.And(y > 0, y < 1) if pos_open else sp.And(y > 0, y <= 1)
        return ((cond, sp.And(m > 0, a >= 0), "Marcum Q survival-function range"),)

    if (
        expr.func is sp.besseli
        and len(expr.args) == 2
        and expr.args[1] == var
        and (pos_open or pos_closed)
    ):
        nu = expr.args[0]
        rg = _real_guard(nu)
        if pos_open:
            return (
                (y > 1, sp.And(rg, sp.Eq(nu, 0)), "DLMF 10.37, 10.30, 10.40"),
                (y > 0, sp.And(rg, nu > 0), "DLMF 10.37, 10.30, 10.40"),
            )
        return (
            (y >= 1, sp.And(rg, sp.Eq(nu, 0)), "DLMF 10.37, 10.30, 10.40"),
            (y >= 0, sp.And(rg, nu > 0), "DLMF 10.37, 10.30, 10.40"),
        )

    if (
        expr.func in (sp.lowergamma, sp.uppergamma)
        and len(expr.args) == 2
        and expr.args[1] == var
        and (pos_open or pos_closed)
    ):
        a = expr.args[0]
        G = sp.gamma(a)
        if expr.func is sp.lowergamma:
            cond = sp.And(y > 0, y < G) if pos_open else sp.And(y >= 0, y < G)
        else:
            cond = sp.And(y > 0, y < G) if pos_open else sp.And(y > 0, y <= G)
        return ((cond, a > 0, "DLMF 8.2.1-8.2.2 + endpoint limits"),)

    if (
        expr.func in (sp.betainc, sp.betainc_regularized)
        and len(expr.args) == 4
        and expr.args[2] == 0
        and expr.args[3] == var
        and (unit_open or unit_closed)
    ):
        a, b = expr.args[:2]
        guard = sp.And(a > 0, b > 0)
        endpoint = S.One if expr.func is sp.betainc_regularized else sp.beta(a, b)
        cond = (
            sp.And(y > 0, y < endpoint) if unit_open else sp.And(y >= 0, y <= endpoint)
        )
        return ((cond, guard, "DLMF 8.17.1-8.17.3"),)
    return ()


def positive_axis_monotonicity_strata(expr, var, domain=None):
    from ..results import Monotonicity
    from ..statistical_special_functions import MittagLefflerE, OwenT

    expr = sp.sympify(expr)
    var = sp.sympify(var)
    nu_ratio = _besseli_successive_ratio_order(expr, var)
    if nu_ratio is not None:
        return (
            (
                Monotonicity.STRICTLY_INCREASING,
                sp.And(_real_guard(nu_ratio), nu_ratio >= -sp.Rational(1, 2)),
                "classical modified-Bessel quotient theorem",
            ),
        )
    from ..statistical_special_functions import RegularizedGammaP, RegularizedGammaQ

    if expr.func is RegularizedGammaP and len(expr.args) == 2 and expr.args[1] == var:
        return (
            (
                Monotonicity.STRICTLY_INCREASING,
                expr.args[0] > 0,
                "regularized lower gamma derivative",
            ),
        )
    if expr.func is RegularizedGammaQ and len(expr.args) == 2 and expr.args[1] == var:
        return (
            (
                Monotonicity.STRICTLY_DECREASING,
                expr.args[0] > 0,
                "regularized upper gamma derivative",
            ),
        )
    if expr.func is sp.polygamma and len(expr.args) == 2 and expr.args[1] == var:
        n = expr.args[0]
        if n == 0:
            return ((Monotonicity.STRICTLY_INCREASING, S.true, "trigamma positivity"),)
        odd = sp.And(sp.Contains(n, S.Integers), n >= 1, sp.Eq(sp.Mod(n, 2), 1))
        even = sp.And(sp.Contains(n, S.Integers), n >= 1, sp.Eq(sp.Mod(n, 2), 0))
        return (
            (Monotonicity.STRICTLY_DECREASING, odd, "polygamma sign theorem"),
            (Monotonicity.STRICTLY_INCREASING, even, "polygamma sign theorem"),
        )
    if expr.func is sp.erf and expr.args == (var,):
        return ((Monotonicity.STRICTLY_INCREASING, S.true, "erf derivative"),)
    if expr.func is sp.erfc and expr.args == (var,):
        return ((Monotonicity.STRICTLY_DECREASING, S.true, "erfc derivative"),)
    if expr.func is OwenT and len(expr.args) == 2 and expr.args[1] == var:
        return ((Monotonicity.STRICTLY_INCREASING, S.true, "Owen T derivative"),)
    if expr.func is sp.marcumq and len(expr.args) == 3:
        m, a, b = expr.args
        if b == var:
            return (
                (
                    Monotonicity.STRICTLY_DECREASING,
                    sp.And(m > 0, a >= 0),
                    "Marcum Q threshold monotonicity",
                ),
            )
        if a == var:
            return (
                (
                    Monotonicity.STRICTLY_INCREASING,
                    sp.And(m > 0, b > 0),
                    "Marcum Q noncentrality monotonicity",
                ),
            )
    if expr.func is MittagLefflerE and len(expr.args) == 2 and expr.args[1] == -var:
        alpha = expr.args[0]
        return (
            (
                Monotonicity.NONINCREASING,
                sp.And(alpha > 0, alpha <= 1),
                "Pollard complete monotonicity",
            ),
        )
    if (
        expr.func in (sp.betainc, sp.betainc_regularized)
        and len(expr.args) == 4
        and expr.args[2] == 0
        and expr.args[3] == var
        and domain is not None
    ):
        a, b = expr.args[:2]
        try:
            on_unit = (
                sp.simplify_logic(
                    sp.Implies(semantic_condition(domain), sp.And(var > 0, var < 1))
                )
                is S.true
            )
        except SYMBOLIC_QUERY_ERRORS:
            on_unit = False
        if on_unit:
            return (
                (
                    Monotonicity.STRICTLY_INCREASING,
                    sp.And(a > 0, b > 0),
                    "incomplete beta integrand positivity",
                ),
            )
    if expr.func is sp.besselk and len(expr.args) == 2 and expr.args[1] == var:
        return (
            (
                Monotonicity.STRICTLY_DECREASING,
                _real_guard(expr.args[0]),
                "DLMF 10.37 + 10.27.3",
            ),
        )
    if expr.func is sp.besseli and len(expr.args) == 2 and expr.args[1] == var:
        nu = expr.args[0]
        return (
            (
                Monotonicity.STRICTLY_INCREASING,
                sp.And(_real_guard(nu), nu >= 0),
                "DLMF 10.37",
            ),
        )
    if expr.func is sp.beta and len(expr.args) == 2 and expr.args[1] == var:
        return (
            (
                Monotonicity.STRICTLY_DECREASING,
                expr.args[0] > 0,
                "derived from digamma monotonicity / beta quotient",
            ),
        )
    if expr.func is sp.lowergamma and len(expr.args) == 2 and expr.args[1] == var:
        return (
            (
                Monotonicity.STRICTLY_INCREASING,
                expr.args[0] > 0,
                "DLMF 8.2.1 integrand > 0",
            ),
        )
    if expr.func is sp.uppergamma and len(expr.args) == 2 and expr.args[1] == var:
        return (
            (
                Monotonicity.STRICTLY_DECREASING,
                expr.args[0] > 0,
                "DLMF 8.2.2 integrand > 0",
            ),
        )
    return ()


def positive_axis_sign_strata(expr, var):
    from ..results import Sign

    expr = sp.sympify(expr)
    var = sp.sympify(var)
    from ..statistical_special_functions import (
        MittagLefflerE,
        OwenT,
        RegularizedGammaP,
        RegularizedGammaQ,
    )

    nu_ratio = _besseli_successive_ratio_order(expr, var)
    if nu_ratio is not None:
        return (
            (
                Sign.POSITIVE,
                sp.And(_real_guard(nu_ratio), nu_ratio >= -sp.Rational(1, 2)),
                "classical modified-Bessel quotient theorem",
            ),
        )

    if (
        expr.func in (RegularizedGammaP, RegularizedGammaQ)
        and len(expr.args) == 2
        and expr.args[1] == var
    ):
        return (
            (
                Sign.POSITIVE,
                expr.args[0] > 0,
                "regularized incomplete gamma positivity",
            ),
        )
    if expr.func is OwenT and len(expr.args) == 2 and expr.args[1] == var:
        return (
            (Sign.POSITIVE, _real_guard(expr.args[0]), "Owen T integral positivity"),
        )
    if expr.func is sp.marcumq and len(expr.args) == 3 and expr.args[2] == var:
        m, a, _ = expr.args
        return (
            (
                Sign.POSITIVE,
                sp.And(m > 0, a >= 0),
                "Marcum Q survival-function positivity",
            ),
        )
    if expr.func is sp.polygamma and len(expr.args) == 2 and expr.args[1] == var:
        n = expr.args[0]
        pos = sp.And(sp.Contains(n, S.Integers), n >= 1, sp.Eq(sp.Mod(n, 2), 1))
        neg = sp.And(sp.Contains(n, S.Integers), n >= 1, sp.Eq(sp.Mod(n, 2), 0))
        return (
            (Sign.POSITIVE, pos, "polygamma series"),
            (Sign.NEGATIVE, neg, "polygamma series"),
        )
    if expr.func is MittagLefflerE and len(expr.args) == 2 and expr.args[1] == -var:
        alpha = expr.args[0]
        return (
            (
                Sign.POSITIVE,
                sp.And(alpha > 0, alpha <= 1),
                "Pollard complete monotonicity",
            ),
        )
    if expr.func is sp.besselk and len(expr.args) == 2 and expr.args[1] == var:
        return ((Sign.POSITIVE, _real_guard(expr.args[0]), "DLMF 10.37 + 10.27.3"),)
    if expr.func is sp.besseli and len(expr.args) == 2 and expr.args[1] == var:
        nu = expr.args[0]
        return ((Sign.POSITIVE, sp.And(_real_guard(nu), nu >= 0), "DLMF 10.37"),)
    if expr.func is sp.beta and len(expr.args) == 2 and expr.args[1] == var:
        return ((Sign.POSITIVE, expr.args[0] > 0, "DLMF 5.12.1"),)
    if (
        expr.func in (sp.lowergamma, sp.uppergamma)
        and len(expr.args) == 2
        and expr.args[1] == var
    ):
        return ((Sign.POSITIVE, expr.args[0] > 0, "DLMF 8.2.1-8.2.2"),)
    return ()


def positive_axis_convexity_strata(expr, var):
    from ..results import Convexity

    expr = sp.sympify(expr)
    var = sp.sympify(var)
    from ..statistical_special_functions import MittagLefflerE

    nu_ratio = _besseli_successive_ratio_order(expr, var)
    if nu_ratio is not None:
        return (
            (
                Convexity.CONCAVE,
                sp.And(_real_guard(nu_ratio), nu_ratio >= -sp.Rational(1, 2)),
                "classical modified-Bessel quotient theorem",
            ),
        )
    if expr.func is sp.loggamma and len(expr.args) == 1 and expr.args[0] == var:
        return ((Convexity.STRICTLY_CONVEX, S.true, "trigamma positivity"),)
    if (
        expr.func is sp.polygamma
        and len(expr.args) == 2
        and expr.args[1] == var
        and expr.args[0] == 0
    ):
        return ((Convexity.STRICTLY_CONCAVE, S.true, "polygamma sign theorem"),)
    if expr.func is MittagLefflerE and len(expr.args) == 2 and expr.args[1] == -var:
        alpha = expr.args[0]
        return (
            (Convexity.CONVEX, sp.And(alpha > 0, alpha <= 1), "complete monotonicity"),
        )
    if expr.func is sp.besselk and len(expr.args) == 2 and expr.args[1] == var:
        return (
            (
                Convexity.STRICTLY_CONVEX,
                _real_guard(expr.args[0]),
                "DLMF-backed integral/complete monotonicity",
            ),
        )
    if expr.func is sp.beta and len(expr.args) == 2 and expr.args[1] == var:
        return (
            (Convexity.STRICTLY_CONVEX, expr.args[0] > 0, "derived from beta quotient"),
        )
    return ()


# Additional theorem-oriented catalogue facts used by cross-property reasoning.
EXPANDED_SPECIAL_FUNCTION_FACTS = (
    DLMFFact(
        "LambertW",
        "principal branch positive/increasing",
        "4.13",
        "https://dlmf.nist.gov/4.13",
        note="real x>0 on W_0; derivative is positive",
    ),
    DLMFFact(
        "zeta",
        "positive/decreasing",
        "25.2",
        "https://dlmf.nist.gov/25.2",
        note="real s>1 from the Dirichlet series and its differentiated series",
    ),
    DLMFFact(
        "struveh",
        "differential-equation/real-axis structure",
        "11.2",
        "https://dlmf.nist.gov/11.2",
        status="derived",
        note="catalogue hook; property claims require guarded theorem records",
    ),
    DLMFFact(
        "polylog",
        "principal branch/series domain",
        "25.12",
        "https://dlmf.nist.gov/25.12",
        note="series facts are used only inside their convergence domain",
    ),
)
DLMF_FACTS = DLMF_FACTS + EXPANDED_SPECIAL_FUNCTION_FACTS
