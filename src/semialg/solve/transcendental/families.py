from __future__ import annotations

from collections.abc import Callable, Iterable, Sequence
from dataclasses import dataclass, field

import sympy as sp


@dataclass(frozen=True)
class TransFamDetection:
    family_name: str
    matches: tuple[sp.Expr, ...]
    metadata: dict = field(default_factory=dict)


@dataclass(frozen=True)
class TransFamHandler:
    family_name: str
    detector: Callable[[sp.Expr], tuple[sp.Expr, ...]]
    rewrite_builder: Callable[[sp.Expr, Sequence[sp.Symbol]], sp.Expr | None] | None = None


def _sorted_matches(exprs: Iterable[sp.Expr]) -> tuple[sp.Expr, ...]:
    return tuple(sorted(set(exprs), key=sp.default_sort_key))


def detect_trigonometric_fam(expr: sp.Expr) -> tuple[sp.Expr, ...]:
    return _sorted_matches(expr.atoms(sp.sin, sp.cos, sp.tan, sp.cot, sp.sec, sp.csc))


def detect_hyperbolic_family(expr: sp.Expr) -> tuple[sp.Expr, ...]:
    return _sorted_matches(expr.atoms(sp.sinh, sp.cosh, sp.tanh, sp.coth, sp.sech, sp.csch))


def detect_exponential_fam(expr: sp.Expr) -> tuple[sp.Expr, ...]:
    out = list(expr.atoms(sp.exp))
    out.extend([p for p in expr.atoms(sp.Pow) if not p.exp.is_rational and not p.base.is_number])
    return _sorted_matches(out)


def detect_inverse_family(expr: sp.Expr) -> tuple[sp.Expr, ...]:
    return _sorted_matches(
        expr.atoms(sp.log, sp.asin, sp.acos, sp.atan, sp.acot, sp.asinh, sp.acosh, sp.atanh)
    )


def detect_productlog_family(expr: sp.Expr) -> tuple[sp.Expr, ...]:
    return _sorted_matches(expr.atoms(sp.LambertW))


def detect_statistical_fam(expr: sp.Expr) -> tuple[sp.Expr, ...]:
    atoms = []
    for cls in [sp.erf, sp.erfc]:
        atoms.extend(expr.atoms(cls))
    return _sorted_matches(atoms)


def detect_arg_family(expr: sp.Expr) -> tuple[sp.Expr, ...]:
    return _sorted_matches(expr.atoms(sp.arg))


def detect_spec_function_fam(expr: sp.Expr) -> tuple[sp.Expr, ...]:
    atoms = []
    for cls in [sp.gamma, sp.zeta, sp.besselj, sp.bessely, sp.airyai, sp.airybi]:
        atoms.extend(expr.atoms(cls))
    return _sorted_matches(atoms)


def classify_trans_fams(expr: sp.Expr) -> tuple[TransFamDetection, ...]:
    families = [
        ("trigonometric", detect_trigonometric_fam(expr)),
        ("hyperbolic", detect_hyperbolic_family(expr)),
        ("exponential", detect_exponential_fam(expr)),
        ("inverse", detect_inverse_family(expr)),
        ("productlog", detect_productlog_family(expr)),
        ("statistical", detect_statistical_fam(expr)),
        ("arg", detect_arg_family(expr)),
        ("special", detect_spec_function_fam(expr)),
    ]
    return tuple(
        TransFamDetection(name, matches, {"count": len(matches)})
        for name, matches in families
        if matches
    )


def _rewrite_simple_equality(expr: sp.Expr, variables: Sequence[sp.Symbol]) -> sp.Expr | None:
    if not isinstance(expr, sp.Equality) or len(variables) != 1:
        return None
    x = variables[0]
    lhs = sp.simplify(expr.lhs - expr.rhs)
    expo = list(lhs.atoms(sp.exp))
    if len(expo) != 1:
        return None
    e = expo[0]
    rest = sp.simplify(lhs - e)
    try:
        sol = sp.solve(sp.Eq(e, -rest), x)
    except Exception:
        return None
    if len(sol) == 1:
        return sp.Eq(x, sp.simplify(sol[0]))
    return None


def _rewrite_simple_equalit2(expr: sp.Expr, variables: Sequence[sp.Symbol]) -> sp.Expr | None:
    if not isinstance(expr, sp.Equality) or len(variables) != 1:
        return None
    x = variables[0]
    logs = list((expr.lhs - expr.rhs).atoms(sp.log))
    if len(logs) != 1:
        return None
    try:
        sol = sp.solve(expr, x)
    except Exception:
        return None
    if len(sol) == 1:
        return sp.Eq(x, sp.simplify(sol[0]))
    return None


def _rewrite_trig_equal_zero(expr: sp.Expr, variables: Sequence[sp.Symbol]) -> sp.Expr | None:
    if not isinstance(expr, sp.Equality) or len(variables) != 1:
        return None
    x = variables[0]
    diff = sp.simplify(expr.lhs - expr.rhs)
    k = sp.Symbol(f"k_{x.name}", integer=True)
    if diff == sp.sin(x):
        return sp.Exists(k, sp.Eq(x, sp.pi * k))
    if diff == sp.cos(x):
        return sp.Exists(k, sp.Eq(x, sp.pi / 2 + sp.pi * k))
    if diff == sp.tan(x):
        return sp.Exists(k, sp.Eq(x, sp.pi * k))
    return None


def _rewrite_productlog_simp(expr: sp.Expr, variables: Sequence[sp.Symbol]) -> sp.Expr | None:
    if not isinstance(expr, sp.Equality) or len(variables) != 1:
        return None
    x = variables[0]
    try:
        sols = sp.solve(expr, x)
    except Exception:
        return None
    if len(sols) == 1:
        return sp.Eq(x, sp.simplify(sols[0]))
    return None


def _rewrite_erf_simple(expr: sp.Expr, variables: Sequence[sp.Symbol]) -> sp.Expr | None:
    if not isinstance(expr, sp.Equality) or len(variables) != 1:
        return None
    x = variables[0]
    diff = sp.simplify(expr.lhs - expr.rhs)
    erfs = list(diff.atoms(sp.erf))
    if len(erfs) != 1:
        return None
    target = erfs[0]
    rest = sp.simplify(diff - target)
    if target == sp.erf(x):
        return sp.Eq(x, sp.erfinv(-rest))
    return None


def _rewrite_erfc_simple(expr: sp.Expr, variables: Sequence[sp.Symbol]) -> sp.Expr | None:
    if not isinstance(expr, sp.Equality) or len(variables) != 1:
        return None
    x = variables[0]
    diff = sp.simplify(expr.lhs - expr.rhs)
    erfc_terms = list(diff.atoms(sp.erfc))
    if len(erfc_terms) != 1:
        return None
    target = erfc_terms[0]
    rest = sp.simplify(diff - target)
    if target == sp.erfc(x):
        return sp.Eq(x, sp.erfcinv(-rest))
    return None


def _rewrite_arg_real(expr: sp.Expr, variables: Sequence[sp.Symbol]) -> sp.Expr | None:
    if not isinstance(expr, sp.Equality) or len(variables) != 1:
        return None
    x = variables[0]
    if expr.lhs == sp.arg(x) and expr.rhs == 0:
        return x > 0
    return None


def default_trans_handlers() -> tuple[TransFamHandler, ...]:
    return (
        TransFamHandler("trigonometric", detect_trigonometric_fam, _rewrite_trig_equal_zero),
        TransFamHandler("hyperbolic", detect_hyperbolic_family, None),
        TransFamHandler("exponential", detect_exponential_fam, _rewrite_simple_equality),
        TransFamHandler("inverse", detect_inverse_family, _rewrite_simple_equalit2),
        TransFamHandler("productlog", detect_productlog_family, _rewrite_productlog_simp),
        TransFamHandler("statistical", detect_statistical_fam, _rewrite_erf_simple),
        TransFamHandler("arg", detect_arg_family, _rewrite_arg_real),
        TransFamHandler("special", detect_spec_function_fam, None),
    )


__all__ = [
    "TransFamDetection",
    "TransFamHandler",
    "detect_trigonometric_fam",
    "detect_hyperbolic_family",
    "detect_exponential_fam",
    "detect_inverse_family",
    "detect_productlog_family",
    "detect_statistical_fam",
    "detect_arg_family",
    "detect_spec_function_fam",
    "classify_trans_fams",
    "default_trans_handlers",
]
