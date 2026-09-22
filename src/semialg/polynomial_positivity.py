"""Exact specialized decision procedures for polynomial negativity/nonnegativity.

The backend implements the critical-value core of the Zeng decision strategy:
cheap exact witnesses are tried first, then coercive polynomials are reduced to
an exact zero-dimensional gradient system.  A result is marked complete only
when these hypotheses prove that every global minimum is represented by the
critical system. Unsupported cases remain explicitly incomplete.

References: Zeng--Zeng (2004) [ZZ2004] motivates the ``zeng``
naming; the present attained-minimum/critical-value route is closer in spirit
to Zeng--Xiao (2012) [ZX2012].  See ``docs/references.md``.  This module does
not claim a verbatim implementation of either paper.
"""

from __future__ import annotations

from collections.abc import Callable, Sequence
from dataclasses import dataclass

import sympy as sp

from ._linear_relations import certified_sign
from .algebraic.rational_univariate import RationalUnivariateError
from .decision_portfolio import CertifiedDecisionResult
from .solve.zero_dimensional import is_zero_dimensional, solve_zero_dimensional_system
from .witness_heuristics import _certified_negative_point, find_negative_witness_fast


@dataclass(frozen=True)
class PolynomialNegativityResult:
    variables: tuple[sp.Symbol, ...]
    has_negative_point: bool | None
    point: tuple[tuple[sp.Symbol, sp.Expr], ...] | None
    complete: bool
    method: str
    notes: tuple[str, ...] = ()

    @property
    def assignment(self) -> dict[sp.Symbol, sp.Expr] | None:
        return dict(self.point) if self.point is not None else None

    @property
    def nonnegative(self) -> bool | None:
        if not self.complete or self.has_negative_point is None:
            return None
        return not self.has_negative_point

    def require_decision(self) -> bool:
        if not self.complete or self.has_negative_point is None:
            raise NotImplementedError("polynomial negativity decision is incomplete")
        return self.has_negative_point


@dataclass(frozen=True)
class _CoercivityCertificate:
    certified: bool
    method: str
    leading_form: sp.Expr
    margin: sp.Rational | None = None
    notes: tuple[str, ...] = ()


def _leading_homogeneous_form(
    polynomial: sp.Expr, variables: tuple[sp.Symbol, ...]
) -> tuple[sp.Expr, int]:
    poly = sp.Poly(sp.expand(polynomial), *variables, domain=sp.QQ)
    degree = int(poly.total_degree())
    leading = sp.Add(
        *(
            sp.Rational(coeff)
            * sp.prod(variable**power for variable, power in zip(variables, exponent, strict=True))
            for exponent, coeff in poly.terms()
            if sum(exponent) == degree
        )
    )
    return sp.expand(leading), degree


def _coercive_even_monomial_certificate(
    leading: sp.Expr, degree: int, variables: tuple[sp.Symbol, ...]
) -> bool:
    """Cheap sufficient positivity test retained from the original Zeng backend."""
    if degree <= 0 or degree % 2:
        return False
    poly = sp.Poly(leading, *variables, domain=sp.QQ)
    top = [(exp, sp.Rational(coeff)) for exp, coeff in poly.terms()]
    if not top or any(coeff < 0 for _, coeff in top):
        return False
    if any(any(int(power) % 2 for power in exponent) for exponent, _ in top):
        return False
    for index in range(len(variables)):
        pure = tuple(degree if j == index else 0 for j in range(len(variables)))
        if not any(exp == pure and coeff > 0 for exp, coeff in top):
            return False
    return True


def _obvious_nonpositive_leading_direction(
    leading: sp.Expr, variables: tuple[sp.Symbol, ...]
) -> tuple[sp.Expr, ...] | None:
    """Return a cheap exact nonzero direction where the form is nonpositive."""
    nvars = len(variables)
    directions: list[tuple[sp.Expr, ...]] = []
    for index in range(nvars):
        directions.append(tuple(sp.Integer(1 if j == index else 0) for j in range(nvars)))
    for left in range(nvars):
        for right in range(left + 1, nvars):
            for sign in (sp.Integer(1), sp.Integer(-1)):
                directions.append(
                    tuple(
                        sp.Integer(1) if j == left else sign if j == right else sp.Integer(0)
                        for j in range(nvars)
                    )
                )
    for direction in directions:
        value = sp.expand(leading.subs(dict(zip(variables, direction, strict=True))))
        sign = certified_sign(value)
        if sign in {-1, 0}:
            return direction
    return None


def _sos_positive_leading_form_certificate(
    leading: sp.Expr, degree: int, variables: tuple[sp.Symbol, ...], *, searcher=None
) -> _CoercivityCertificate | None:
    """Try to prove a strict radial SOS lower bound for a homogeneous form.

    A certificate for ``H - eps*(sum(x_i**2))**(d/2)`` being SOS, with
    exact rational ``eps > 0``, proves ``H(x) >= eps*||x||**d`` and hence
    positive definiteness of the leading form.  Merely proving ``H`` SOS would
    be insufficient because an SOS form may vanish on a nonzero direction.
    """
    if searcher is None:
        from .sos_certificates import search_sos_certificate as searcher

    radial = sp.expand(sum(variable**2 for variable in variables) ** (degree // 2))
    poly = sp.Poly(leading, *variables, domain=sp.QQ)
    scale = max((abs(sp.Rational(coeff)) for coeff in poly.coeffs()), default=sp.Integer(1))
    if scale == 0:
        return None
    margins = tuple(sp.Rational(scale, 2**power) for power in range(1, 9))
    for margin in margins:
        shifted = sp.expand(leading - margin * radial)
        searched = searcher(shifted, variables, backend="auto", use_planner=True)
        if searched.certified:
            return _CoercivityCertificate(
                True,
                "leading_form_sos_margin",
                leading,
                margin,
                (f"sos_margin={margin}",),
            )
    return None


def _sphere_positive_leading_form_certificate(
    leading: sp.Expr, degree: int, variables: tuple[sp.Symbol, ...]
) -> _CoercivityCertificate | None:
    """Decide positive definiteness exactly on the real unit sphere.

    For an even homogeneous form H, positivity on ``sum(x_i**2)=1`` is
    equivalent to positive definiteness on every nonzero real direction.
    Unsatisfiability of ``sphere AND H <= 0`` is therefore an exact coercivity
    certificate for any polynomial whose leading homogeneous form is H.
    """
    from .decision import is_satisfiable

    sphere = sp.Eq(sum(variable**2 for variable in variables), 1)
    nonpositive = sp.Le(leading, 0)
    try:
        result = is_satisfiable(
            sp.And(sphere, nonpositive),
            variables,
            strategy="cad",
            return_result=True,
        )
    except (ArithmeticError, TypeError, ValueError, NotImplementedError, sp.PolynomialError):
        return None
    if result.satisfiable:
        return _CoercivityCertificate(
            False,
            "leading_form_sphere",
            leading,
            notes=("nonpositive_direction_exists", f"sphere_method={result.method}"),
        )
    return _CoercivityCertificate(
        True,
        "leading_form_sphere",
        leading,
        notes=(f"sphere_method={result.method}",),
    )


def _coercivity_certificate(
    polynomial: sp.Expr, variables: tuple[sp.Symbol, ...], *, sos_searcher=None
) -> _CoercivityCertificate:
    """Certify coercivity from exact positivity of the leading homogeneous form."""
    leading, degree = _leading_homogeneous_form(polynomial, variables)
    if degree <= 0 or degree % 2:
        return _CoercivityCertificate(False, "leading_form_degree", leading)
    if _coercive_even_monomial_certificate(leading, degree, variables):
        return _CoercivityCertificate(True, "leading_form_even_monomials", leading)

    bad_direction = _obvious_nonpositive_leading_direction(leading, variables)
    if bad_direction is not None:
        rendered = ",".join(sp.sstr(value) for value in bad_direction)
        return _CoercivityCertificate(
            False,
            "leading_form_direction",
            leading,
            notes=(f"nonpositive_direction=({rendered})",),
        )

    sos = _sos_positive_leading_form_certificate(leading, degree, variables, searcher=sos_searcher)
    if sos is not None:
        return sos

    sphere = _sphere_positive_leading_form_certificate(leading, degree, variables)
    if sphere is not None:
        return sphere
    return _CoercivityCertificate(False, "leading_form_unresolved", leading)


def zeng_negative_point(
    polynomial: sp.Expr,
    variables: Sequence[sp.Symbol],
    *,
    random_lines: int = 8,
    seed: int = 1234,
) -> PolynomialNegativityResult:
    """Decide whether a rational polynomial is negative somewhere when certified.

    The implementation is a Zeng-family inspired exact critical-value reduction;
    see [ZZ2004] and [ZX2012] in ``docs/references.md``.  It is not a literal
    transcription of the 2004 semidefiniteness algorithm. On the certified
    coercive fragment, global minima are attained and the gradient ideal is
    solved exactly by the package's RUR backend.  Fast odd-degree/ray and
    random-line searches may prove the positive (negative-point) side earlier;
    they never prove nonnegativity.
    """
    vars_ = tuple(variables)
    if random_lines < 0:
        raise ValueError("random_lines must be nonnegative")
    if not vars_:
        sign = certified_sign(sp.sympify(polynomial))
        return PolynomialNegativityResult(
            vars_, sign == -1, tuple() if sign == -1 else None, sign is not None, "zeng_constant"
        )
    expr = sp.expand(polynomial)
    try:
        poly = sp.Poly(expr, *vars_, domain=sp.QQ)
    except (sp.PolynomialError, ValueError, TypeError) as exc:
        raise RationalUnivariateError("Zeng decision requires a rational polynomial") from exc

    if poly.is_zero:
        return PolynomialNegativityResult(vars_, False, None, True, "zeng", ("zero_polynomial",))
    constant_sign = certified_sign(expr) if not expr.free_symbols.intersection(vars_) else None
    if constant_sign is not None:
        point = tuple((v, sp.Integer(0)) for v in vars_) if constant_sign == -1 else None
        return PolynomialNegativityResult(vars_, constant_sign == -1, point, True, "zeng_constant")

    fast = find_negative_witness_fast(expr, vars_, random_lines=random_lines, seed=seed)
    if fast.found:
        return PolynomialNegativityResult(vars_, True, fast.point, True, f"zeng+{fast.method}")

    coercivity = _coercivity_certificate(expr, vars_)
    if not coercivity.certified:
        return PolynomialNegativityResult(
            vars_,
            None,
            None,
            False,
            "zeng",
            (
                "coercivity_not_certified",
                f"coercivity_method={coercivity.method}",
                *coercivity.notes,
            ),
        )
    coercivity_notes = (f"coercivity_method={coercivity.method}", *coercivity.notes)

    gradient = tuple(sp.expand(sp.diff(expr, variable)) for variable in vars_)
    gradient = tuple(g for g in gradient if g != 0)
    if not gradient:
        sign = certified_sign(expr.subs({v: 0 for v in vars_}))
        return PolynomialNegativityResult(vars_, sign == -1, None, sign is not None, "zeng")
    if not is_zero_dimensional(gradient, vars_):
        # For a coercive polynomial every global minimum is critical.  On the
        # critical locus, f < 0 is equivalent to the existence of a real slack
        # variable u satisfying f*u**2 + 1 = 0.  This turns negativity into a
        # pure polynomial-equality feasibility problem, exactly the fragment
        # handled by the ARS positive-dimensional backend.
        from .real_algebraic import real_algebraic_feasibility

        slack = sp.Dummy("zeng_slack", real=True)
        ars_vars = (*vars_, slack)
        negativity_equations = (*gradient, sp.expand(expr * slack**2 + 1))
        try:
            ars = real_algebraic_feasibility(negativity_equations, ars_vars, return_result=True)
        except (RationalUnivariateError, sp.PolynomialError, ValueError, TypeError):
            ars = None
        if ars is None or not ars.complete or ars.satisfiable is None:
            notes = (*coercivity_notes, "positive_dimensional_gradient_locus")
            if ars is not None:
                notes = (*notes, *ars.notes)
            return PolynomialNegativityResult(vars_, None, None, False, "zeng+ars", notes)
        if ars.satisfiable:
            assignment = ars.assignment or {}
            point_values = tuple(assignment[v] for v in vars_)
            verified = _certified_negative_point(expr, vars_, point_values)
            if verified is None:
                return PolynomialNegativityResult(
                    vars_,
                    None,
                    None,
                    False,
                    "zeng+ars",
                    (*coercivity_notes, "ars_witness_verification_failed"),
                )
            return PolynomialNegativityResult(
                vars_,
                True,
                verified,
                True,
                "zeng+ars",
                (*coercivity_notes, "positive_dimensional_gradient_locus"),
            )
        return PolynomialNegativityResult(
            vars_,
            False,
            None,
            True,
            "zeng+ars",
            (*coercivity_notes, "positive_dimensional_gradient_locus"),
        )
    try:
        points = solve_zero_dimensional_system(gradient, variables=vars_, real=True)
    except RationalUnivariateError:
        return PolynomialNegativityResult(
            vars_, None, None, False, "zeng", (*coercivity_notes, "critical_rur_failure")
        )

    undecidable = False
    for point in points:
        verified = _certified_negative_point(expr, vars_, point)
        if verified is not None:
            return PolynomialNegativityResult(
                vars_,
                True,
                verified,
                True,
                "zeng_critical_values",
                (*coercivity_notes, f"critical_points={len(points)}"),
            )
        assignment = dict(zip(vars_, point, strict=True))
        if certified_sign(expr.subs(assignment)) is None:
            undecidable = True
    if undecidable:
        return PolynomialNegativityResult(
            vars_,
            None,
            None,
            False,
            "zeng_critical_values",
            (*coercivity_notes, "critical_value_sign_undecidable"),
        )
    return PolynomialNegativityResult(
        vars_,
        False,
        None,
        True,
        "zeng_critical_values",
        (*coercivity_notes, f"critical_points={len(points)}"),
    )


def polynomial_nonnegative(
    polynomial: sp.Expr,
    variables: Sequence[sp.Symbol],
    *,
    strategy: str = "auto",
    sos_backend: str | Callable[[sp.Expr, list[sp.Symbol]], object] = "auto",
    return_result: bool = False,
    random_lines: int = 8,
    seed: int = 1234,
) -> bool | CertifiedDecisionResult:
    """Decide certified global polynomial nonnegativity.

    The default return is the mathematical Boolean. Set ``return_result=True``
    to obtain the structured certified portfolio trace, including the selected
    backend, ordered attempts, witness, and certificate.

    ``strategy`` may be ``"auto"``, ``"sos"``, ``"zeng"``, ``"ars"``,
    or ``"cad"``. Automatic mode tries exact-verified SOS search, then the
    Zeng/ARS specialized route, then complete CAD. ``random_lines`` and ``seed``
    tune only the Zeng witness-search stage.
    """
    from .decision_portfolio import _polynomial_nonnegative_portfolio

    return _polynomial_nonnegative_portfolio(
        polynomial,
        variables,
        strategy=strategy,
        sos_backend=sos_backend,
        return_result=return_result,
        random_lines=random_lines,
        seed=seed,
    )


def find_negative_point(
    polynomial: sp.Expr,
    variables: Sequence[sp.Symbol],
    **kwargs,
) -> dict[sp.Symbol, sp.Expr] | None:
    """Return a certified negative point or ``None`` for certified nonnegativity."""
    result = zeng_negative_point(polynomial, variables, **kwargs)
    negative = result.require_decision()
    return result.assignment if negative else None


__all__ = [
    "PolynomialNegativityResult",
    "find_negative_point",
    "polynomial_nonnegative",
    "zeng_negative_point",
]
