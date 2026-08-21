from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
from functools import cmp_to_key

import sympy as sp

from ..algebraic.cache import CACHE
from ..algebraic.roots import isolate_real_roots
from ..algebraic.samples import AlgebraicRoot, RationalSample, Sample, sample_to_expr
from ..errors import ExactEvaluationFailure
from ..exact_arithmetic import compare_exact_reals


class CADBound:
    """Typed endpoint of a cylindrical CAD coordinate constraint."""

    closed: bool

    def as_expr(self) -> sp.Expr:
        raise NotImplementedError

    def evaluate(self, base_point: Mapping[sp.Symbol, object] | None = None) -> object:
        return self.as_expr()

    @property
    def is_infinite(self) -> bool:
        return False

    def __sympy__(self):
        return self.as_expr()


@dataclass(frozen=True)
class ExplicitCADBound(CADBound):
    value: sp.Expr
    closed: bool = False

    def as_expr(self) -> sp.Expr:
        return sp.sympify(self.value)


@dataclass(frozen=True)
class InfiniteCADBound(CADBound):
    sign: int
    closed: bool = False

    def __post_init__(self) -> None:
        if self.sign not in (-1, 1):
            raise ValueError("infinite bound sign must be -1 or 1")

    def as_expr(self) -> sp.Expr:
        return -sp.oo if self.sign < 0 else sp.oo

    @property
    def is_infinite(self) -> bool:
        return True


@dataclass(frozen=True)
class AlgebraicNumberBound(CADBound):
    sample: Sample
    closed: bool = False

    def as_expr(self) -> sp.Expr:
        return sample_to_expr(self.sample)

    def evaluate(self, base_point: Mapping[sp.Symbol, object] | None = None) -> Sample:
        return self.sample


@dataclass(frozen=True)
class DelineabilityCertificate:
    polynomial: sp.Expr
    fiber_variable: sp.Symbol
    root_index: int
    base_variables: tuple[sp.Symbol, ...] = ()
    base_index: tuple[int, ...] | None = None
    section_index: tuple[int, ...] | None = None
    defining_polynomial_key: str | None = None
    stack_root_index: int | None = None
    sign_invariant: bool = False
    stack_order_verified: bool = False
    sample_root_verified: bool = False
    sample_root_value: sp.Expr | None = None
    radical_branch_index: int | None = None
    representation_verified: bool = False
    regular_section_verified: bool = False
    notes: tuple[str, ...] = ()

    @property
    def certified(self) -> bool:
        return self.sign_invariant and self.stack_order_verified and self.sample_root_verified

    def verify(self) -> bool:
        return (
            self.certified
            and self.root_index >= 0
            and self.fiber_variable in sp.sympify(self.polynomial).free_symbols
        )

    @property
    def regular(self) -> bool:
        return self.verify() and self.regular_section_verified

    def verify_regularity(self) -> bool:
        return self.regular


@dataclass(frozen=True)
class RootOrderCertificate:
    fiber_variable: sp.Symbol
    base_index: tuple[int, ...] | None
    lower_root_index: int | None
    upper_root_index: int | None
    adjacent: bool
    order_verified: bool
    notes: tuple[str, ...] = ()

    @property
    def certified(self) -> bool:
        return self.adjacent and self.order_verified

    def verify(self) -> bool:
        return self.certified


@dataclass(frozen=True)
class CertifiedRootComparison:
    """Result of a conservative certified comparison of algebraic roots.

    ``relation`` is ``-1``, ``0`` or ``1`` when the left operand is proven
    smaller than, equal to, or larger than the right operand.  ``None`` means
    that semialg deliberately declined to infer a global ordering.
    """

    relation: int | None
    certified: bool
    scope: str = "unknown"
    reason: str = ""

    def verify(self) -> bool:
        return self.certified and self.relation in (-1, 0, 1)


def _exact_expr_comparison(left: sp.Expr, right: sp.Expr) -> int | None:
    diff = sp.simplify(sp.sympify(left) - sp.sympify(right))
    if diff == 0 or diff.is_zero is True:
        return 0
    if diff.is_negative is True:
        return -1
    if diff.is_positive is True:
        return 1
    try:
        sign = sp.sign(diff)
        if sign in (-1, 0, 1):
            return int(sign)
    except (TypeError, ValueError, NotImplementedError):
        pass
    return None


@dataclass(frozen=True)
class AlgebraicRootFunction(CADBound):
    """A delineable real root function over a CAD base cell.

    ``root_index`` is zero-based in the sorted real-root order used by the CAD
    lifting stack.  The polynomial coefficients may depend on ``base_variables``.
    """

    polynomial: sp.Expr
    fiber_variable: sp.Symbol
    root_index: int
    base_variables: tuple[sp.Symbol, ...] = ()
    base_index: tuple[int, ...] | None = None
    certificate: DelineabilityCertificate | None = None
    stack_root_index: int | None = None
    closed: bool = False

    def __post_init__(self) -> None:
        if self.root_index < 0:
            raise ValueError("root_index must be nonnegative")
        if self.fiber_variable not in sp.sympify(self.polynomial).free_symbols:
            raise ValueError("fiber variable must occur in the defining polynomial")

    def as_expr(self) -> sp.Expr:
        from ..reconstruct.radicals import fiber_root_candidates
        from ..reconstruct.root_functions import root_of

        poly = sp.expand(self.polynomial)
        # For coefficient-dependent roots, the CAD root index is an ordered-root
        # identity over a base cell.  A quadratic formula branch is only a safe
        # presentation when lifting certified which branch has that identity.
        cert = self.certificate
        if (
            cert is not None
            and cert.representation_verified
            and cert.radical_branch_index is not None
        ):
            candidates = fiber_root_candidates(poly, self.fiber_variable, ordered=False)
            if 0 <= cert.radical_branch_index < len(candidates):
                return candidates[cert.radical_branch_index]
        if poly.free_symbols <= {self.fiber_variable}:
            try:
                roots = tuple(sp.real_roots(poly))
                if 0 <= self.root_index < len(roots):
                    return roots[self.root_index]
            except (sp.PolynomialError, ValueError, TypeError, NotImplementedError):
                pass
        return root_of(poly, self.fiber_variable, sp.Integer(self.root_index))

    def evaluate(self, base_point: Mapping[sp.Symbol, object] | None = None) -> Sample:
        substitutions: dict[sp.Symbol, sp.Expr] = {}
        for var in self.base_variables:
            if base_point is None or var not in base_point:
                raise ValueError(f"missing base coordinate {var}")
            value = base_point[var]
            if isinstance(value, (RationalSample, AlgebraicRoot)):
                substitutions[var] = sample_to_expr(value)
            else:
                substitutions[var] = sp.sympify(value)
        eval_key = (
            "root-evaluate",
            sp.srepr(sp.expand(self.polynomial)),
            sp.srepr(self.fiber_variable),
            int(self.root_index),
            tuple((sp.srepr(v), sp.srepr(substitutions[v])) for v in self.base_variables),
        )
        cached_eval = CACHE.specializations.get(eval_key)
        if cached_eval is not None:
            CACHE.stats.specialization_hits += 1
            return cached_eval  # type: ignore[return-value]
        CACHE.stats.specialization_misses += 1
        specialized = sp.expand(sp.sympify(self.polynomial).subs(substitutions))
        if specialized.free_symbols - {self.fiber_variable}:
            raise ValueError("base point did not specialize all coefficient variables")
        try:
            roots = isolate_real_roots(sp.Poly(specialized, self.fiber_variable, domain="EX"))
            if 0 <= self.root_index < len(roots):
                result = roots[self.root_index]
                CACHE.specializations.put(eval_key, result)
                return result
        except (sp.PolynomialError, ValueError, TypeError, NotImplementedError):
            pass

        # RUR is a second exact backend for fixed-coefficient zero-dimensional
        # fiber problems.  It returns exact expressions, which we then isolate.
        try:
            from ..algebraic.rational_univariate import solve_zero_dimensional_system_with_rur

            points = solve_zero_dimensional_system_with_rur(
                (specialized,), (self.fiber_variable,), real=True
            )
            values = sorted((p[0] for p in points), key=cmp_to_key(compare_exact_reals))
            if 0 <= self.root_index < len(values):
                value = values[self.root_index]
                roots = isolate_real_roots(
                    sp.Poly(sp.minpoly(value, self.fiber_variable), self.fiber_variable)
                )
                for root in roots:
                    if sp.simplify(root.as_expr() - value) == 0:
                        CACHE.specializations.put(eval_key, root)
                        return root
        except (
            sp.PolynomialError,
            ValueError,
            TypeError,
            NotImplementedError,
            RuntimeError,
        ) as exc:
            raise ExactEvaluationFailure(
                "could not evaluate algebraic root function exactly"
            ) from exc
        raise IndexError("root index is not present in the specialized real fiber")

    @property
    def is_regular(self) -> bool:
        """Whether regularity is certified over the entire CAD base cell."""

        return self.certificate is not None and self.certificate.verify_regularity()

    def fiber_derivative_expr(self) -> sp.Expr:
        """Return ``d polynomial / d fiber`` evaluated on this root function."""

        root = self.as_expr()
        return sp.simplify(
            sp.diff(sp.sympify(self.polynomial), self.fiber_variable).subs(
                self.fiber_variable, root
            )
        )

    def derivative_expr(self, variable: sp.Symbol, *, require_regular: bool = True) -> sp.Expr:
        """Return the certified implicit derivative with respect to a base variable.

        A nonzero fiber derivative over the base cell is required by default;
        callers may request the formal quotient with ``require_regular=False``.
        """

        variable = sp.sympify(variable)
        if variable not in self.base_variables:
            return sp.Integer(0)
        if require_regular and not self.is_regular:
            raise ValueError(
                "implicit derivative requires certified regularity on the CAD base cell"
            )
        root = self.as_expr()
        p = sp.sympify(self.polynomial)
        numerator = sp.diff(p, variable).subs(self.fiber_variable, root)
        denominator = sp.diff(p, self.fiber_variable).subs(self.fiber_variable, root)
        if require_regular and sp.simplify(denominator) == 0:
            raise ValueError("fiber derivative vanishes on the algebraic section")
        return sp.cancel(-numerator / denominator)

    def specialize(self, base_point: Mapping[sp.Symbol, object]) -> CADBound:
        """Specialize coefficient variables exactly.

        Full specialization returns an :class:`AlgebraicNumberBound`.  Partial
        specialization returns a new root function and intentionally drops the
        old cell certificate, because its domain certificate belongs to the
        original CAD base cell.
        """

        substitutions: dict[sp.Symbol, sp.Expr] = {}
        for var, value in base_point.items():
            if var not in self.base_variables:
                continue
            substitutions[var] = (
                sample_to_expr(value)
                if isinstance(value, (RationalSample, AlgebraicRoot))
                else sp.sympify(value)
            )
        spec_key = (
            "root-specialize",
            sp.srepr(sp.expand(self.polynomial)),
            sp.srepr(self.fiber_variable),
            int(self.root_index),
            tuple(sp.srepr(v) for v in self.base_variables),
            tuple(
                (sp.srepr(v), sp.srepr(substitutions[v]))
                for v in self.base_variables
                if v in substitutions
            ),
            bool(self.closed),
        )
        cached_spec = CACHE.specializations.get(spec_key)
        if cached_spec is not None:
            CACHE.stats.specialization_hits += 1
            return cached_spec  # type: ignore[return-value]
        CACHE.stats.specialization_misses += 1
        remaining = tuple(var for var in self.base_variables if var not in substitutions)
        if not remaining:
            result: CADBound = AlgebraicNumberBound(self.evaluate(base_point), closed=self.closed)
        else:
            specialized = sp.expand(sp.sympify(self.polynomial).subs(substitutions))
            result = AlgebraicRootFunction(
                specialized,
                self.fiber_variable,
                self.root_index,
                remaining,
                self.base_index,
                None,
                self.stack_root_index,
                self.closed,
            )
        CACHE.specializations.put(spec_key, result)
        return result

    def compare_certified(
        self,
        other: object,
        *,
        base_point: Mapping[sp.Symbol, object] | None = None,
    ) -> CertifiedRootComparison:
        """Conservatively compare this root with another bound or exact value.

        Global root-order comparisons are certified only when both root
        functions belong to the same delineable polynomial stack/base cell.
        With ``base_point`` supplied, an exact pointwise comparison is attempted.
        """

        if isinstance(other, AlgebraicRootFunction):
            same_stack = (
                sp.expand(self.polynomial - other.polynomial) == 0
                and self.fiber_variable == other.fiber_variable
                and self.base_variables == other.base_variables
                and self.base_index == other.base_index
            )
            certs_ok = (
                self.certificate is not None
                and other.certificate is not None
                and self.certificate.verify()
                and other.certificate.verify()
            )
            if same_stack and certs_ok:
                relation = (
                    -1
                    if self.root_index < other.root_index
                    else (1 if self.root_index > other.root_index else 0)
                )
                return CertifiedRootComparison(
                    relation, True, "base-cell", "ordered roots of one certified CAD stack"
                )
        if base_point is not None:
            try:
                left = sample_to_expr(self.evaluate(base_point))
                if isinstance(other, AlgebraicRootFunction):
                    right = sample_to_expr(other.evaluate(base_point))
                elif isinstance(other, AlgebraicNumberBound):
                    right = sample_to_expr(other.sample)
                elif isinstance(other, CADBound):
                    right = sp.sympify(other.evaluate(base_point))
                else:
                    right = sp.sympify(other)
                relation = _exact_expr_comparison(left, right)
                return CertifiedRootComparison(
                    relation, relation is not None, "point", "exact specialization comparison"
                )
            except (
                ValueError,
                IndexError,
                TypeError,
                sp.PolynomialError,
                NotImplementedError,
                ExactEvaluationFailure,
            ):
                pass
        return CertifiedRootComparison(None, False, "unknown", "global ordering is not certified")


def as_cad_bound(value: object, *, closed: bool = False) -> CADBound:
    if isinstance(value, CADBound):
        if value.closed == closed:
            return value
        if isinstance(value, ExplicitCADBound):
            return ExplicitCADBound(value.value, closed)
        if isinstance(value, InfiniteCADBound):
            return InfiniteCADBound(value.sign, closed)
        if isinstance(value, AlgebraicNumberBound):
            return AlgebraicNumberBound(value.sample, closed)
        if isinstance(value, AlgebraicRootFunction):
            return AlgebraicRootFunction(
                value.polynomial,
                value.fiber_variable,
                value.root_index,
                value.base_variables,
                value.base_index,
                value.certificate,
                value.stack_root_index,
                closed,
            )
    if value == -sp.oo:
        return InfiniteCADBound(-1, closed)
    if value == sp.oo:
        return InfiniteCADBound(1, closed)
    if isinstance(value, (RationalSample, AlgebraicRoot)):
        return AlgebraicNumberBound(value, closed)
    return ExplicitCADBound(sp.sympify(value), closed)


def bound_expr(bound: CADBound | sp.Expr) -> sp.Expr:
    return bound.as_expr() if isinstance(bound, CADBound) else sp.sympify(bound)


__all__ = [
    "CADBound",
    "ExplicitCADBound",
    "InfiniteCADBound",
    "AlgebraicNumberBound",
    "AlgebraicRootFunction",
    "CertifiedRootComparison",
    "DelineabilityCertificate",
    "RootOrderCertificate",
    "as_cad_bound",
    "bound_expr",
]


@dataclass(frozen=True)
class CADBoundLevelVerification:
    level: int
    variable: sp.Symbol
    dependencies_valid: bool
    section_verified: bool
    root_order_verified: bool
    repr_consistent: bool
    sample_contained: bool
    openness_valid: bool
    notes: tuple[str, ...] = ()

    def verify(self) -> bool:
        return (
            self.dependencies_valid
            and self.section_verified
            and self.root_order_verified
            and self.repr_consistent
            and self.sample_contained
            and self.openness_valid
        )


@dataclass(frozen=True)
class CADCellBoundsCertificate:
    cell_index: tuple[int, ...]
    levels: tuple[CADBoundLevelVerification, ...]
    source_path_consistent: bool
    certified: bool
    notes: tuple[str, ...] = ()

    def verify(self) -> bool:
        return (
            self.certified
            and self.source_path_consistent
            and all(level.verify() for level in self.levels)
        )


def _eval_bound_expr(bound: CADBound, previous: Mapping[sp.Symbol, object]) -> sp.Expr:
    if isinstance(bound, AlgebraicRootFunction):
        cert = bound.certificate
        if cert is not None and cert.sample_root_verified and cert.sample_root_value is not None:
            return sp.sympify(cert.sample_root_value)
        return sample_to_expr(bound.evaluate(previous))
    return sp.simplify(
        bound.as_expr().subs(
            {
                k: sample_to_expr(v)
                if isinstance(v, (RationalSample, AlgebraicRoot))
                else sp.sympify(v)
                for k, v in previous.items()
            }
        )
    )


def verify_cad_cell_bounds(cell: object) -> CADCellBoundsCertificate:
    """Verify typed cylindrical bounds for a reconstructed CAD cell.

    The verifier is intentionally local and auditable: it checks triangular
    dependencies, section/root certificates, adjacent-root ordering, sample
    containment, and open/closed semantics at every coordinate level.
    """

    variables = tuple(getattr(cell, "variables", ()))
    levels = tuple(getattr(cell, "levels", ()))
    sample_map = dict(getattr(cell, "sample", {}))
    checks: list[CADBoundLevelVerification] = []
    previous: dict[sp.Symbol, object] = {}
    for i, level in enumerate(levels):
        var = level.variable
        allowed = set(variables[:i])
        lower = (
            level.typed_lower
            if hasattr(level, "typed_lower")
            else as_cad_bound(level.lower, closed=getattr(level, "lower_closed", False))
        )
        upper = (
            level.typed_upper
            if hasattr(level, "typed_upper")
            else as_cad_bound(level.upper, closed=getattr(level, "upper_closed", False))
        )
        lower_expr = lower.as_expr()
        upper_expr = upper.as_expr()
        dependencies_valid = (lower_expr.free_symbols - allowed - {var} == set()) and (
            upper_expr.free_symbols - allowed - {var} == set()
        )
        # Algebraic root-function placeholders contain the fiber variable in
        # their display expression; the underlying coefficients must depend
        # only on preceding variables.
        for bound in (lower, upper):
            if isinstance(bound, AlgebraicRootFunction):
                coefficient_symbols = sp.sympify(bound.polynomial).free_symbols - {
                    bound.fiber_variable
                }
                dependencies_valid = dependencies_valid and coefficient_symbols <= allowed
                if bound.certificate is not None:
                    dependencies_valid = dependencies_valid and bound.certificate.verify()
        section_verified = True
        if getattr(level, "is_section", False):
            section_verified = sp.simplify(lower_expr - upper_expr) == 0
            if getattr(level, "delineability", None) is not None:
                section_verified = section_verified and level.delineability.verify()
        root_order_verified = True
        if getattr(level, "root_order", None) is not None:
            root_order_verified = level.root_order.verify()
        repr_consistent = True
        for bound in (lower, upper):
            if not isinstance(bound, AlgebraicRootFunction):
                continue
            cert = bound.certificate
            if cert is None:
                repr_consistent = False
                continue
            display = bound.as_expr()
            # root_of is the typed root identity itself.  For a readable radical,
            # additionally check that its specialization equals the certified root.
            if getattr(display.func, "__name__", "") != "root_of":
                try:
                    substitutions = {
                        k: sample_to_expr(v)
                        if isinstance(v, (RationalSample, AlgebraicRoot))
                        else sp.sympify(v)
                        for k, v in previous.items()
                    }
                    displayed_value = sp.simplify(display.subs(substitutions))
                    if cert.sample_root_value is None:
                        repr_consistent = False
                    else:
                        repr_consistent = (
                            repr_consistent
                            and sp.simplify(displayed_value - cert.sample_root_value) == 0
                        )
                except (ValueError, IndexError, sp.PolynomialError, TypeError, NotImplementedError):
                    repr_consistent = False
            else:
                repr_consistent = repr_consistent and cert.verify()
        sample_contained = True
        try:
            value = sp.sympify(sample_map[var])
            lo = _eval_bound_expr(lower, previous)
            hi = _eval_bound_expr(upper, previous)
            if getattr(level, "is_section", False):
                sample_contained = sp.simplify(value - lo) == 0
            else:
                if lo == -sp.oo:
                    left_ok = True
                else:
                    left_cmp = compare_exact_reals(value, lo)
                    left_ok = left_cmp > 0 or (
                        left_cmp == 0 and getattr(level, "lower_closed", False)
                    )
                if hi == sp.oo:
                    right_ok = True
                else:
                    right_cmp = compare_exact_reals(value, hi)
                    right_ok = right_cmp < 0 or (
                        right_cmp == 0 and getattr(level, "upper_closed", False)
                    )
                sample_contained = left_ok and right_ok
        except (ValueError, IndexError, TypeError, sp.PolynomialError, NotImplementedError):
            sample_contained = False
        openness_valid = True
        if getattr(level, "is_section", False):
            openness_valid = getattr(level, "lower_closed", True) and getattr(
                level, "upper_closed", True
            )
        checks.append(
            CADBoundLevelVerification(
                level=i + 1,
                variable=var,
                dependencies_valid=dependencies_valid,
                section_verified=section_verified,
                root_order_verified=root_order_verified,
                repr_consistent=repr_consistent,
                sample_contained=sample_contained,
                openness_valid=openness_valid,
            )
        )
        previous[var] = sample_map.get(var)
    source = getattr(cell, "source_cell", None)
    source_path_consistent = source is None or getattr(
        source, "index", getattr(cell, "index", None)
    ) == getattr(cell, "index", None)
    certified = source_path_consistent and all(check.verify() for check in checks)
    return CADCellBoundsCertificate(
        tuple(getattr(cell, "index", ())), tuple(checks), source_path_consistent, certified
    )


__all__.extend(
    [
        "CADBoundLevelVerification",
        "CADCellBoundsCertificate",
        "verify_cad_cell_bounds",
    ]
)
