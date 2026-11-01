"""Exact real coverage of birational curve normalization maps."""

from dataclasses import dataclass

import sympy as sp
from algroots.curves import CurveParametrization

from .decision import is_satisfiable
from .decision.solution import EquivalenceResult
from .function_verification import _symbols
from .qe import QuantifierEliminationResult, quantifier_eliminate
from .simplify.atoms import normalize_atoms


@dataclass(frozen=True)
class RealCurveParametrizationVerification:
    """Real image of P1(R), including its finite image at the parameter infinity.

    A birational normalization need not cover isolated real singularities: their
    geometric branches can be complex conjugate. Coverage is an independent
    exact decision, rather than a consequence of generic map identities.
    """

    parametrization: CurveParametrization
    coordinates: tuple[sp.Symbol, sp.Symbol]
    domain: sp.Expr
    finite_parameter_graph: sp.Expr
    infinity_point: tuple[sp.Expr, sp.Expr] | None
    image_formula: sp.Expr
    curve_formula: sp.Expr
    elimination: QuantifierEliminationResult
    comparison: EquivalenceResult

    @property
    def covers_curve(self):
        return self.comparison.equivalent

    @property
    def certified(self):
        return self.elimination.certified

    @property
    def missing_condition(self):
        if self.covers_curve:
            return sp.false
        return normalize_atoms(sp.And(self.domain, self.curve_formula, sp.Not(self.image_formula)))


def verify_real_curve_parametrization(parametrization, *, coordinates=None, domain=sp.true):
    """Compare the complete real normalization image with a real plane curve.

    Both finite and infinite real parameter charts are included. The returned
    formulas describe the actual real image; ``covers_curve`` may be false even
    for a valid birational map. ``domain`` may restrict real branches or exclude
    singularities. Exact solver support/resource failures propagate.
    """
    if not isinstance(parametrization, CurveParametrization) or not parametrization.verify():
        raise ValueError("a verified algroots CurveParametrization is required")
    constants = parametrization.domain
    if constants != sp.QQ and constants.ext.as_expr().is_real is not True:
        raise ValueError("real curve coverage requires a real embedding of the constants")
    coordinates = parametrization.curve.variables if coordinates is None else tuple(coordinates)
    if len(coordinates) != 2:
        raise ValueError("two output coordinates are required")
    _symbols((*coordinates, parametrization.parameter), "curve and parameter coordinates")
    domain = sp.sympify(domain)
    if domain.free_symbols - set(coordinates):
        raise ValueError("domain may use only the output coordinates")
    x, y = coordinates
    a, b, c = parametrization.projective_coordinates
    graph = sp.And(sp.Eq(x * c - a, 0), sp.Eq(y * c - b, 0), sp.Ne(c, 0))
    inverse_denominator = sp.fraction(sp.cancel(parametrization.inverse))[1].xreplace(
        dict(zip(parametrization.curve.variables, coordinates, strict=True))
    )
    # Birational replay proves coverage wherever the real inverse is defined.
    # Only its finite exceptional locus needs QE; this preserves isolated real
    # singularities whose normalization preimages are nonreal.
    exceptional_graph = sp.And(graph, sp.Eq(inverse_denominator, 0))
    projected = quantifier_eliminate(
        exceptional_graph,
        [("exists", parametrization.parameter)],
        variables=(*coordinates, parametrization.parameter),
        return_result=True,
    )
    at_infinity = parametrization.point_at_infinity
    point = None
    extra = sp.false
    if at_infinity[2] != 0:
        point = tuple(sp.cancel(v / at_infinity[2]) for v in at_infinity[:2])
        extra = sp.And(sp.Eq(x, point[0]), sp.Eq(y, point[1]))
    curve = normalize_atoms(
        sp.Eq(
            parametrization.curve.polynomial.as_expr().xreplace(
                dict(zip(parametrization.curve.variables, coordinates, strict=True))
            ),
            0,
        )
    )
    image = normalize_atoms(
        sp.Or(sp.And(curve, sp.Ne(inverse_denominator, 0)), projected.formula, extra)
    )
    # The verified maps prove image => curve and coverage on D != 0.
    # The remaining decision is finite; expose D=0 explicitly so exact RUR
    # dispatch can handle algebraic constants without an ambient CAD.
    missing = normalize_atoms(
        sp.And(
            domain,
            curve,
            sp.Eq(inverse_denominator, 0),
            sp.Not(sp.Or(projected.formula, extra)),
        )
    )
    decision = is_satisfiable(missing, coordinates, return_result=True)
    comparison = EquivalenceResult(
        not bool(decision),
        normalize_atoms(sp.And(domain, image)),
        normalize_atoms(sp.And(domain, curve)),
        coordinates,
        counterexample=decision.witness,
        failed_direction="rhs_implies_lhs" if decision else None,
        method="birational_inverse_chart_and_exceptional_projection",
        diagnostics={"exceptional_locus_decision": decision},
    )
    return RealCurveParametrizationVerification(
        parametrization, coordinates, domain, graph, point, image, curve, projected, comparison
    )


__all__ = ["RealCurveParametrizationVerification", "verify_real_curve_parametrization"]
