"""Checked equality preprocessing with explicit exceptional parameter branches."""

import sympy as sp
from algroots import groebner_basis
from algroots.parametric_groebner import parametric_groebner_basis


def equality_domain(formula):
    """Retain visible poles of outer equalities before coefficient-field cancellation."""
    atoms = sp.And.make_args(formula)
    # A coefficient field cancels removable factors. Real predicates must
    # retain poles that are still present in the supplied expression tree.
    presentation_domain = sp.And(
        *(
            sp.Ne(sp.together(power.base).as_numer_denom()[0], 0)
            for equation in atoms
            if isinstance(equation, sp.Equality)
            for power in equation.atoms(sp.Pow)
            if power.exp.is_Integer and power.exp < 0
        )
    )
    return presentation_domain


def equality_numerators(formula):
    """Clear outer equality denominators; pair the result with equality_domain."""
    return sp.And(
        *(
            sp.Eq(sp.together(atom.lhs - atom.rhs).as_numer_denom()[0], 0)
            if isinstance(atom, sp.Equality)
            else atom
            for atom in sp.And.make_args(formula)
        )
    )


def guarded_equality_preprocessing(formula, variables, parameters, *, backend, limits=None):
    """Replace conjunctive equalities only where exact specialization identities apply.

    The exceptional branch retains every original equality and inequality.
    The caller must eliminate the complete returned formula, not only its
    generic branch. General disjunctions without common outer equalities are
    returned unchanged for the complete real-elimination backend.
    """
    atoms = sp.And.make_args(formula)
    equalities = tuple(
        sp.together(atom.lhs - atom.rhs).as_numer_denom()[0]
        for atom in atoms
        if isinstance(atom, sp.Equality)
    )
    if not equalities or not variables:
        return formula, {"groebner_preprocessing": False}
    remaining = tuple(atom for atom in atoms if not isinstance(atom, sp.Equality))
    presentation_domain = equality_domain(formula)
    remaining += (presentation_domain,)
    if parameters:
        proof = parametric_groebner_basis(
            equalities, variables, parameters, backend=backend, limits=limits
        )
        ground = proof.domain.domain
        if ground.is_AlgebraicField and ground.ext.as_expr().is_real is not True:
            raise ValueError("real elimination requires a real algebraic coefficient embedding")
        return sp.And(proof.as_formula(), *remaining), {
            "groebner_preprocessing": True,
            "groebner_backend": backend,
            "groebner_guard": proof.guard,
            "groebner_input_guard": proof.input_guard,
            "groebner_presentation_domain": presentation_domain,
            "generic_basis": tuple(proof.basis),
            "exceptional_condition": proof.exceptional_condition,
            "exceptional_branch_preserved": True,
        }
    basis = groebner_basis(equalities, variables, backend=backend, limits=limits)
    return sp.And(*(sp.Eq(p, 0) for p in basis), *remaining), {
        "groebner_preprocessing": True,
        "groebner_backend": backend,
        "generic_basis": tuple(basis),
    }
