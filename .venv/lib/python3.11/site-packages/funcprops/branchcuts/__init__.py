"""Symbolic principal-branch cuts, boundary jumps, and restrictions.

The low-level API computes cut conditions and jumps.  The property-level
:func:`funcprops.branchcuts` query converts those loci into the common
funcprops singularity vocabulary.
"""

from ..engine import simplify_condition
from .computation import (
    branch_cut_jumps,
    branch_cut_jumps_pairs,
    condition_false_q,
    condition_simplify,
    differentiate_pairs,
    jump_pairs_from_table,
    numeric_branch_cut_jumps,
    unsatisfiable_q,
)
from .core import (
    IMPLEMENTED_RULE_NAMES,
    BranchRule,
    ConditionalJump,
    JumpTerm,
    RestrictedExpr,
    Undefined,
    apply_restrictions,
    branch_jump,
    get_branch_rule,
    propagate_undefined,
    restrict,
    simplify_restricted,
)
from .post_simplify import (
    merge_duplicate_pairs,
    post_simplify_pairs,
    post_simplify_restricted,
    rewrite_branch_expr,
    simplify_branch_condition,
)


def branch_cut_conditions(expr, variable, *, assumptions=True):
    """Return unique symbolic principal-cut conditions for ``expr``.

    This property query is intentionally cheaper than full jump propagation.
    It pulls primitive cut rules back through arguments and recursively includes
    inner cuts, but does not form plus/minus boundary-value products.
    """
    from sympy import Dummy, S, sympify

    expr = sympify(expr)
    variable = sympify(variable)
    out = []
    seen = set()

    def add(cond):
        cond = simplify_condition(sympify(cond), assumptions=assumptions)
        cond = simplify_branch_condition(cond, variable)
        if cond is S.false:
            return
        key = repr(cond)
        if key not in seen:
            seen.add(key)
            out.append(cond)

    def visit(e):
        e = sympify(e)
        if not e.has(variable):
            return
        # Exact primitive rule first.
        try:
            rule = get_branch_rule(e, variable)
        except (NotImplementedError, ValueError, TypeError, AttributeError):
            rule = None
        if rule is not None:
            for term in rule.jump_terms:
                add(term.cond)
        else:
            # Pull a one-variable primitive rule back through each varying argument.
            for i, arg in enumerate(getattr(e, "args", ())):
                if not getattr(arg, "has", lambda *_: False)(variable):
                    continue
                u = Dummy("funcprops_branch_arg")
                args = list(e.args)
                args[i] = u
                try:
                    template = e.func(*args)
                    rule = get_branch_rule(template, u)
                except (NotImplementedError, ValueError, TypeError, AttributeError):
                    continue
                for term in rule.jump_terms:
                    add(term.cond.xreplace({u: arg}))
        # Branches inside arguments remain branches of the composite expression.
        for arg in getattr(e, "args", ()):
            if getattr(arg, "has", lambda *_: False)(variable):
                visit(arg)

    visit(expr)
    return tuple(out)


def coverage_summary():
    """Return BranchCuts-rule coverage against the funcprops catalogue."""
    from ..registry import registry

    families = list(registry())
    return {
        "implemented_rules": {"total": len(IMPLEMENTED_RULE_NAMES)},
        "sympy": {"total": sum(bool(f.sympy_names) for f in families)},
        "mpmath": {"total": sum(bool(f.mpmath_names) for f in families)},
    }


__all__ = [
    "IMPLEMENTED_RULE_NAMES",
    "BranchRule",
    "ConditionalJump",
    "JumpTerm",
    "RestrictedExpr",
    "Undefined",
    "apply_restrictions",
    "branch_cut_conditions",
    "branch_cut_jumps",
    "branch_cut_jumps_pairs",
    "branch_jump",
    "condition_false_q",
    "condition_simplify",
    "coverage_summary",
    "differentiate_pairs",
    "get_branch_rule",
    "jump_pairs_from_table",
    "merge_duplicate_pairs",
    "numeric_branch_cut_jumps",
    "post_simplify_pairs",
    "post_simplify_restricted",
    "propagate_undefined",
    "restrict",
    "rewrite_branch_expr",
    "simplify_branch_condition",
    "simplify_restricted",
    "unsatisfiable_q",
]
