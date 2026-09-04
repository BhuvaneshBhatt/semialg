from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field

import sympy as sp

from ..algebraic.samples import sample_to_expr
from ..decomposition.cylindrical import CADResult, CellSet, cad
from ..exact_arithmetic import compare_exact_reals
from ..normalization import normalize_formula, normalize_problem_variables
from ..reconstruct.cylindrical import path_condition, section_value_bound
from ..reconstruct.radicals import fiber_root_candidates
from ..structural_keys import symbol_identity_key
from .bounds import (
    AlgebraicRootFunction,
    CADBound,
    DelineabilityCertificate,
    RootOrderCertificate,
    as_cad_bound,
    bound_expr,
)
from .lifting.stack import CADCell
from .polynomial_utils import exact_univariate_poly


@dataclass(frozen=True)
class CADBoundaryDescriptor:
    """Exact provenance for one finite delineable CAD-cell boundary.

    ``root_index`` is 1-based to match the public curve/surface evaluation
    APIs.  Explicit non-root boundaries intentionally return no descriptor.
    """

    polynomial: sp.Expr
    fiber_variable: sp.Symbol
    root_index: int
    level: int
    side: str
    base_variables: tuple[sp.Symbol, ...] = ()
    certificate: DelineabilityCertificate | None = None

    @property
    def pair(self) -> tuple[sp.Expr, int]:
        return (self.polynomial, self.root_index)

    @property
    def certified(self) -> bool:
        return self.certificate is not None and self.certificate.verify()


@dataclass(frozen=True)
class StructuredCADLevel:
    """One coordinate-level condition along a CAD cell path."""

    variable: sp.Symbol
    level: int
    kind: str
    lower: sp.Expr
    upper: sp.Expr
    sample: sp.Expr
    condition: sp.Expr
    index: tuple[int, ...]
    lower_bound: CADBound | None = None
    upper_bound: CADBound | None = None
    lower_closed: bool = False
    upper_closed: bool = False
    delineability: DelineabilityCertificate | None = None
    root_order: RootOrderCertificate | None = None

    def boundary_descriptor(self, side: str) -> CADBoundaryDescriptor | None:
        """Return the delineable algebraic boundary descriptor for ``side``.

        The descriptor is extracted from the typed CAD bound, so callers do
        not need to rediscover which projection polynomial/root generated a
        reconstructed radical or ``RootOf`` expression.
        """
        if side not in {"lower", "upper"}:
            raise ValueError("side must be 'lower' or 'upper'")
        bound = self.lower_bound if side == "lower" else self.upper_bound
        if not isinstance(bound, AlgebraicRootFunction):
            return None
        return CADBoundaryDescriptor(
            polynomial=sp.sympify(bound.polynomial),
            fiber_variable=bound.fiber_variable,
            root_index=bound.root_index + 1,
            level=self.level,
            side=side,
            base_variables=tuple(bound.base_variables),
            certificate=bound.certificate,
        )

    @property
    def boundary_descriptors(self) -> tuple[CADBoundaryDescriptor, ...]:
        out = (self.boundary_descriptor("lower"), self.boundary_descriptor("upper"))
        return tuple(item for item in out if item is not None)

    @property
    def is_section(self) -> bool:
        return self.kind == "section"

    @property
    def is_sector(self) -> bool:
        return self.kind == "sector"

    @property
    def dimension(self) -> int:
        return 0 if self.is_section else 1


@dataclass(frozen=True)
class StructuredCADCell:
    """Structured, user-facing representation of one CAD leaf cell."""

    variables: tuple[sp.Symbol, ...]
    index: tuple[int, ...]
    levels: tuple[StructuredCADLevel, ...]
    sample: Mapping[sp.Symbol, sp.Expr]
    source_cell: CADCell | None = None
    selected: bool = True
    _cell_levels: Mapping[int, Sequence[CADCell]] = field(
        default_factory=dict, repr=False, compare=False, hash=False
    )

    @property
    def dimension(self) -> int:
        return sum(level.dimension for level in self.levels)

    @property
    def bounded(self) -> bool:
        return all(level.lower != -sp.oo and level.upper != sp.oo for level in self.levels)

    @property
    def is_full_dimensional(self) -> bool:
        return self.dimension == len(self.variables)

    @property
    def boundary_descriptors(self) -> tuple[CADBoundaryDescriptor, ...]:
        """All finite algebraic root-function boundaries along this cell path."""
        return tuple(desc for level in self.levels for desc in level.boundary_descriptors)

    def boundary_pairs(
        self, *, variable: sp.Symbol | None = None
    ) -> tuple[tuple[sp.Expr, int], ...]:
        """Return unique ``(polynomial, 1-based root-index)`` boundary pairs."""
        seen: set[tuple[str, int]] = set()
        out: list[tuple[sp.Expr, int]] = []
        for desc in self.boundary_descriptors:
            if variable is not None and desc.fiber_variable != variable:
                continue
            key = (sp.srepr(sp.expand(desc.polynomial)), desc.root_index)
            if key not in seen:
                seen.add(key)
                out.append(desc.pair)
        return tuple(out)

    def as_formula(self, *, closed: bool = False) -> sp.Expr:
        if self.source_cell is not None:
            # The reconstructed path formula preserves algebraic section bounds
            # from the underlying CAD and is more faithful than recomposing from
            # the cached scalar fields when root functions are involved.
            return path_condition(
                self.source_cell,
                self.variables,
                self._cell_levels,
                closed=closed,
            )  # type: ignore[arg-type]
        parts = [level.condition for level in self.levels]
        return sp.And(*parts) if parts else sp.true


@dataclass(frozen=True)
class StructuredCADCellDecomposition:
    """Structured extraction of selected cells from a CAD computation."""

    variables: tuple[sp.Symbol, ...]
    cells: tuple[StructuredCADCell, ...]
    formula: sp.Expr
    cad_result: CADResult | None = None

    @property
    def full_dimensional_cells(self) -> tuple[StructuredCADCell, ...]:
        n = len(self.variables)
        return tuple(cell for cell in self.cells if cell.dimension == n)

    @property
    def sample_points(self) -> tuple[Mapping[sp.Symbol, sp.Expr], ...]:
        return tuple(cell.sample for cell in self.cells)


def _siblings(
    cells_by_level: Mapping[int, Sequence[CADCell]],
    level: int,
    parent_index: tuple[int, ...] | None,
) -> tuple[CADCell, ...]:
    return tuple(
        sorted(
            (cell for cell in cells_by_level.get(level, ()) if cell.parent_index == parent_index),
            key=lambda c: c.stack_position,
        )
    )


def _section_at_position(
    cells_by_level: Mapping[int, Sequence[CADCell]],
    level: int,
    parent_index: tuple[int, ...] | None,
    position: int,
) -> CADCell | None:
    for cell in _siblings(cells_by_level, level, parent_index):
        if cell.stack_position == position and cell.kind == "section":
            return cell
    return None


def _section_delineability_certificate(
    section: CADCell,
    variable: sp.Symbol,
    variables: Sequence[sp.Symbol],
    cells_by_level: Mapping[int, Sequence[CADCell]],
    *,
    sign_invariant: bool = True,
) -> DelineabilityCertificate:
    """Build a certificate that identifies and orders a CAD section over its base cell."""
    poly = sp.sympify(section.section_polynomial)
    base_subs = {
        variables[i]: sample_to_expr(section.sample[i])
        for i in range(min(section.level - 1, len(variables)))
    }
    sample_subs = dict(base_subs)
    sample_subs[variable] = sample_to_expr(section.sample[section.level - 1])
    sample_root_verified = bool(sp.simplify(poly.subs(sample_subs)) == 0)
    local_root_index = 0
    try:
        specialized = exact_univariate_poly(
            sp.expand(poly.subs(base_subs)), variable, algebraic_extension=False
        )
        roots = tuple(sp.real_roots(specialized.as_expr()))
        sample_value = sample_to_expr(section.sample[section.level - 1])
        matches = [i for i, root in enumerate(roots) if sp.simplify(root - sample_value) == 0]
        if matches:
            local_root_index = matches[0]
    except (sp.PolynomialError, ValueError, TypeError, NotImplementedError):
        local_root_index = 0
    siblings = _siblings(cells_by_level, section.level, section.parent_index)
    sections = [cell for cell in siblings if cell.kind == "section"]
    ordered = [cell.root_index for cell in sections]
    stack_order_verified = all(idx is not None for idx in ordered) and tuple(ordered) == tuple(
        sorted(ordered)
    )

    radical_branch_index = None
    representation_verified = False
    try:
        candidates = fiber_root_candidates(poly, variable, ordered=False)
        sample_value = sample_to_expr(section.sample[section.level - 1])
        matches = []
        for idx, candidate in enumerate(candidates):
            specialized_candidate = sp.simplify(candidate.subs(base_subs))
            if sp.simplify(specialized_candidate - sample_value) == 0:
                matches.append(idx)
        if len(matches) == 1:
            radical_branch_index = matches[0]
            # CAD delineability plus sign-invariance prevents a branch identity
            # from swapping without a projection/root event inside the base cell.
            representation_verified = bool(
                sign_invariant and stack_order_verified and sample_root_verified
            )
    except (sp.PolynomialError, ValueError, TypeError, NotImplementedError):
        radical_branch_index = None
        representation_verified = False

    regular_section_verified = False
    try:
        fiber_derivative = sp.diff(poly, variable).subs(sample_subs)
        sample_deriv_nonzero = bool(sp.simplify(fiber_derivative) != 0)
        # In a certified CAD lifting stack, delineability plus invariant root
        # order prevents a simple section from becoming multiple without a
        # projection event.  The nonzero sample derivative therefore upgrades
        # to a cell-wide regularity certificate only when the CAD certificate
        # itself is valid; it is not treated as a stand-alone proof.
        regular_section_verified = bool(
            sign_invariant
            and stack_order_verified
            and sample_root_verified
            and sample_deriv_nonzero
        )
    except (ValueError, TypeError, NotImplementedError):
        regular_section_verified = False

    return DelineabilityCertificate(
        polynomial=poly,
        fiber_variable=variable,
        root_index=int(local_root_index),
        base_variables=tuple(variables[: section.level - 1]),
        base_index=section.parent_index,
        section_index=section.index,
        defining_polynomial_key=section.defining_polynomial_key,
        stack_root_index=section.root_index,
        sign_invariant=sign_invariant,
        stack_order_verified=stack_order_verified,
        sample_root_verified=sample_root_verified,
        sample_root_value=sample_to_expr(section.sample[section.level - 1]),
        radical_branch_index=radical_branch_index,
        representation_verified=representation_verified,
        regular_section_verified=regular_section_verified,
        notes=("derived from CAD lifting stack order",),
    )


def _root_order_certificate(
    level_cell: CADCell,
    variable: sp.Symbol,
    cells_by_level: Mapping[int, Sequence[CADCell]],
) -> RootOrderCertificate | None:
    if level_cell.kind != "sector":
        return None
    left = _section_at_position(
        cells_by_level, level_cell.level, level_cell.parent_index, level_cell.stack_position - 1
    )
    right = _section_at_position(
        cells_by_level, level_cell.level, level_cell.parent_index, level_cell.stack_position + 1
    )
    lower_idx = None if left is None else left.root_index
    upper_idx = None if right is None else right.root_index
    adjacent = True
    if lower_idx is not None and upper_idx is not None:
        adjacent = upper_idx == lower_idx + 1
    order_verified = True
    if left is not None and right is not None:
        lv = sample_to_expr(left.sample[level_cell.level - 1])
        rv = sample_to_expr(right.sample[level_cell.level - 1])
        try:
            order_verified = compare_exact_reals(rv, lv) > 0
        except (TypeError, ValueError, NotImplementedError, sp.PolynomialError):
            order_verified = False
    return RootOrderCertificate(
        fiber_variable=variable,
        base_index=level_cell.parent_index,
        lower_root_index=lower_idx,
        upper_root_index=upper_idx,
        adjacent=adjacent,
        order_verified=order_verified,
        notes=("adjacent section roots in lifting stack",),
    )


def _level_typed_bounds(
    level_cell: CADCell,
    variable: sp.Symbol,
    variables: Sequence[sp.Symbol],
    cells_by_level: Mapping[int, Sequence[CADCell]],
    *,
    sign_invariant: bool = True,
) -> tuple[CADBound, CADBound, DelineabilityCertificate | None, RootOrderCertificate | None]:
    """Recover typed cylindrical bounds and their root certificates for one level.

    Sections become equal algebraic-root bounds. Sectors obtain their adjacent
    section roots when available, otherwise their stored finite or infinite
    bounds. Root-order metadata is attached only when the lifting stack proves it.
    """
    base_vars = tuple(variables[: level_cell.level - 1])
    if level_cell.kind == "section":
        cert = _section_delineability_certificate(
            level_cell, variable, variables, cells_by_level, sign_invariant=sign_invariant
        )
        bound = section_value_bound(
            level_cell, variable, base_variables=base_vars, certificate=cert, closed=True
        )
        return bound, bound, cert, None
    left_section = _section_at_position(
        cells_by_level, level_cell.level, level_cell.parent_index, level_cell.stack_position - 1
    )
    right_section = _section_at_position(
        cells_by_level, level_cell.level, level_cell.parent_index, level_cell.stack_position + 1
    )
    left_cert = None
    right_cert = None
    if left_section is not None:
        left_cert = _section_delineability_certificate(
            left_section, variable, variables, cells_by_level, sign_invariant=sign_invariant
        )
        lower = section_value_bound(
            left_section, variable, base_variables=base_vars, certificate=left_cert, closed=False
        )
    elif level_cell.lower_bound is not None:
        lower = as_cad_bound(level_cell.lower_bound, closed=False)
    else:
        lower = as_cad_bound(-sp.oo, closed=False)
    if right_section is not None:
        right_cert = _section_delineability_certificate(
            right_section, variable, variables, cells_by_level, sign_invariant=sign_invariant
        )
        upper = section_value_bound(
            right_section, variable, base_variables=base_vars, certificate=right_cert, closed=False
        )
    elif level_cell.upper_bound is not None:
        upper = as_cad_bound(level_cell.upper_bound, closed=False)
    else:
        upper = as_cad_bound(sp.oo, closed=False)
    cert = left_cert or right_cert
    return lower, upper, cert, _root_order_certificate(level_cell, variable, cells_by_level)


def _level_condition(
    variable: sp.Symbol, kind: str, lower: sp.Expr, upper: sp.Expr, *, closed: bool = False
) -> sp.Expr:
    if kind == "section" or sp.simplify(upper - lower) == 0:
        return sp.Eq(variable, lower)
    parts: list[sp.Expr] = []
    if lower != -sp.oo:
        parts.append(variable >= lower if closed else variable > lower)
    if upper != sp.oo:
        parts.append(variable <= upper if closed else variable < upper)
    return sp.And(*parts) if parts else sp.true


def _structured_cell_from_leaf(
    leaf: CADCell,
    variables: Sequence[sp.Symbol],
    cells_by_level: Mapping[int, Sequence[CADCell]],
    *,
    selected: bool = True,
    sign_invariant: bool = True,
) -> StructuredCADCell:
    """Reconstruct nested typed bounds for a full CAD leaf.

    Each prefix of the leaf index identifies the corresponding lifting cell.
    Reusing those prefixes preserves the exact algebraic branch selected during
    CAD construction instead of rediscovering roots numerically.
    """
    levels: list[StructuredCADLevel] = []
    sample_map: dict[sp.Symbol, sp.Expr] = {}
    for level, variable in enumerate(variables, start=1):
        prefix = leaf.index[:level]
        level_cell = next(cell for cell in cells_by_level[level] if cell.index == prefix)
        lower_bound, upper_bound, delineability, root_order = _level_typed_bounds(
            level_cell, variable, variables, cells_by_level, sign_invariant=sign_invariant
        )
        lower, upper = bound_expr(lower_bound), bound_expr(upper_bound)
        sample = sample_to_expr(level_cell.sample[level - 1])
        condition = _level_condition(variable, level_cell.kind, lower, upper)
        sample_map[variable] = sample
        levels.append(
            StructuredCADLevel(
                variable=variable,
                level=level,
                kind=level_cell.kind,
                lower=lower,
                upper=upper,
                sample=sample,
                condition=condition,
                index=level_cell.index,
                lower_bound=lower_bound,
                upper_bound=upper_bound,
                lower_closed=level_cell.kind == "section",
                upper_closed=level_cell.kind == "section",
                delineability=delineability,
                root_order=root_order,
            )
        )
    return StructuredCADCell(
        variables=tuple(variables),
        index=leaf.index,
        levels=tuple(levels),
        sample=sample_map,
        source_cell=leaf,
        selected=selected,
        _cell_levels=cells_by_level,
    )


def extract_structured_cad_cells(
    condition_or_cad: object,
    variables: Sequence[sp.Symbol | str] | None = None,
    *,
    selected_only: bool = True,
) -> StructuredCADCellDecomposition:
    """Extract structured cells from a CAD result or from a formula.

    The representation is built directly from CAD lifting data and therefore
    supports algebraic stack bounds such as ``-sqrt(x) < y < sqrt(x)`` from a
    cell over ``y**2 - x``.
    """

    if isinstance(condition_or_cad, CADResult):
        result = condition_or_cad
    elif isinstance(condition_or_cad, CellSet):
        vars_ = tuple(condition_or_cad.variables)
        cells = tuple(
            _structured_cell_from_leaf(cell, vars_, condition_or_cad.cells_by_level, selected=True)
            for cell in condition_or_cad.cells
        )
        return StructuredCADCellDecomposition(vars_, cells, condition_or_cad.formula)
    else:
        formula = normalize_formula(condition_or_cad)
        if variables is None:
            variables = tuple(sorted(formula.free_symbols, key=symbol_identity_key))
        result = cad(
            formula,
            normalize_problem_variables(variables, formula),
            output="cells",
            return_result=True,
        )

    vars_ = tuple(result.variables)
    if selected_only:
        leaves = tuple(result.cell_set.cells)
    else:
        leaves = tuple(result.cad.cells_by_level.get(len(vars_), ()))
    selected_indices = {cell.index for cell in result.cell_set.cells}
    diagnostics = result.cad.diagnostics
    sign_invariant = (
        not diagnostics.invariant_failures
        if diagnostics is not None
        else result.cad.verify_sign_invariance().ok
    )
    structured = tuple(
        _structured_cell_from_leaf(
            cell,
            vars_,
            result.cad.cells_by_level,
            selected=cell.index in selected_indices,
            sign_invariant=sign_invariant,
        )
        for cell in leaves
    )
    return StructuredCADCellDecomposition(vars_, structured, result.formula, result)


__all__ = [
    "CADBoundaryDescriptor",
    "StructuredCADLevel",
    "StructuredCADCell",
    "StructuredCADCellDecomposition",
    "extract_structured_cad_cells",
]
