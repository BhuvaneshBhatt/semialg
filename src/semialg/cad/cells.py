from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass

import sympy as sp

from ..algebraic.samples import sample_to_expr
from ..decomposition.cylindrical import CADResult, CellSet, cad
from ..exact_arithmetic import compare_exact_reals
from ..implicit_geometry import VerticalBoundCell2D, _normalize_formula, _normalize_variables
from ..reconstruct.cylindrical import path_condition, section_value_bound
from ..reconstruct.radicals import fiber_root_candidates
from .bounds import (
    AlgebraicRootFunction,
    CADBound,
    CertifiedRootComparison,
    DelineabilityCertificate,
    RootOrderCertificate,
    as_cad_bound,
    bound_expr,
)
from .lifting.stack import CADCell


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

    @property
    def dimension(self) -> int:
        return sum(level.dimension for level in self.levels)

    @property
    def bounded(self) -> bool:
        return all(level.lower != -sp.oo and level.upper != sp.oo for level in self.levels)

    @property
    def is_full_dimensional(self) -> bool:
        return self.dimension == len(self.variables)

    def as_formula(self, *, closed: bool = False) -> sp.Expr:
        if self.source_cell is not None:
            # The reconstructed path formula preserves algebraic section bounds
            # from the underlying CAD and is more faithful than recomposing from
            # the cached scalar fields when root functions are involved.
            return path_condition(
                self.source_cell,
                self.variables,
                _cells_by_level_from_source(self.source_cell),
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


# Weak map from underlying CADCell identity to its cells-by-level mapping. This
# keeps StructuredCADCell.as_formula useful even though the public dataclass keeps
# only the leaf cell by default.
_CELLS_BY_LEVEL_REGISTRY: dict[int, Mapping[int, Sequence[CADCell]]] = {}


def _cells_by_level_from_source(cell: CADCell) -> Mapping[int, Sequence[CADCell]]:
    return _CELLS_BY_LEVEL_REGISTRY.get(id(cell), {})


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
        specialized = sp.Poly(sp.expand(poly.subs(base_subs)), variable, domain="EX")
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
    _CELLS_BY_LEVEL_REGISTRY[id(leaf)] = cells_by_level
    return StructuredCADCell(
        variables=tuple(variables),
        index=leaf.index,
        levels=tuple(levels),
        sample=sample_map,
        source_cell=leaf,
        selected=selected,
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
        formula = _normalize_formula(condition_or_cad)
        if variables is None:
            variables = tuple(sorted(formula.free_symbols, key=lambda s: s.name))
        result = cad(
            formula, _normalize_variables(variables, formula), output="cells", return_result=True
        )

    vars_ = tuple(result.variables)
    if selected_only:
        leaves = tuple(result.cell_set.cells)
    else:
        leaves = tuple(result.cad.cells_by_level.get(len(vars_), ()))
    selected_indices = {cell.index for cell in result.cell_set.cells}
    sign_invariance = result.cad.verify_sign_invariance()
    structured = tuple(
        _structured_cell_from_leaf(
            cell,
            vars_,
            result.cad.cells_by_level,
            selected=cell.index in selected_indices,
            sign_invariant=sign_invariance.ok,
        )
        for cell in leaves
    )
    return StructuredCADCellDecomposition(vars_, structured, result.formula, result)


@dataclass(frozen=True)
class CylindricalCoordinateConstraint:
    """One coordinate constraint in a cylindrical solution cell.

    The constraint is interpreted relative to the preceding coordinates in the
    same cell. For example, in variables ``(x, y, z)`` the third level may have
    bounds involving both ``x`` and ``y``. Sections represent equalities
    ``variable == lower == upper``; sectors represent open intervals by default
    and may be rendered with weak inequalities when ``closed=True`` is requested.
    """

    variable: sp.Symbol
    level: int
    kind: str
    lower: sp.Expr
    upper: sp.Expr
    sample: sp.Expr
    index: tuple[int, ...]
    lower_bound: CADBound | None = None
    upper_bound: CADBound | None = None
    lower_closed: bool = False
    upper_closed: bool = False
    delineability: DelineabilityCertificate | None = None
    root_order: RootOrderCertificate | None = None

    @property
    def typed_lower(self) -> CADBound:
        return self.lower_bound or as_cad_bound(self.lower, closed=self.lower_closed)

    @property
    def typed_upper(self) -> CADBound:
        return self.upper_bound or as_cad_bound(self.upper, closed=self.upper_closed)

    @property
    def is_section(self) -> bool:
        return self.kind == "section"

    @property
    def is_sector(self) -> bool:
        return self.kind == "sector"

    @property
    def dimension(self) -> int:
        return 0 if self.is_section else 1

    def as_formula(self, *, closed: bool = False) -> sp.Expr:
        if self.is_section:
            return sp.Eq(self.variable, self.lower)
        parts: list[sp.Expr] = []
        if self.lower != -sp.oo:
            lower_closed = closed or self.lower_closed
            parts.append(
                self.variable >= self.lower if lower_closed else self.variable > self.lower
            )
        if self.upper != sp.oo:
            upper_closed = closed or self.upper_closed
            parts.append(
                self.variable <= self.upper if upper_closed else self.variable < self.upper
            )
        return sp.And(*parts) if parts else sp.true

    def as_limit(self) -> tuple[sp.Symbol, sp.Expr, sp.Expr]:
        """Return ``(variable, lower, upper)`` for iterated-integral style use."""

        return (self.variable, self.lower, self.upper)


@dataclass(frozen=True)
class CylindricalSolutionCell:
    """One nested cylindrical solution cell.

    A cell is an ordered path through a CAD: first a constraint on ``x1``, then
    a constraint on ``x2`` whose bounds may depend on ``x1``, and so on. This
    is the user-facing representation for solutions such as ``x`` in a base
    cell and ``y`` between algebraic functions over that base cell.
    """

    variables: tuple[sp.Symbol, ...]
    levels: tuple[CylindricalCoordinateConstraint, ...]
    sample: Mapping[sp.Symbol, sp.Expr]
    index: tuple[int, ...]
    selected: bool = True
    source_cell: StructuredCADCell | None = None

    @property
    def dimension(self) -> int:
        return sum(level.dimension for level in self.levels)

    @property
    def codimension(self) -> int:
        return len(self.variables) - self.dimension

    @property
    def bounded(self) -> bool:
        return all(level.lower != -sp.oo and level.upper != sp.oo for level in self.levels)

    @property
    def is_full_dimensional(self) -> bool:
        return self.dimension == len(self.variables)

    def as_formula(self, *, closed: bool = False) -> sp.Expr:
        parts = [level.as_formula(closed=closed) for level in self.levels]
        return sp.And(*parts) if parts else sp.true

    def sample_point(self) -> dict[sp.Symbol, sp.Expr]:
        return {sym: sp.simplify(value) for sym, value in self.sample.items()}

    def iterated_limits(self) -> tuple[tuple[sp.Symbol, sp.Expr, sp.Expr], ...]:
        """Return nested expression limits in CAD variable order."""

        return tuple(level.as_limit() for level in self.levels)

    def cylindrical_bounds(self) -> tuple[tuple[sp.Symbol, CADBound, CADBound], ...]:
        """Return typed nested bounds preserving algebraic root functions."""

        return tuple(
            (level.variable, level.typed_lower, level.typed_upper) for level in self.levels
        )

    def verify_bounds(self):
        from .bounds import verify_cad_cell_bounds

        return verify_cad_cell_bounds(self)


@dataclass(frozen=True)
class CylindricalDecompositionCertificate:
    """Certificate that solution cells form a complete disjoint decomposition."""

    coverage_verified: bool
    pairwise_disjoint: bool
    cells_verified: bool
    source: str = "unknown"
    notes: tuple[str, ...] = ()

    @property
    def certified(self) -> bool:
        return self.coverage_verified and self.pairwise_disjoint and self.cells_verified

    def verify(self) -> bool:
        return self.certified


@dataclass(frozen=True)
class CylindricalSolution:
    """Full cylindrical solution representation for a semialgebraic formula."""

    variables: tuple[sp.Symbol, ...]
    cells: tuple[CylindricalSolutionCell, ...]
    formula: sp.Expr
    source_decomposition: StructuredCADCellDecomposition | None = None
    decomposition_cert: CylindricalDecompositionCertificate | None = None

    @property
    def dimension(self) -> int | None:
        if not self.cells:
            return None
        return max(cell.dimension for cell in self.cells)

    @property
    def bounded(self) -> bool | None:
        if not self.cells:
            return True
        return all(cell.bounded for cell in self.cells)

    @property
    def full_dimensional_cells(self) -> tuple[CylindricalSolutionCell, ...]:
        n = len(self.variables)
        return tuple(cell for cell in self.cells if cell.dimension == n)

    @property
    def sample_points(self) -> tuple[Mapping[sp.Symbol, sp.Expr], ...]:
        return tuple(cell.sample_point() for cell in self.cells)

    def cells_by_dimension(self, dimension: int) -> tuple[CylindricalSolutionCell, ...]:
        return tuple(cell for cell in self.cells if cell.dimension == dimension)

    def as_formula(self, *, closed: bool = False) -> sp.Expr:
        if not self.cells:
            return sp.false
        pieces = [cell.as_formula(closed=closed) for cell in self.cells]
        return sp.Or(*pieces) if len(pieces) > 1 else pieces[0]


def _cylindrical_cell_from_structured(cell: StructuredCADCell) -> CylindricalSolutionCell:
    levels = tuple(
        CylindricalCoordinateConstraint(
            variable=level.variable,
            level=level.level,
            kind=level.kind,
            lower=level.lower,
            upper=level.upper,
            sample=level.sample,
            index=level.index,
            lower_bound=level.lower_bound,
            upper_bound=level.upper_bound,
            lower_closed=level.lower_closed,
            upper_closed=level.upper_closed,
            delineability=level.delineability,
            root_order=level.root_order,
        )
        for level in cell.levels
    )
    return CylindricalSolutionCell(
        variables=cell.variables,
        levels=levels,
        sample=cell.sample,
        index=cell.index,
        selected=cell.selected,
        source_cell=cell,
    )


def cylindrical_solution_from_structured(
    decomposition: StructuredCADCellDecomposition,
) -> CylindricalSolution:
    """Convert structured CAD cells to the public cylindrical solution form."""

    cells = tuple(_cylindrical_cell_from_structured(cell) for cell in decomposition.cells)
    from .bounds import verify_cad_cell_bounds

    cert = CylindricalDecompositionCertificate(
        coverage_verified=True,
        pairwise_disjoint=True,
        cells_verified=all(verify_cad_cell_bounds(cell).verify() for cell in cells),
        source="cad",
        notes=("CAD leaf cells form a cylindrical partition of the selected formula",),
    )
    return CylindricalSolution(
        variables=decomposition.variables,
        cells=cells,
        formula=decomposition.formula,
        source_decomposition=decomposition,
        decomposition_cert=cert,
    )


def extract_cylindrical_solution(
    condition_or_cad: object,
    variables: Sequence[sp.Symbol | str] | None = None,
    *,
    selected_only: bool = True,
) -> CylindricalSolution:
    """Return a nested cylindrical representation of the selected solution cells.

    This is the high-level public form of the CAD cell extraction layer. It can
    represent arbitrary CAD paths in any dimension: each cell stores coordinate
    bounds in variable order, and bounds at level ``k`` may involve variables
    from levels ``< k``.
    """

    if not isinstance(condition_or_cad, (CADResult, CellSet)):
        formula = _normalize_formula(condition_or_cad)
        vars_ = (
            _normalize_variables(variables, formula)
            if variables is not None
            else tuple(sorted(formula.free_symbols, key=lambda symbol: symbol.name))
        )
        explicit = extract_explicit_cylindrical_solution(formula, vars_)
        if explicit is not None and any(
            isinstance(bound, AlgebraicRootFunction)
            for cell in explicit.cells
            for level in cell.levels
            for bound in (level.lower_bound, level.upper_bound)
            if bound is not None
        ):
            return explicit
    decomp = extract_structured_cad_cells(condition_or_cad, variables, selected_only=selected_only)
    return cylindrical_solution_from_structured(decomp)


def _sample_between_bounds(
    lower: sp.Expr, upper: sp.Expr, previous: Mapping[sp.Symbol, sp.Expr]
) -> sp.Expr:
    lo = sp.simplify(lower.subs(previous)) if hasattr(lower, "subs") else lower
    hi = sp.simplify(upper.subs(previous)) if hasattr(upper, "subs") else upper
    if lo == -sp.oo and hi == sp.oo:
        return sp.Integer(0)
    if lo == -sp.oo:
        return sp.simplify(hi - 1)
    if hi == sp.oo:
        return sp.simplify(lo + 1)
    if sp.simplify(hi - lo) == 0:
        return lo
    return sp.simplify((lo + hi) / 2)


def _relational_bound_for_variable(
    atom: sp.Expr, variable: sp.Symbol
) -> tuple[str, sp.Expr, bool] | None:
    """Return ("lower"|"upper"|"equal", bound) for simple explicit atoms."""

    if not isinstance(atom, sp.core.relational.Relational):
        return None
    lhs, rhs = atom.lhs, atom.rhs
    if lhs == variable and not rhs.has(variable):
        if isinstance(atom, (sp.GreaterThan, sp.StrictGreaterThan)):
            return ("lower", rhs, isinstance(atom, sp.GreaterThan))
        if isinstance(atom, (sp.LessThan, sp.StrictLessThan)):
            return ("upper", rhs, isinstance(atom, sp.LessThan))
        if isinstance(atom, sp.Equality):
            return ("equal", rhs, True)
    if rhs == variable and not lhs.has(variable):
        if isinstance(atom, (sp.GreaterThan, sp.StrictGreaterThan)):
            return ("upper", lhs, isinstance(atom, sp.GreaterThan))
        if isinstance(atom, (sp.LessThan, sp.StrictLessThan)):
            return ("lower", lhs, isinstance(atom, sp.LessThan))
        if isinstance(atom, sp.Equality):
            return ("equal", lhs, True)
    return None


def _and_atoms(expr: sp.Expr) -> tuple[sp.Expr, ...]:
    if expr is sp.true or expr == sp.true:
        return ()
    if isinstance(expr, sp.And):
        out: list[sp.Expr] = []
        for arg in expr.args:
            out.extend(_and_atoms(arg))
        return tuple(out)
    return (expr,)


def _or_pieces(expr: sp.Expr) -> tuple[sp.Expr, ...]:
    if isinstance(expr, sp.Or):
        return tuple(expr.args)
    return (expr,)


def _merge_lower(a: sp.Expr, b: sp.Expr) -> sp.Expr:
    if a == -sp.oo:
        return b
    if b == -sp.oo:
        return a
    return sp.Max(a, b)


def _merge_upper(a: sp.Expr, b: sp.Expr) -> sp.Expr:
    if a == sp.oo:
        return b
    if b == sp.oo:
        return a
    return sp.Min(a, b)


def _compare_expr(a: sp.Expr, b: sp.Expr) -> int | None:
    """Return -1/0/1 when ``a`` is provably below/equal/above ``b``."""

    diff = sp.simplify(sp.sympify(a) - sp.sympify(b))
    if diff == 0:
        return 0
    if diff.is_positive is True:
        return 1
    if diff.is_negative is True:
        return -1
    return None


def _affine_extreme_over_levels(
    expr: sp.Expr,
    levels: Sequence[CylindricalCoordinateConstraint],
    *,
    minimum: bool,
) -> sp.Expr | None:
    """Return an exact affine extreme over an established triangular cell.

    The routine deliberately handles only expressions affine in each previous
    coordinate with a coefficient whose sign is already decidable.  It is used
    by the explicit cylindrical fast path to prove dependent bounds such as
    ``0 <= x + y`` from earlier bounds ``x >= 0`` and ``y >= 0``.  Unsupported
    nonlinear or sign-indeterminate cases return ``None`` and therefore fall
    back to full CAD.
    """

    value = sp.expand(sp.sympify(expr))
    for level in reversed(tuple(levels)):
        var = level.variable
        if var not in value.free_symbols:
            continue
        try:
            poly = sp.Poly(value, var, domain="EX")
        except sp.PolynomialError:
            return None
        if poly.degree() > 1:
            return None
        coeff = sp.simplify(poly.coeff_monomial(var))
        if var in coeff.free_symbols:
            return None
        sign = _compare_expr(coeff, sp.Integer(0))
        if sign is None:
            return None
        if sign == 0:
            value = sp.simplify(value.subs(var, 0))
            continue
        choose_lower = (minimum and sign > 0) or (not minimum and sign < 0)
        bound = level.lower if choose_lower else level.upper
        if bound in (-sp.oo, sp.oo):
            return None
        value = sp.expand(value.subs(var, bound))
    return sp.simplify(value)


def _compare_over_base_cell(
    left: sp.Expr,
    right: sp.Expr,
    levels: Sequence[CylindricalCoordinateConstraint],
) -> int | None:
    """Compare bounds exactly, using prior triangular bounds when needed."""

    direct = _compare_expr(left, right)
    if direct is not None:
        return direct
    diff = sp.expand(sp.sympify(right) - sp.sympify(left))
    minimum = _affine_extreme_over_levels(diff, levels, minimum=True)
    if minimum is not None:
        min_cmp = _compare_expr(minimum, 0)
        if min_cmp is not None and min_cmp >= 0:
            maximum = _affine_extreme_over_levels(diff, levels, minimum=False)
            if maximum is not None and _compare_expr(maximum, 0) == 0:
                return 0
            return -1
    maximum = _affine_extreme_over_levels(diff, levels, minimum=False)
    if maximum is not None:
        max_cmp = _compare_expr(maximum, 0)
        if max_cmp is not None and max_cmp < 0:
            return 1
    return None


def _cell_pair_provably_disjoint(a: CylindricalSolutionCell, b: CylindricalSolutionCell) -> bool:
    # A single coordinate separation is enough.  This deliberately declines
    # dependent-bound cases whose separation is not globally provable.
    for la, lb in zip(a.levels, b.levels, strict=True):
        if la.variable != lb.variable:
            return False
        cmp_ab = _compare_expr(la.upper, lb.lower)
        if cmp_ab == -1 or (cmp_ab == 0 and not (la.upper_closed and lb.lower_closed)):
            return True
        cmp_ba = _compare_expr(lb.upper, la.lower)
        if cmp_ba == -1 or (cmp_ba == 0 and not (lb.upper_closed and la.lower_closed)):
            return True
    return False


def _cells_pairwise_disjoint(cells: Sequence[CylindricalSolutionCell]) -> bool:
    return all(
        _cell_pair_provably_disjoint(cells[i], cells[j])
        for i in range(len(cells))
        for j in range(i + 1, len(cells))
    )


def _coefficient_sign_over_levels(
    expr: sp.Expr,
    levels: Sequence[CylindricalCoordinateConstraint],
) -> int | None:
    """Prove the sign of an expression over an established triangular base cell.

    The helper is intentionally conservative.  It first uses SymPy's exact sign
    information, then the existing affine-extreme machinery.  Unsupported
    nonlinear coefficient dependencies are declined rather than approximated.
    """

    value = sp.simplify(sp.sympify(expr))
    direct = _compare_expr(value, sp.Integer(0))
    if direct is not None:
        return direct
    minimum = _affine_extreme_over_levels(value, levels, minimum=True)
    maximum = _affine_extreme_over_levels(value, levels, minimum=False)
    if minimum is not None:
        minimum_sign = _compare_expr(minimum, sp.Integer(0))
        if minimum_sign is not None and minimum_sign >= 0:
            if maximum is not None and _compare_expr(maximum, sp.Integer(0)) == 0:
                return 0
            return 1
    if maximum is not None:
        maximum_sign = _compare_expr(maximum, sp.Integer(0))
        if maximum_sign is not None and maximum_sign <= 0:
            if minimum is not None and _compare_expr(minimum, sp.Integer(0)) == 0:
                return 0
            return -1
    return None


def _global_monotonicity_sign(
    polynomial: sp.Poly,
    variable: sp.Symbol,
    levels: Sequence[CylindricalCoordinateConstraint],
) -> int | None:
    """Return ``1``/``-1`` for globally increasing/decreasing odd polynomials.

    We only certify derivatives whose nonzero powers are even and whose
    coefficients have a common proven sign over the preceding cylindrical
    levels.  This covers forms such as ``z**3 + x*z + y`` on ``x >= 0`` while
    deliberately declining sign-indefinite or more complicated cases.
    """

    if polynomial.degree() < 1 or polynomial.degree() % 2 == 0:
        return None
    derivative = sp.Poly(sp.diff(polynomial.as_expr(), variable), variable, domain="EX")
    if derivative.is_zero:
        return None
    signs: list[int] = []
    for (power,), coefficient in derivative.terms():
        if power % 2:
            return None
        sign = _coefficient_sign_over_levels(coefficient, levels)
        if sign is None:
            return None
        if sign != 0:
            signs.append(sign)
    if not signs:
        return None
    if all(sign > 0 for sign in signs):
        return 1
    if all(sign < 0 for sign in signs):
        return -1
    return None


def _implicit_monotone_root_bound(
    atom: sp.Expr,
    variable: sp.Symbol,
    levels: Sequence[CylindricalCoordinateConstraint],
    sample_map: Mapping[sp.Symbol, sp.Expr],
) -> tuple[str, AlgebraicRootFunction, bool] | None:
    """Convert a certified monotone polynomial relation into one root bound.

    The polynomial must have odd degree in ``variable`` and a derivative whose
    sign is globally certified by :func:`_global_monotonicity_sign`.  Such a
    polynomial has exactly one real root in every base fiber, so its relational
    atom is equivalent to a single lower/upper bound (or section).
    """

    if not isinstance(
        atom,
        (
            sp.Equality,
            sp.LessThan,
            sp.StrictLessThan,
            sp.GreaterThan,
            sp.StrictGreaterThan,
        ),
    ):
        return None
    residual = sp.expand(atom.lhs - atom.rhs)
    try:
        polynomial = sp.Poly(residual, variable, domain="EX")
    except sp.PolynomialError:
        return None
    if polynomial.degree() <= 1:
        return None
    previous_variables = tuple(level.variable for level in levels)
    if (residual.free_symbols - {variable}) - set(previous_variables):
        return None
    monotonicity = _global_monotonicity_sign(polynomial, variable, levels)
    if monotonicity is None:
        return None

    specialized = sp.expand(residual.subs(dict(sample_map)))
    try:
        sample_roots = tuple(sp.real_roots(sp.Poly(specialized, variable).as_expr()))
    except (sp.PolynomialError, ValueError, TypeError, NotImplementedError):
        return None
    if len(sample_roots) != 1:
        return None
    sample_root = sp.sympify(sample_roots[0])
    sample_verified = sp.simplify(specialized.subs(variable, sample_root)) == 0
    if not sample_verified:
        return None

    base_variables = tuple(
        level.variable for level in levels if level.variable in (residual.free_symbols - {variable})
    )
    certificate = DelineabilityCertificate(
        polynomial=residual,
        fiber_variable=variable,
        root_index=0,
        base_variables=base_variables,
        sign_invariant=True,
        stack_order_verified=True,
        sample_root_verified=True,
        sample_root_value=sample_root,
        representation_verified=True,
        notes=("unique real root certified by global monotonicity over the explicit base cell",),
    )
    root = AlgebraicRootFunction(
        polynomial=residual,
        fiber_variable=variable,
        root_index=0,
        base_variables=base_variables,
        certificate=certificate,
    )

    if isinstance(atom, sp.Equality):
        return ("equal", root, True)
    relation_is_less = isinstance(atom, (sp.LessThan, sp.StrictLessThan))
    closed = isinstance(atom, (sp.LessThan, sp.GreaterThan))
    if monotonicity > 0:
        kind = "upper" if relation_is_less else "lower"
    else:
        kind = "lower" if relation_is_less else "upper"
    return (kind, root, closed)


def _explicit_cylindrical_cell_from_conjunction(
    expr: sp.Expr,
    variables: Sequence[sp.Symbol | str],
    *,
    cell_index: int = 0,
) -> CylindricalSolutionCell | None:
    """Convert a provably cylindrical conjunction into one typed solution cell.

    The shortcut declines when symbolic bound ordering cannot be certified, so
    callers can fall back to the general CAD decomposition.
    """
    vars_ = tuple(sp.Symbol(v, real=True) if isinstance(v, str) else v for v in variables)
    atoms = list(_and_atoms(expr))
    levels: list[CylindricalCoordinateConstraint] = []
    sample_map: dict[sp.Symbol, sp.Expr] = {}
    used: set[int] = set()
    for level, var in enumerate(vars_, start=1):
        later = set(vars_[level:])
        lower: sp.Expr = -sp.oo
        upper: sp.Expr = sp.oo
        lower_closed = False
        upper_closed = False
        lower_typed: CADBound | None = None
        upper_typed: CADBound | None = None
        for i, atom in enumerate(atoms):
            bound = _relational_bound_for_variable(atom, var)
            if bound is None:
                continue
            kind, value, is_closed = bound
            if any(value.has(sym) for sym in later):
                continue
            new_value = sp.sympify(value)
            if kind == "lower":
                if lower == -sp.oo:
                    lower, lower_closed = new_value, is_closed
                else:
                    cmp = _compare_over_base_cell(lower, new_value, levels)
                    if cmp is not None:
                        cmp = -cmp
                    if cmp is None:
                        return None
                    if cmp > 0:
                        lower, lower_closed = new_value, is_closed
                    elif cmp == 0:
                        lower_closed = lower_closed and is_closed
            elif kind == "upper":
                if upper == sp.oo:
                    upper, upper_closed = new_value, is_closed
                else:
                    cmp = _compare_over_base_cell(new_value, upper, levels)
                    if cmp is None:
                        return None
                    if cmp < 0:
                        upper, upper_closed = new_value, is_closed
                    elif cmp == 0:
                        upper_closed = upper_closed and is_closed
            else:
                if lower != -sp.oo:
                    cmp = _compare_over_base_cell(lower, new_value, levels)
                    if cmp is not None:
                        cmp = -cmp
                    if cmp is None or cmp < 0 or (cmp == 0 and not lower_closed):
                        return None
                if upper != sp.oo:
                    cmp = _compare_over_base_cell(new_value, upper, levels)
                    if cmp is None or cmp > 0 or (cmp == 0 and not upper_closed):
                        return None
                lower = upper = new_value
                lower_closed = upper_closed = True
            used.add(i)
        for i, atom in enumerate(atoms):
            if i in used:
                continue
            if not atom.has(var) or any(atom.has(sym) for sym in later):
                continue
            implicit = _implicit_monotone_root_bound(atom, var, levels, sample_map)
            if implicit is None:
                return None
            kind, root, is_closed = implicit
            root_expr = root.as_expr()
            if kind == "lower":
                if lower != -sp.oo:
                    return None
                lower = root_expr
                lower_closed = is_closed
                lower_typed = as_cad_bound(root, closed=is_closed)
            elif kind == "upper":
                if upper != sp.oo:
                    return None
                upper = root_expr
                upper_closed = is_closed
                upper_typed = as_cad_bound(root, closed=is_closed)
            else:
                if lower != -sp.oo or upper != sp.oo:
                    return None
                lower = upper = root_expr
                lower_closed = upper_closed = True
                lower_typed = as_cad_bound(root, closed=True)
                upper_typed = as_cad_bound(root, closed=True)
            used.add(i)

        if lower != -sp.oo and upper != sp.oo:
            ordering = _compare_over_base_cell(lower, upper, levels)
            if ordering is None:
                return None
            if ordering > 0:
                return None
            if ordering == 0 and not (lower_closed and upper_closed):
                return None

        if isinstance(lower_typed, AlgebraicRootFunction) and upper == sp.oo:
            sample = sp.simplify(lower_typed.certificate.sample_root_value + 1)
        elif isinstance(upper_typed, AlgebraicRootFunction) and lower == -sp.oo:
            sample = sp.simplify(upper_typed.certificate.sample_root_value - 1)
        elif (
            isinstance(lower_typed, AlgebraicRootFunction)
            and isinstance(upper_typed, AlgebraicRootFunction)
            and lower == upper
        ):
            sample = sp.sympify(lower_typed.certificate.sample_root_value)
        else:
            sample = _sample_between_bounds(lower, upper, sample_map)
        sample_map[var] = sample
        cell_kind = "section" if _compare_over_base_cell(lower, upper, levels) == 0 else "sector"
        lower_expr = sp.simplify(lower)
        upper_expr = sp.simplify(upper)
        levels.append(
            CylindricalCoordinateConstraint(
                variable=var,
                level=level,
                kind=cell_kind,
                lower=lower_expr,
                upper=upper_expr,
                sample=sp.simplify(sample),
                index=tuple([cell_index + 1] * level),
                lower_bound=lower_typed or as_cad_bound(lower_expr, closed=lower_closed),
                upper_bound=upper_typed or as_cad_bound(upper_expr, closed=upper_closed),
                lower_closed=lower_closed,
                upper_closed=upper_closed,
            )
        )
    return CylindricalSolutionCell(
        variables=vars_,
        levels=tuple(levels),
        sample=sample_map,
        index=tuple([cell_index + 1] * len(vars_)),
        selected=True,
        source_cell=None,
    )


def extract_explicit_cylindrical_solution(
    condition: object,
    variables: Sequence[sp.Symbol | str],
) -> CylindricalSolution | None:
    """Extract explicit nested bounds without invoking CAD.

    This recognizes cylindrical formulas such as ``0 <= x <= 1``,
    ``0 <= y <= x``, ``0 <= z <= y`` in arbitrary dimension. It is a fast
    companion to the complete CAD extractor and lets high-dimensional APIs use
    nested bounds when the user already supplied them explicitly.
    """

    formula = _normalize_formula(condition)
    vars_ = tuple(sp.Symbol(v, real=True) if isinstance(v, str) else v for v in variables)
    if formula is sp.false or formula == sp.false:
        cert = CylindricalDecompositionCertificate(
            coverage_verified=True,
            pairwise_disjoint=True,
            cells_verified=True,
            source="explicit-cylindrical",
            notes=("normalized formula is empty",),
        )
        return CylindricalSolution(vars_, (), formula, None, cert)
    cells: list[CylindricalSolutionCell] = []
    pieces = _or_pieces(formula)
    for i, piece in enumerate(pieces):
        cell = _explicit_cylindrical_cell_from_conjunction(piece, vars_, cell_index=i)
        if cell is None:
            return None
        cells.append(cell)
    if len(cells) > 1 and not _cells_pairwise_disjoint(cells):
        return None
    from .bounds import verify_cad_cell_bounds

    cells_verified = all(verify_cad_cell_bounds(cell).verify() for cell in cells)
    cert = CylindricalDecompositionCertificate(
        coverage_verified=True,
        pairwise_disjoint=True,
        cells_verified=cells_verified,
        source="explicit-cylindrical",
        notes=("each Or piece was preserved exactly and pairwise disjointness was proven",),
    )
    return CylindricalSolution(
        variables=vars_,
        cells=tuple(cells),
        formula=formula,
        source_decomposition=None,
        decomposition_cert=cert,
    )


def structured_cad_cells_to_vertical_bounds_2d(
    cells: Sequence[StructuredCADCell],
) -> tuple[VerticalBoundCell2D, ...]:
    """Convert two-dimensional structured CAD cells to vertical bounds."""

    out: list[VerticalBoundCell2D] = []
    for cell in cells:
        if len(cell.variables) != 2:
            raise ValueError("vertical-bound conversion requires two-dimensional cells")
        x, y = cell.variables
        if len(cell.levels) != 2:
            continue
        x_level, y_level = cell.levels
        y_bounds = ((y_level.lower, y_level.upper),)
        out.append(
            VerticalBoundCell2D(
                x_variable=x,
                y_variable=y,
                x_interval=(x_level.lower, x_level.upper),
                y_bounds=y_bounds,
                x_condition=x_level.condition,
                source_formula=cell.as_formula(closed=False)
                if cell.source_cell is not None
                else None,
            )
        )
    return tuple(out)


def extract_vertical_bounds_from_cad_2d(
    condition_or_cad: object,
    variables: Sequence[sp.Symbol | str] | None = None,
    *,
    full_dimensional_only: bool = False,
) -> tuple[VerticalBoundCell2D, ...]:
    """Extract 2D vertical-bound cells from a full CAD decomposition."""

    decomposition = extract_structured_cad_cells(condition_or_cad, variables, selected_only=True)
    cells = decomposition.full_dimensional_cells if full_dimensional_only else decomposition.cells
    return structured_cad_cells_to_vertical_bounds_2d(cells)


__all__ = [
    "CADBound",
    "AlgebraicRootFunction",
    "CertifiedRootComparison",
    "DelineabilityCertificate",
    "RootOrderCertificate",
    "CADCellBoundsCertificate",
    "verify_cad_cell_bounds",
    "CADCellIntegral",
    "IntrinsicCellStratum",
    "IntrinsicStratification",
    "stratify_intrinsic_solution",
    "full_dimensional_cell_integral",
    "full_dimensional_solution_integrals",
    "intrinsic_cell_integral",
    "intrinsic_solution_integrals",
    "CylindricalCoordinateConstraint",
    "CylindricalSolutionCell",
    "CylindricalSolution",
    "CylindricalDecompositionCertificate",
    "cylindrical_solution_from_structured",
    "extract_cylindrical_solution",
    "extract_explicit_cylindrical_solution",
    "StructuredCADLevel",
    "StructuredCADCell",
    "StructuredCADCellDecomposition",
    "extract_structured_cad_cells",
    "structured_cad_cells_to_vertical_bounds_2d",
    "extract_vertical_bounds_from_cad_2d",
]

# Public convenience re-exports kept here so the package-level lazy CAD-cell
# namespace can expose typed bounds and cell integration without extra dispatch
# tables.
from .bounds import CADCellBoundsCertificate, verify_cad_cell_bounds  # noqa: E402
from .integration import (  # noqa: E402
    CADCellIntegral,
    IntrinsicCellStratum,
    IntrinsicStratification,
    full_dimensional_cell_integral,
    full_dimensional_solution_integrals,
    intrinsic_cell_integral,
    intrinsic_solution_integrals,
    stratify_intrinsic_solution,
)
