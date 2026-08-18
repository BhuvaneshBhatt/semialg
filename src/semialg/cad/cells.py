from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass

import sympy as sp

from ..algebraic.samples import sample_to_expr
from ..decomposition.cylindrical import CADResult, CellSet, cad
from ..implicit_utils import VerticalBoundCell2D, _normalize_formula, _normalize_variables
from ..reconstruct.cylindrical import path_condition, section_value_expr
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


def _level_bounds(
    level_cell: CADCell, variable: sp.Symbol, cells_by_level: Mapping[int, Sequence[CADCell]]
) -> tuple[sp.Expr, sp.Expr]:
    if level_cell.kind == "section":
        value = section_value_expr(level_cell, variable)
        return sp.simplify(value), sp.simplify(value)
    left_section = _section_at_position(
        cells_by_level, level_cell.level, level_cell.parent_index, level_cell.stack_position - 1
    )
    right_section = _section_at_position(
        cells_by_level, level_cell.level, level_cell.parent_index, level_cell.stack_position + 1
    )
    if left_section is not None:
        lower = section_value_expr(left_section, variable)
    elif level_cell.lower_bound is not None:
        lower = sample_to_expr(level_cell.lower_bound)
    else:
        lower = -sp.oo
    if right_section is not None:
        upper = section_value_expr(right_section, variable)
    elif level_cell.upper_bound is not None:
        upper = sample_to_expr(level_cell.upper_bound)
    else:
        upper = sp.oo
    return sp.simplify(lower), sp.simplify(upper)


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
) -> StructuredCADCell:
    levels: list[StructuredCADLevel] = []
    sample_map: dict[sp.Symbol, sp.Expr] = {}
    for level, variable in enumerate(variables, start=1):
        prefix = leaf.index[:level]
        level_cell = next(cell for cell in cells_by_level[level] if cell.index == prefix)
        lower, upper = _level_bounds(level_cell, variable, cells_by_level)
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

    Unlike the earlier vertical-slice parser, this function uses the package's
    actual CAD lifting data. It therefore works for arbitrary formulas that the
    complete CAD engine can decompose, including algebraic stack bounds such as
    ``-sqrt(x) < y < sqrt(x)`` produced by a cell over ``y**2 - x``.
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
    structured = tuple(
        _structured_cell_from_leaf(
            cell, vars_, result.cad.cells_by_level, selected=cell.index in selected_indices
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
        return _level_condition(self.variable, self.kind, self.lower, self.upper, closed=closed)

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
        """Return nested ``(var, lower, upper)`` limits in CAD variable order."""

        return tuple(level.as_limit() for level in self.levels)


@dataclass(frozen=True)
class CylindricalSolution:
    """Full cylindrical solution representation for a semialgebraic formula."""

    variables: tuple[sp.Symbol, ...]
    cells: tuple[CylindricalSolutionCell, ...]
    formula: sp.Expr
    source_decomposition: StructuredCADCellDecomposition | None = None

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

    return CylindricalSolution(
        variables=decomposition.variables,
        cells=tuple(_cylindrical_cell_from_structured(cell) for cell in decomposition.cells),
        formula=decomposition.formula,
        source_decomposition=decomposition,
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
) -> tuple[str, sp.Expr] | None:
    """Return ("lower"|"upper"|"equal", bound) for simple explicit atoms."""

    if not isinstance(atom, sp.core.relational.Relational):
        return None
    lhs, rhs = atom.lhs, atom.rhs
    if lhs == variable and not rhs.has(variable):
        if isinstance(atom, (sp.GreaterThan, sp.StrictGreaterThan)):
            return ("lower", rhs)
        if isinstance(atom, (sp.LessThan, sp.StrictLessThan)):
            return ("upper", rhs)
        if isinstance(atom, sp.Equality):
            return ("equal", rhs)
    if rhs == variable and not lhs.has(variable):
        if isinstance(atom, (sp.GreaterThan, sp.StrictGreaterThan)):
            return ("upper", lhs)
        if isinstance(atom, (sp.LessThan, sp.StrictLessThan)):
            return ("lower", lhs)
        if isinstance(atom, sp.Equality):
            return ("equal", lhs)
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


def _explicit_cylindrical_cell_from_conjunction(
    expr: sp.Expr,
    variables: Sequence[sp.Symbol | str],
    *,
    cell_index: int = 0,
) -> CylindricalSolutionCell | None:
    vars_ = tuple(sp.Symbol(v, real=True) if isinstance(v, str) else v for v in variables)
    atoms = list(_and_atoms(expr))
    levels: list[CylindricalCoordinateConstraint] = []
    sample_map: dict[sp.Symbol, sp.Expr] = {}
    used: set[int] = set()
    for level, var in enumerate(vars_, start=1):
        later = set(vars_[level:])
        lower: sp.Expr = -sp.oo
        upper: sp.Expr = sp.oo
        for i, atom in enumerate(atoms):
            bound = _relational_bound_for_variable(atom, var)
            if bound is None:
                continue
            kind, value = bound
            # Bounds may depend only on earlier variables for a cylindrical
            # cell. If the apparent bound contains a later variable, this atom
            # actually belongs to a later coordinate level (e.g. y <= x while
            # processing x), so leave it for that level instead of rejecting the
            # whole explicit decomposition.
            if any(value.has(sym) for sym in later):
                continue
            if kind == "lower":
                lower = _merge_lower(lower, sp.sympify(value))
            elif kind == "upper":
                upper = _merge_upper(upper, sp.sympify(value))
            else:
                lower = upper = sp.sympify(value)
            used.add(i)
        # Unsupported atoms involving this variable cannot be represented by
        # this explicit syntactic path; leave them to the complete CAD extractor.
        for i, atom in enumerate(atoms):
            if i in used:
                continue
            if atom.has(var) and not any(atom.has(sym) for sym in later):
                return None
        sample = _sample_between_bounds(lower, upper, sample_map)
        sample_map[var] = sample
        cell_kind = "section" if lower == upper else "sector"
        levels.append(
            CylindricalCoordinateConstraint(
                variable=var,
                level=level,
                kind=cell_kind,
                lower=sp.simplify(lower),
                upper=sp.simplify(upper),
                sample=sp.simplify(sample),
                index=tuple([cell_index + 1] * level),
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
    cells: list[CylindricalSolutionCell] = []
    for i, piece in enumerate(_or_pieces(formula)):
        cell = _explicit_cylindrical_cell_from_conjunction(piece, vars_, cell_index=i)
        if cell is None:
            return None
        cells.append(cell)
    return CylindricalSolution(
        variables=vars_, cells=tuple(cells), formula=formula, source_decomposition=None
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
    "CylindricalCoordinateConstraint",
    "CylindricalSolutionCell",
    "CylindricalSolution",
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
