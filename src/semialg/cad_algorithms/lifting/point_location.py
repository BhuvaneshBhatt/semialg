"""Exact point location in cylindrical stacks with moving sections."""

from __future__ import annotations

import sympy as sp

from ...algebraic.comparison import compare_samples
from ...algebraic.samples import AlgebraicRoot, RationalSample, sample_to_expr
from ...exact_arithmetic import compare_exact_reals


class _CellLocator:
    """Index a CAD's stacks and lazily recover its section root selectors."""

    def __init__(self, variables, cells_by_level):
        self.variables = tuple(variables)
        self.stacks = {}
        for level, cells in cells_by_level.items():
            for cell in cells:
                self.stacks.setdefault((level, cell.parent_index), []).append(cell)
        for cells in self.stacks.values():
            cells.sort(key=lambda cell: cell.stack_position)
        self.section_bounds = {}
        self._last_location = ((), ())

    def locate(self, point):
        """Return the cell containing an exact point, including prefix points.

        Stored isolating intervals belong to the CAD sample fibers. At a new
        base point, delineability preserves section order but their values must
        be obtained by specializing the defining root selectors.
        """
        if len(point) > len(self.variables):
            raise ValueError("point has more coordinates than the CAD")
        if not self.variables:
            cells = self.stacks.get((0, None), ())
            return cells[0] if cells else None
        values = tuple(
            value if isinstance(value, (RationalSample, AlgebraicRoot)) else sp.sympify(value)
            for value in point
        )
        base_point = dict(zip(self.variables, values, strict=False))
        last_point, last_path = self._last_location
        shared = 0
        for value, previous in zip(values, last_point, strict=False):
            if not _same_coordinate(value, previous):
                break
            shared += 1
        # Refinement cells arrive in cylindrical order and share their parent
        # coordinates. Retain just the previous path, avoiding repeated root
        # comparisons without keeping an unbounded point-location cache.
        path = list(last_path[:shared])
        found = path[-1] if path else None
        parent_index = found.index if found is not None else None
        for level, value in enumerate(values[shared:], start=shared + 1):
            cells = self.stacks.get((level, parent_index), ())
            found = self._locate_in_stack(cells, value, base_point)
            if found is None:
                return None
            parent_index = found.index
            path.append(found)
        self._last_location = (values, tuple(path))
        return found

    def _locate_in_stack(self, cells, value, point):
        """Compare against each section in fiber order and retain intervening sectors."""
        from ...reconstruct.cylindrical import section_value_bound

        at_sample_base = bool(cells) and all(
            _same_coordinate(point[var], sample)
            for var, sample in zip(
                self.variables[: cells[0].level - 1], cells[0].sample[:-1], strict=True
            )
        )
        preceding_sector = None
        for cell in cells:
            if cell.kind == "sector":
                preceding_sector = cell
                continue
            boundary = cell.sample[-1]
            variable = self.variables[cell.level - 1]
            polynomial = cell.section_polynomial
            if (
                not at_sample_base
                and polynomial is not None
                and polynomial.free_symbols - {variable}
            ):
                bound = self.section_bounds.get(cell.index)
                if bound is None:
                    bound = section_value_bound(
                        cell, variable, base_variables=self.variables[: cell.level - 1]
                    )
                    self.section_bounds[cell.index] = bound
                boundary = bound.evaluate(point)
            comparison = _compare_coordinate(value, boundary)
            if comparison == 0:
                return cell
            if comparison < 0:
                return preceding_sector
            preceding_sector = None
        return preceding_sector


def _same_coordinate(value, sample):
    if value is sample:
        return True
    if isinstance(sample, RationalSample):
        return (value.value if isinstance(value, RationalSample) else value) == sample.value
    return False


def _compare_coordinate(value, sample):
    if isinstance(value, (RationalSample, AlgebraicRoot)):
        return compare_samples(value, sample)
    if value.is_Rational:
        return compare_samples(RationalSample(value), sample)
    return compare_exact_reals(value, sample_to_expr(sample))
