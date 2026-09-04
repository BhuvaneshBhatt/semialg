from __future__ import annotations

import pytest
import sympy as sp

from semialg import SemialgebraicRegion
from semialg.cad_algorithms.cells import extract_structured_cad_cells
from semialg.cad_region import triangulate_cad_cell


def test_contains_rejects_unresolved_symbolic_point():
    x, a = sp.symbols("x a", real=True)
    region = SemialgebraicRegion(sp.And(x >= 0, x <= 1), (x,))
    with pytest.raises(ValueError, match="membership remains symbolic"):
        region.contains((a,))


def test_point_location_rejects_unresolved_symbolic_coordinates():
    x, a = sp.symbols("x a", real=True)
    cad = SemialgebraicRegion(sp.And(x >= 0, x <= 1), (x,)).as_cad_region()
    with pytest.raises(ValueError, match="exact resolved coordinates"):
        cad.locate_point((a,))


def test_unbounded_structured_cell_triangulation_fails_conservatively():
    x, y = sp.symbols("x y", real=True)
    cad = SemialgebraicRegion(sp.And(x >= 0, y >= 0), (x, y)).as_cad_region()
    decomp = extract_structured_cad_cells(cad.result)
    unbounded = next(cell for cell in decomp.full_dimensional_cells if not cell.bounded)
    with pytest.raises(ValueError, match="finite bounds"):
        triangulate_cad_cell(unbounded)


def test_delineable_curve_rejects_disappearing_root_branch():
    from semialg.cad_region import evaluate_delineable_curve

    x, y, z = sp.symbols("x y z", real=True)
    with pytest.raises(ValueError, match="only 0 real roots"):
        evaluate_delineable_curve(
            ((z - y, 1),),
            (y**2 - x, 1),
            (x, y, z),
            (-1,),
        )


def test_delineable_surface_rejects_root_outside_requested_range():
    from semialg.cad_region import evaluate_delineable_surfaces

    x, y, z = sp.symbols("x y z", real=True)
    with pytest.raises(ValueError, match="outside requested range"):
        evaluate_delineable_surfaces(
            ((z - 2, 1),),
            (x, y, z),
            ((0, 0),),
            z_range=(-1, 1),
        )


def test_nonconforming_shared_face_is_detected():
    from semialg.cad_region import _mesh_is_conforming

    # Cell A uses one shared edge (0, 1); cell B splits that same sampled
    # interface at vertex 2.  The provenance marks 0 and 1 as shared while
    # vertex 2 belongs only to B, so the induced shared-face complexes differ.
    simplices = ((0, 1, 3), (0, 2, 4), (2, 1, 4))
    owners = ((1,), (2,), (2,))
    vertex_cells = (
        ((1,), (2,)),
        ((1,), (2,)),
        ((2,),),
        ((1,),),
        ((2,),),
    )
    assert not _mesh_is_conforming(simplices, owners, vertex_cells)
