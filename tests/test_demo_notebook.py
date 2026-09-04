from __future__ import annotations

from pathlib import Path

import matplotlib
import matplotlib.pyplot as plt
import nbformat

matplotlib.use("Agg", force=True)

ROOT = Path(__file__).resolve().parents[1]
NOTEBOOK = ROOT / "notebooks" / "semialg_demo.ipynb"


def test_demo_notebook_executes_all_code_cells() -> None:
    """The shipped demo notebook is an executable public-API contract."""

    notebook = nbformat.read(NOTEBOOK, as_version=4)
    namespace: dict[str, object] = {"__name__": "__main__"}

    for index, cell in enumerate(notebook.cells):
        if cell.cell_type != "code":
            continue
        filename = f"{NOTEBOOK}:cell-{index}"
        exec(compile(cell.source, filename, "exec"), namespace)
        plt.close("all")


def test_demo_notebook_covers_main_capabilities() -> None:
    """Keep the package's main capabilities visible in the demonstration notebook."""

    notebook = nbformat.read(NOTEBOOK, as_version=4)
    source = "\n".join(cell.source for cell in notebook.cells)

    required_fragments = (
        "return_result=True",
        "sp.real_root(x, 3)",
        "sp.sin(x) + sp.cos(2*x)",
        "sp.exp(x) + sp.exp(-x)",
        "simplify_semialgebraic_formula",
        'piece.diagnostics["integration_variable_order"]',
        "suggest_cad_variable_order",
    )
    for fragment in required_fragments:
        assert fragment in source, f"demo notebook is missing {fragment!r}"
