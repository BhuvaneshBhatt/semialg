from __future__ import annotations

import copy
import os
import subprocess
import sys
from pathlib import Path

import matplotlib
import nbformat
import pytest
from nbclient import NotebookClient

matplotlib.use("Agg", force=True)
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
NOTEBOOK = ROOT / "notebooks" / "semialg_demo.ipynb"


def test_demo_notebook_executes_in_a_real_kernel() -> None:
    """The shipped demo must execute with notebook/kernel semantics."""

    socket_probe = subprocess.run(
        [
            sys.executable,
            "-c",
            "import zmq; c=zmq.Context(); s=c.socket(zmq.ROUTER); "
            "s.bind_to_random_port('tcp://127.0.0.1'); s.close(); c.term()",
        ],
        capture_output=True,
        check=False,
    )
    if socket_probe.returncode != 0:
        pytest.skip("this sandbox forbids the local sockets required by Jupyter")

    notebook = nbformat.read(NOTEBOOK, as_version=4)
    kernel_env = os.environ.copy()
    source_path = str(ROOT / "src")
    kernel_env["PYTHONPATH"] = os.pathsep.join(
        part for part in (source_path, kernel_env.get("PYTHONPATH", "")) if part
    )
    executed = NotebookClient(
        copy.deepcopy(notebook),
        timeout=180,
        kernel_name="python3",
        resources={"metadata": {"path": str(ROOT)}},
    ).execute(env=kernel_env)

    code_cells = [cell for cell in executed.cells if cell.cell_type == "code"]
    assert [cell.execution_count for cell in code_cells] == list(range(1, len(code_cells) + 1))
    assert not [
        output
        for cell in code_cells
        for output in cell.get("outputs", ())
        if output.output_type == "error"
    ]
    plt.close("all")


def test_demo_notebook_code_is_valid_as_one_progressive_session() -> None:
    """Retain executable coverage in environments that prohibit Jupyter sockets."""

    notebook = nbformat.read(NOTEBOOK, as_version=4)
    namespace: dict[str, object] = {"__name__": "__main__"}
    for index, cell in enumerate(notebook.cells):
        if cell.cell_type == "code":
            filename = f"{NOTEBOOK}:cell-{index}"
            exec(compile(cell.source, filename, "exec"), namespace)
            plt.close("all")


def test_demo_notebook_has_stable_capability_tags() -> None:
    """Track tutorial coverage without coupling the test to exact source spelling."""

    notebook = nbformat.read(NOTEBOOK, as_version=4)
    tags = {tag for cell in notebook.cells for tag in cell.metadata.get("tags", ())}
    required = {
        "decision",
        "parameters",
        "roots",
        "integration",
        "certification",
        "performance",
    }
    assert required <= tags, f"demo notebook is missing capability tags: {required - tags}"
