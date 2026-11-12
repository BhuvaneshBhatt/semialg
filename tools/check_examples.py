"""Execute the advertised examples and marked tutorials against installed packages."""

from __future__ import annotations

import argparse
import importlib
import importlib.util
import json
import os
import re
import subprocess
import sys
import tempfile
import tomllib
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
PACKAGE = tomllib.loads((ROOT / "pyproject.toml").read_text(encoding="utf-8"))["project"]["name"]


def catalog():
    """Validate that the catalog accounts for every standalone application script."""
    entries = json.loads((ROOT / "examples/catalog.json").read_text(encoding="utf-8"))
    paths = [row["file"] for row in entries]
    actual = {
        p.relative_to(ROOT / "examples").as_posix()
        for p in (ROOT / "examples").rglob("*.py")
        if p.name != "run_all.py"
    }
    if len(paths) != len(set(paths)) or set(paths) != actual:
        raise ValueError(
            f"example catalog mismatch: missing={actual - set(paths)}, stale={set(paths) - actual}"
        )
    if any(row.get("profile", "core") not in {"core", "real"} for row in entries):
        raise ValueError("example profiles must be core or real")
    return entries


def snippets(path):
    """Select executable Python fences while leaving displayed signatures as text."""
    text = path.read_text(encoding="utf-8")
    all_blocks = PACKAGE == "alggeom" or "<!-- examples: executable -->" in text
    pattern = r"(?P<marker><!--\s*(?:algroots: execute|semialg-exec|example: executable)\s*-->\s*)?^```python[^\S\n]*\n(?P<code>.*?)^```[^\S\n]*$"
    for match in re.finditer(pattern, text, re.M | re.S):
        if all_blocks or match.group("marker"):
            yield "\n" * text.count("\n", 0, match.start("code")) + match.group("code")


def installed():
    """Reject editable/source imports, including sibling backend checkouts."""
    for name in (
        PACKAGE,
        "algroots",
        *(("semialg",) if PACKAGE == "alggeom" and importlib.util.find_spec("semialg") else ()),
    ):
        module = importlib.import_module(name)
        origin = Path(module.__file__).resolve()
        if origin.is_relative_to(ROOT / "src") or any(
            origin.is_relative_to(ROOT.parent / p / "src")
            for p in ("alggeom", "algroots", "semialg")
        ):
            raise RuntimeError(
                f"{name} imported source instead of an installed distribution: {origin}"
            )


def execute(code, path, timeout):
    """Isolate each workflow from the checkout and from previous example state."""
    # Explicit interpreter paths preserve the installation under test without
    # evaluating editable-install hooks from the user's normal site directory.
    bootstrap = (
        "import sys; sys.path[:]= "
        + repr([str(Path(p).resolve()) for p in sys.path if p])
        + "; "
        + code
    )
    with tempfile.TemporaryDirectory(prefix=PACKAGE + "-example-") as temp:
        environment = {
            **os.environ,
            "MPLBACKEND": "Agg",
            "MPLCONFIGDIR": temp,
            "PYTHONDONTWRITEBYTECODE": "1",
        }
        result = subprocess.run(
            [sys.executable, "-S", "-B", "-c", bootstrap],
            cwd=temp,
            env=environment,
            capture_output=True,
            text=True,
            timeout=timeout,
            check=False,
        )
    if result.returncode:
        raise RuntimeError(f"{path}:\n{result.stdout}\n{result.stderr}")


def main():
    """Check the selected dependency profile with real mathematical assertions enabled."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--profile", choices=("core", "real"), default="real")
    parser.add_argument("--timeout", type=float, default=180)
    args = parser.parse_args()
    installed()
    entries = catalog()
    selected = [
        row for row in entries if args.profile == "real" or row.get("profile", "core") == "core"
    ]
    for row in selected:
        path = (ROOT / "examples" / row["file"]).resolve()
        execute(
            "import runpy; runpy.run_path(" + repr(str(path)) + ',run_name="__main__")',
            path,
            args.timeout,
        )
        print("Example:", row["file"], flush=True)
    blocks = 0
    for path in (ROOT / "README.md", *sorted((ROOT / "docs").rglob("*.md"))):
        text = path.read_text(encoding="utf-8")
        if args.profile == "core" and (
            ("<!-- examples: real -->" in text)
            or (
                PACKAGE == "alggeom"
                and path.relative_to(ROOT).as_posix() in {"docs/applications.md"}
            )
        ):
            continue
        source = list(snippets(path))
        if source:
            execute(
                "exec(compile(" + repr("\n".join(source)) + "," + repr(str(path)) + ',"exec"),{})',
                path,
                args.timeout,
            )
            blocks += len(source)
    print(
        f"{PACKAGE}: {len(selected)} installed examples and {blocks} executable documentation blocks passed",
        flush=True,
    )


if __name__ == "__main__":
    main()
