"""Compare source, installed package, wheel, sdist and a wheel rebuilt from the sdist."""

from __future__ import annotations

import argparse
import importlib.util
import json
import subprocess
import sys
import tarfile
import tempfile
import zipfile
from pathlib import Path


def _modules(root):
    return {p.relative_to(root).as_posix(): p.read_bytes() for p in root.rglob("*.py")}


def _contract_files(root):
    extensions = {
        "docs": (".md", ".toml", ".svg", ".json"),
        "examples": (".py", ".md", ".json"),
        "notebooks": (".ipynb", ".md"),
        "scripts": (".py", ".sh"),
        "tests": (".py", ".toml", ".json", ".md"),
        "tools": (".py", ".md"),
        ".github": (".yml", ".yaml"),
        "ci": (".toml", ".txt", ".tar.gz", ".patch"),
    }
    paths = [
        root / name
        for name in (
            "README.md",
            "RELEASE_NOTES.md",
            "LICENSE",
            "pyproject.toml",
            "mkdocs.yml",
            "MANIFEST.in",
        )
    ]
    for directory, suffixes in extensions.items():
        paths.extend(
            p for p in (root / directory).rglob("*") if p.is_file() and p.name.endswith(suffixes)
        )
    return {p.relative_to(root).as_posix(): p.read_bytes() for p in paths}


def _wheel_modules(path):
    with zipfile.ZipFile(path) as archive:
        return {
            n.removeprefix("semialg/"): archive.read(n)
            for n in archive.namelist()
            if n.startswith("semialg/") and n.endswith(".py")
        }


def _require_match(expected, actual, label):
    missing = sorted(expected.keys() - actual.keys())
    extra = sorted(actual.keys() - expected.keys())
    changed = sorted(k for k in expected.keys() & actual.keys() if expected[k] != actual[k])
    if missing or extra or changed:
        raise ValueError(f"{label}: missing={missing}, extra={extra}, changed={changed}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("wheel", type=Path)
    parser.add_argument("sdist", type=Path)
    parser.add_argument("--output", type=Path, default=Path("distribution-check.json"))
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[1]
    expected = _modules(root / "src/semialg")
    contracts = _contract_files(root)
    _require_match(expected, _wheel_modules(args.wheel), "wheel")
    installed = Path(importlib.util.find_spec("semialg").origin).parent
    if installed == root / "src/semialg":
        raise ValueError("install the wheel before checking distribution contents")
    _require_match(expected, _modules(installed), "installed package")
    with tempfile.TemporaryDirectory(prefix="semialg-sdist-") as temp:
        folder = Path(temp)
        with tarfile.open(args.sdist) as archive:
            members = {
                m.name.split("/src/semialg/", 1)[1]: archive.extractfile(m).read()
                for m in archive.getmembers()
                if "/src/semialg/" in m.name and m.name.endswith(".py")
            }
            _require_match(expected, members, "sdist")
            payload = {
                m.name.partition("/")[2]: archive.extractfile(m).read()
                for m in archive.getmembers()
                if m.isfile()
            }
            _require_match(
                contracts,
                {name: payload[name] for name in contracts if name in payload},
                "sdist contracts",
            )
            archive.extractall(folder, filter="data")
        (project,) = [p for p in folder.iterdir() if p.is_dir()]
        subprocess.run(
            [
                sys.executable,
                "-m",
                "build",
                "--wheel",
                "--no-isolation",
                "--outdir",
                str(folder / "rebuilt"),
            ],
            cwd=project,
            check=True,
        )
        (rebuilt,) = (folder / "rebuilt").glob("*.whl")
        _require_match(expected, _wheel_modules(rebuilt), "rebuilt wheel")
    args.output.write_text(
        json.dumps(
            {
                "modules": len(expected),
                "contract_files": len(contracts),
                "verified": True,
                "formats": ["source", "wheel", "installed", "sdist", "rebuilt wheel"],
            },
            indent=2,
        )
        + "\n"
    )
    print(
        f"{len(expected)} package modules match in all five forms; {len(contracts)} sdist contract files match"
    )


if __name__ == "__main__":
    main()
