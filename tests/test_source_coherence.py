import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
TEXT_ROOTS = ("src", "tests", "docs", "examples")
CAS_TERMS = (
    "mathe" + "matica",
    "wolf" + "ram",
    "power" + "expand",
    "full" + "simplify",
    "generate" + "conditions",
)


def _text_files():
    suffixes = {".py", ".md", ".rst", ".toml", ".yml", ".yaml", ".txt", ".cfg"}
    for root_name in TEXT_ROOTS:
        for path in (ROOT / root_name).rglob("*"):
            if path.is_file() and path.suffix.lower() in suffixes:
                yield path
    yield ROOT / "README.md"


def test_no_external_cas_vocabulary():
    matches = []
    for path in _text_files():
        text = path.read_text(encoding="utf-8").lower()
        for term in CAS_TERMS:
            if re.search(rf"\b{re.escape(term)}\b", text):
                matches.append((path.relative_to(ROOT), term))
    assert matches == []


def test_no_trailing_whitespace():
    matches = []
    for path in _text_files():
        for number, line in enumerate(path.read_text(encoding="utf-8").splitlines(), 1):
            if line != line.rstrip():
                matches.append((path.relative_to(ROOT), number))
    assert matches == []


def test_region_coercion_is_not_runtime_patched():
    source = (ROOT / "src/semialg/region_coercion.py").read_text(encoding="utf-8")
    assert "symbolic_regions.as_semialgebraic_region =" not in source
