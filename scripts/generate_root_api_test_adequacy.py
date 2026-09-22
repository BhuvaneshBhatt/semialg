"""Generate a risk-adjusted, function-by-function root API test adequacy audit."""

from __future__ import annotations

import ast
import inspect
import sys
import tomllib
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
import semialg  # noqa: E402

DIMENSION_KEYWORDS = {
    "boundary-degenerate": (
        "boundary",
        "degenerate",
        "endpoint",
        "singular",
        "empty",
        "zero",
        "lower_dim",
        "repeated",
        "rank_drop",
        "strict",
    ),
    "invalid-failure": (
        "invalid",
        "reject",
        "raises",
        "error",
        "malformed",
        "failure",
        "unsupported",
        "bad_",
    ),
    "metamorphic": (
        "invariant",
        "covariant",
        "translation",
        "scal",
        "rename",
        "permut",
        "idempot",
        "commut",
        "equiv",
        "consistent",
        "symmetr",
        "rotation",
    ),
    "independent-oracle": (
        "oracle",
        "brute",
        "direct",
        "closed_form",
        "known",
        "analytic",
        "exact_value",
        "ground_truth",
    ),
    "differential": (
        "differential",
        "agrees",
        "agreement",
        "compare",
        "same_as",
        "matches",
        "versus",
        "vs_",
    ),
    "parameter-regime": (
        "parameter",
        "regime",
        "conditional",
        "specializ",
        "piecewise",
        "free_parameter",
    ),
    "assumptions": ("assumption", "assumptions"),
    "certificate": ("certificate", "certif", "replay", "witness", "proof"),
    "round-trip": ("round_trip", "roundtrip", "reconstruct", "recovery", "inverse"),
    "representation": (
        "representation",
        "canonical",
        "expanded",
        "factored",
        "equivalent_form",
        "presentation",
    ),
    "negative": ("false", "unsat", "negative", "nonexist", "not_", "excludes", "outside"),
    "property-generated": ("generated", "hypothesis", "property", "random"),
}
INDEPENDENT = {
    "metamorphic",
    "independent-oracle",
    "independent-semantic",
    "differential",
    "certificate",
    "round-trip",
    "round-trip/presentation",
    "representation",
    "property-generated",
}
ADVERSE = {
    "boundary-degenerate",
    "invalid-failure",
    "negative",
    "parameter-regime",
    "assumptions",
    "certificate",
}


def collect():
    meta = tomllib.loads((ROOT / "tests/public_api_coverage.toml").read_text())["apis"]
    prior = tomllib.loads((ROOT / "tests/root_api_adequacy_dimensions.toml").read_text())["apis"]
    public = {n for n in semialg.__all__ if inspect.isfunction(getattr(semialg, n))}
    ev = defaultdict(lambda: defaultdict(set))
    for path in (ROOT / "tests").rglob("test_*.py"):
        tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))

        class V(ast.NodeVisitor):
            def __init__(self, current_path):
                self.current_path = current_path
                self.stack = []

            def visit_FunctionDef(self, node):
                self.stack.append(node.name)
                self.generic_visit(node)
                self.stack.pop()

            visit_AsyncFunctionDef = visit_FunctionDef

            def visit_Call(self, node):
                name = None
                if isinstance(node.func, ast.Name):
                    name = node.func.id
                elif (
                    isinstance(node.func, ast.Attribute)
                    and isinstance(node.func.value, ast.Name)
                    and node.func.value.id in {"semialg", "sa"}
                ):
                    name = node.func.attr
                if name in public and self.stack and self.stack[-1].startswith("test_"):
                    ref = f"{self.current_path.relative_to(ROOT).as_posix()}::{self.stack[-1]}"
                    ev[name]["nominal"].add(ref)
                    hay = (self.stack[-1] + " " + self.current_path.as_posix()).lower()
                    for dim, keys in DIMENSION_KEYWORDS.items():
                        if any(k in hay for k in keys):
                            ev[name][dim].add(ref)
                self.generic_visit(node)

        V(path).visit(tree)
    for name, item in prior.items():
        if name in public:
            for dim in item["dimensions"]:
                ev[name][dim].update(item.get("evidence", []))
    return public, meta, ev


def verdict(risk, dims):
    sub = dims - {"nominal"}
    independent = sub & INDEPENDENT
    adverse = sub & ADVERSE
    if risk == "critical":
        missing = []
        if len(dims) < 3:
            missing.append("at least three semantic dimensions")
        if not independent:
            missing.append("an independent mathematical detector")
        if not adverse:
            missing.append("an adverse/boundary/conditional contract")
    elif risk == "high":
        missing = []
        if len(dims) < 2:
            missing.append("a second semantic dimension")
        if not (
            independent
            or sub & {"boundary-degenerate", "invalid-failure", "parameter-regime", "assumptions"}
        ):
            missing.append("an independent, boundary, failure, or parameter-sensitive detector")
    else:
        missing = [] if len(dims) >= 2 else ["a second semantic dimension"]
    return ("adequate", []) if not missing else ("needs-deepening", missing)


def render():
    public, meta, ev = collect()
    rows = []
    for n in sorted(public):
        dims = set(ev[n])
        status, missing = verdict(meta[n]["risk"], dims)
        rows.append(
            (
                n,
                meta[n]["owner"],
                meta[n]["risk"],
                status,
                sorted(dims),
                missing,
                {d: sorted(v) for d, v in sorted(ev[n].items())},
            )
        )
    out = [
        "# Generated by scripts/generate_root_api_test_adequacy.py; do not edit manually.",
        "schema_version = 1",
        "",
    ]
    for n, owner, risk, status, dims, missing, evidence in rows:

        def q(s):
            return '"' + s.replace("\\", "\\\\").replace('"', '\\"') + '"'

        out += [
            f"[apis.{q(n)}]",
            f"owner = {q(owner)}",
            f"risk = {q(risk)}",
            f"status = {q(status)}",
            "dimensions = [" + ", ".join(map(q, dims)) + "]",
            "missing = [" + ", ".join(map(q, missing)) + "]",
        ]
        for d, refs in evidence.items():
            out += [f"[apis.{q(n)}.evidence.{q(d)}]", "tests = [" + ", ".join(map(q, refs)) + "]"]
        out.append("")
    report = [
        "# Root API test adequacy",
        "",
        "This generated audit evaluates every root-level public function by **semantic test dimensions**, not line coverage or raw test count. Existing explicitly curated multidimensional evidence is combined with directly-called tests whose names identify a semantic contract.",
        "",
        "## Risk-adjusted adequacy floor",
        "",
        "- **Critical:** nominal coverage, at least three semantic dimensions, at least one independent mathematical detector, and at least one adverse/boundary/conditional detector.",
        "- **High:** nominal coverage plus a meaningful independent, boundary, failure, assumptions, or parameter-sensitive detector.",
        "- **Medium:** nominal coverage plus at least one second semantic dimension.",
        "",
        "> This is an adequacy floor, not a proof of correctness. Test names only count when the test directly calls the root API; explicit multidimensional evidence is treated as authoritative.",
        "",
        "## Results",
        "",
    ]
    good = sum(r[3] == "adequate" for r in rows)
    report += [
        f"**{good}/{len(rows)} adequate; {len(rows) - good} need deepening.**",
        "",
        "| API | Owner | Risk | Status | Dimensions | Missing |",
        "|---|---|---|---|---|---|",
    ]
    for n, o, r, s, d, m, _ in rows:
        report.append(f"| `{n}` | {o} | {r} | {s} | {', '.join(d)} | {'; '.join(m) or '—'} |")
    return "\n".join(out) + "\n", "\n".join(report) + "\n"


if __name__ == "__main__":
    t, m = render()
    (ROOT / "tests/root_api_test_adequacy.toml").write_text(t)
    (ROOT / "docs/quality/root-api-test-adequacy.md").write_text(m)
