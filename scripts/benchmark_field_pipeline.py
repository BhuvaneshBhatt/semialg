"""Compare verified algebra and real-solvability pipelines across Gröbner engines."""

import argparse
import hashlib
import json
import os
import platform
import random
import statistics
import sys
from contextlib import contextmanager
from functools import partial
from importlib.metadata import version
from pathlib import Path
from time import perf_counter, process_time

import sympy as sp
from algroots import GroebnerBackendUnavailable, NativeF4Limits
from algroots.profiling import GroebnerProfiler
from sympy.core.random import seed as seed_sympy
from sympy.polys.groebnertools import _buchberger, _f5b
from sympy.polys.polyconfig import query, setup

from semialg import solvability_conditions
from semialg.algebraic.gtz_primary import (
    _compute_fraction_field_zero_dimensional,
    gtz_primary_decomposition,
    verify_fraction_field_zero_dim_certificate,
    verify_gtz_primary_decomposition_certificate,
)
from semialg.algebraic.parameter_groebner import guarded_equality_preprocessing
from semialg.algebraic.rational_univariate import compute_rur
from semialg.cache_control import clear_caches


@contextmanager
def sympy_method(method):
    """Select the process-wide SymPy method for a serial run and restore it on errors."""
    previous = query("groebner")
    setup("groebner", method)
    try:
        yield
    finally:
        setup("groebner", previous)


@contextmanager
def audit_sympy_methods():
    """Count actual algorithm entries in a separate, untimed operation."""
    counts = {"buchberger": 0, "f5b": 0}
    algorithms = {_buchberger.__code__: "buchberger", _f5b.__code__: "f5b"}
    previous = sys.getprofile()

    def observe(frame, event, arg):
        if event == "call" and frame.f_code in algorithms:
            counts[algorithms[frame.f_code]] += 1
        if previous is not None:
            previous(frame, event, arg)

    sys.setprofile(observe)
    try:
        yield counts
    finally:
        sys.setprofile(previous)


def measure(operation, repetitions):
    """Warm once, then time full operations with cache teardown outside each sample."""
    clear_caches(collect=False)
    operation()
    samples, cpu_samples, stages, profiles = [], [], [], []
    for _ in range(repetitions):
        clear_caches(collect=False)
        started, cpu_started = perf_counter(), process_time()
        result, profile, stage_times = operation()
        cpu_samples.append(process_time() - cpu_started)
        samples.append(perf_counter() - started)
        stages.append(stage_times)
        profiles.append(profile)
    return (
        result,
        {
            "samples_seconds": samples,
            "median_seconds": statistics.median(samples),
            "cpu_samples_seconds": cpu_samples,
            "stage_samples_seconds": stages,
        },
        profiles,
    )


def rur_operation(source, variables, domain, backend, limits):
    with GroebnerProfiler() as profile:
        result = compute_rur(
            source, variables, domain=domain, groebner_backend=backend, groebner_limits=limits
        )
        result.normalized_coordinate_polynomials()
    return result, profile.records, {}


def decomposition_operation(source, variables, parameters, backend, limits):
    with GroebnerProfiler() as profile:
        if parameters:
            result = _compute_fraction_field_zero_dimensional(
                source, variables, parameters, groebner_backend=backend, groebner_limits=limits
            )
        else:
            result = gtz_primary_decomposition(
                source, variables, groebner_backend=backend, groebner_limits=limits
            )
    return result, profile.records, {}


def solvability_operation(source, variables, parameters, backend, limits, *, matched):
    """Compare the public dispatch and an identical guarded-preprocessing pipeline.

    The public SymPy selection omits optional equality preprocessing. The
    matched variant explicitly preprocesses with either engine, then calls
    the same public real eliminator with preprocessing disabled. This isolates
    the preprocessing policy while retaining the full exceptional real branches.
    """
    stages = {}
    with GroebnerProfiler() as profile:
        if matched:
            started = perf_counter()
            formula, diagnostics = guarded_equality_preprocessing(
                sp.And(*source), variables, parameters, backend=backend, limits=limits
            )
            stages["equality_preprocessing"] = perf_counter() - started
            preprocessing_calls = len(profile.records)
            started = perf_counter()
            result = solvability_conditions(
                formula, variables, parameters, groebner_backend="sympy", return_result=True
            )
            stages["real_elimination"] = perf_counter() - started
            result.diagnostics.update(diagnostics)
            for index, record in enumerate(profile.records):
                record["phase"] = (
                    "equality_preprocessing" if index < preprocessing_calls else "real_elimination"
                )
        else:
            result = solvability_conditions(
                source,
                variables,
                parameters,
                groebner_backend=backend,
                groebner_limits=limits,
                return_result=True,
            )
    return result, profile.records, stages


def verify_rur(source, variables, domain, length, count, result):
    """Check coordinate identities and prescribed algebra lengths and point counts."""
    if result.dimension != length or result.solution_count != count:
        raise ArithmeticError("RUR quotient length or geometric count is wrong")
    coordinates = result.normalized_coordinate_polynomials()
    maps = dict(zip(variables, (p.as_expr() for p in coordinates), strict=True))
    for polynomial in source:
        value = sp.cancel(polynomial.subs(maps, simultaneous=True))
        if (
            not sp.Poly(value, result.parameter, domain=domain)
            .rem(result.defining_polynomial)
            .is_zero
        ):
            raise ArithmeticError("RUR coordinate identity failed")
    return {
        "dimension": result.dimension,
        "geometric_points": result.solution_count,
        "defining_coefficients": list(map(str, result.defining_polynomial.all_coeffs())),
        "coordinates": list(map(str, coordinates)),
    }


def verify_decomposition(parameters, count, result):
    """Replay every decomposition certificate outside the timed calculation."""
    verifier = (
        verify_fraction_field_zero_dim_certificate
        if parameters
        else verify_gtz_primary_decomposition_certificate
    )
    if not result.complete or not verifier(result.certificate) or len(result.components) != count:
        raise ArithmeticError("decomposition certificate or component count failed")
    return {"components": len(result.components), "complete": True, "certificate_replayed": True}


def real_solvability_cases():
    """Return small real problems with independently derived parameter conditions."""
    x, y, t = sp.symbols("x y t", real=True)
    cases = [
        ("exceptional_linear", (sp.Eq(t * x, 0), x > 1), sp.Eq(t, 0)),
        ("positive_reciprocal", (sp.Eq(t * x, 1), x > 0), t > 0),
        ("exceptional_quadratic", (sp.Eq(t * (x**2 - 1), 0), x > 2), sp.Eq(t, 0)),
        ("bounded_square", (sp.Eq(x**2, t), x > 0, x < 1), sp.And(t > 0, t < 1)),
        (
            "two_exceptional_fibers",
            (sp.Eq(t * (t - 1) * x, 0), x > 1),
            sp.Or(sp.Eq(t, 0), sp.Eq(t, 1)),
        ),
        ("coefficient_pole", (sp.Eq(x / t, 0), x >= 0), sp.Ne(t, 0)),
    ]

    return [(name, source, (x,), expected) for name, source, expected in cases] + [
        (
            "coupled_bounds",
            (sp.Eq(t * x - y, 0), sp.Eq(x**2, 1), y > 1),
            (x, y),
            sp.Or(t < -1, t > 1),
        )
    ]


def verify_solvability(parameter, expected, result):
    """Prove whole-line equivalence using exact univariate inequality reduction.

    Sampling alone would miss an exceptional fiber or a strict endpoint.
    Every disjunct of the symmetric difference must be an empty real set.
    """
    difference = sp.to_dnf(sp.Xor(result.formula, expected), simplify=True)
    for branch in sp.Or.make_args(difference):
        if sp.reduce_inequalities(sp.And.make_args(branch), parameter) is not sp.false:
            raise ArithmeticError(
                f"wrong parameter condition: {result.formula}; expected {expected}"
            )
    return {
        "formula": str(result.formula),
        "expected_formula": str(expected),
        "equivalent_on_all_real_parameters": True,
        "verification_method": "SymPy univariate inequality reduction of symmetric difference",
        "diagnostics": {key: str(value) for key, value in result.diagnostics.items()},
    }


def workloads(backend, limits):
    """Build algebraic workloads and both public and matched real-decision variants."""
    x, y, t, u = sp.symbols("x y t u")
    K = sp.QQ.algebraic_field(sp.sqrt(2))
    rurs = [
        ("rur_rational", (x * x - y, y * y - 1), sp.QQ, 4, 4),
        ("rur_number_repeated", ((x * x - sp.sqrt(2)) ** 2, y - x), K, 4, 2),
        ("rur_parameter_repeated", ((x - t) ** 2 * (x + t) ** 3, y - x), sp.QQ.frac_field(t), 5, 2),
        ("rur_mixed_repeated", ((x * x - t - sp.sqrt(2)) ** 2, y - x), K.frac_field(t), 4, 2),
        ("rur_parameter_coupled", (x * y - t, y * y - x - 1), sp.QQ.frac_field(t), 3, 3),
        ("rur_mixed_coupled", (x * y - t - sp.sqrt(2), y * y - x - 1), K.frac_field(t), 3, 3),
    ]
    cases = {}
    for name, source, domain, length, count in rurs:
        cases[name] = dict(
            kind="rur",
            source=list(map(str, source)),
            variables=[str(x), str(y)],
            domain=str(domain),
            operation=partial(rur_operation, source, (x, y), domain, backend, limits),
            verify=partial(verify_rur, source, (x, y), domain, length, count),
        )
    decompositions = [
        ("gtz_embedded", ((x + y) ** 2, y * (x + y)), (x, y), (), 2),
        ("gtz_residual", (x * y,), (x, y), (), 2),
        ("localized_repeated", ((x - t) ** 2 * (x + t) ** 3, y * y), (x, y), (t,), 2),
        ("localized_two_parameters", (x * x - t, y * y - u), (x, y), (t, u), 1),
    ]
    for name, source, variables, parameters, count in decompositions:
        cases[name] = dict(
            kind="primary_decomposition",
            source=list(map(str, source)),
            variables=list(map(str, variables)),
            parameters=list(map(str, parameters)),
            operation=partial(
                decomposition_operation, source, variables, parameters, backend, limits
            ),
            verify=partial(verify_decomposition, parameters, count),
        )
    x, t = sp.symbols("x t", real=True)
    for name, source, variables, expected in real_solvability_cases():
        for mode in ("api", "matched"):
            cases[f"solvability_{mode}_{name}"] = dict(
                kind="solvability",
                comparison=mode,
                source=list(map(str, source)),
                variables=list(map(str, variables)),
                parameters=[str(t)],
                expected_formula=str(expected),
                operation=partial(
                    solvability_operation,
                    source,
                    variables,
                    (t,),
                    backend,
                    limits,
                    matched=mode == "matched",
                ),
                verify=partial(verify_solvability, t, expected),
            )
    return cases


def run_case(name, case, backend, method, limits, repetitions):
    """Retain per-sample profiles and an untimed audit of the selected SymPy engine."""
    selected_method = method if backend == "sympy" else "buchberger"
    metadata = {key: value for key, value in case.items() if key not in {"operation", "verify"}}
    with sympy_method(selected_method):
        result, timing, profiles = measure(case["operation"], repetitions)
        started = perf_counter()
        check = case["verify"](result)
        replay_seconds = perf_counter() - started
        clear_caches(collect=False)
        with audit_sympy_methods() as audit:
            audited_result, audit_profile, _ = case["operation"]()
        case["verify"](audited_result)
    wrong_method = "buchberger" if selected_method == "f5b" else "f5b"
    if audit[wrong_method]:
        raise ArithmeticError(f"{name}: unexpected SymPy {wrong_method} calls: {audit}")
    return dict(
        case=name,
        **metadata,
        backend=backend,
        construction_method=method if backend == "sympy" else "f4_tiled_finish",
        sympy_method=selected_method,
        native_reference_sympy_method="buchberger" if backend == "native" else None,
        native_limits=repr(limits) if limits is not None else None,
        timing=timing,
        independent_check_seconds=replay_seconds,
        check=check,
        groebner_calls=profiles[-1],
        groebner_profiles=profiles,
        method_audit={
            "timed": False,
            "sympy_kernel_calls": audit,
            "profiled_backends": sorted({row["backend"] for row in audit_profile}),
        },
        status="verified",
    )


def main():
    """Write verified timings, explicit methods, cache policy and independent checks."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--repetitions", type=int, default=3)
    parser.add_argument(
        "--backends", nargs="+", choices=["sympy", "native"], default=["sympy", "native"]
    )
    parser.add_argument(
        "--sympy-method",
        choices=["buchberger", "f5b"],
        default="buchberger",
        help="SymPy construction method; native reference checks stay Buchberger",
    )
    parser.add_argument("--cases", nargs="+", help="Exact workload names; omit to run all")
    parser.add_argument("--seed", type=int, default=1729)
    args = parser.parse_args()
    random.seed(args.seed)
    seed_sympy(args.seed)
    if args.repetitions < 1:
        parser.error("repetitions must be positive")
    limits = NativeF4Limits(collect_batch_metrics=True, max_reduction_steps=2000000)
    report = dict(
        environment=dict(
            python=platform.python_version(),
            platform=platform.platform(),
            packages={
                name: version(name) for name in ("semialg", "algroots", "sympy", "python-flint")
            },
        ),
        sympy_method=args.sympy_method,
        seed=args.seed,
        pythonhashseed=os.environ.get("PYTHONHASHSEED"),
        driver_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        native_reference_sympy_method="buchberger",
        policy={
            "timing": "serial warm measurements, with all semialg caches cleared outside every timed operation; SymPy expression caches stay warm",
            "verification": "all algorithm-internal checks included; independent result checks and method audit are outside timed operations",
            "profiles": "per-sample public algroots calls only; direct SymPy calls appear in the separate method audit; nested timers must not be summed as disjoint stages",
            "api_comparison": "public SymPy dispatch skips optional equality preprocessing; public native dispatch enables it",
            "matched_comparison": "both engines preprocess equalities, then use the same public solvability_conditions with preprocessing disabled; includes both stages",
        },
        records=[],
        unavailable=[],
    )
    for backend in args.backends:
        selected_limits = limits if backend == "native" else None
        cases = workloads(backend, selected_limits)
        if args.cases and set(args.cases) - cases.keys():
            parser.error(f"unknown cases: {sorted(set(args.cases) - cases.keys())}")
        for name in args.cases or cases:
            print(
                f"Running {backend}/{args.sympy_method if backend == 'sympy' else 'F4'}: {name}",
                flush=True,
            )
            try:
                record = run_case(
                    name, cases[name], backend, args.sympy_method, selected_limits, args.repetitions
                )
            except GroebnerBackendUnavailable as error:
                report["unavailable"].append(dict(case=name, backend=backend, reason=str(error)))
            else:
                report["records"].append(record)
            args.output.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(f"Verified {len(report['records'])} complete pipeline measurements", flush=True)


if __name__ == "__main__":
    main()
