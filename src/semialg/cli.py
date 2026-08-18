from __future__ import annotations

import argparse
import json
from pathlib import Path

import sympy as sp

from .decomposition import cad_text, component_instances_text, generic_cad_text
from .solve.find_instance import find_instance_text
from .solve.reduce import reduce_text
from .solve.resolve import resolve_text
from .validation.checkers import SymPyInequalityChecker
from .validation.corpus import built_in_smoke_cases, read_jsonl_cases
from .validation.regression import write_failing_cases
from .validation.runner import run_validation_cases


def _symbol_table(names: list[str] | None) -> dict[str, sp.Symbol]:
    return {name: sp.Symbol(name, real=True) for name in (names or [])}


def _json_default(value):
    return str(value)


def run_reduce(args: argparse.Namespace) -> None:
    symbols = _symbol_table(args.variables)
    result = reduce_text(
        args.formula,
        symbols=symbols or None,
        variable_order=tuple(symbols.values()) or None,
        domain=args.domain,
        return_result=True,
        strategy=args.strategy,
    )
    if args.json:
        payload = {
            "method": result.method,
            "domain": str(result.domain),
            "result": str(result.result),
            "metadata": {key: str(value) for key, value in result.metadata.items()},
        }
        print(json.dumps(payload, indent=2, sort_keys=True))
    else:
        print(result.result)


def run_resolve(args: argparse.Namespace) -> None:
    symbols = _symbol_table(args.variables)
    result = resolve_text(
        args.formula,
        symbols=symbols or None,
        variable_order=tuple(symbols.values()) or None,
        domain=args.domain,
        strategy=args.strategy,
    )
    print(
        json.dumps(result, indent=2, sort_keys=True, default=_json_default) if args.json else result
    )


def run_cad(args: argparse.Namespace) -> None:
    symbols = _symbol_table(args.variables)
    result = cad_text(
        args.formula,
        variables=tuple(symbols.values()) or None,
        symbols=symbols or None,
        output=args.output,
        operation=args.operation,
        domain=args.domain,
        strategy=args.strategy,
        return_result=True,
    )
    if args.json:
        payload = {
            "status": result.status,
            "formula": str(result.formula),
            "variables": [str(var) for var in result.variables],
            "cell_count_by_level": {str(k): v for k, v in result.cell_count_by_level().items()},
            "projection_polynomial_count_by_level": {
                str(k): v for k, v in result.proj_poly_count_by_level().items()
            },
        }
        print(json.dumps(payload, indent=2, sort_keys=True, default=_json_default))
    else:
        print(result.formula)


def run_generic_cad(args: argparse.Namespace) -> None:
    names = list(dict.fromkeys((args.parameters or []) + (args.variables or [])))
    symbols = _symbol_table(names)
    result = generic_cad_text(
        args.formula,
        variables=[symbols[name] for name in (args.variables or [])],
        parameters=[symbols[name] for name in (args.parameters or [])],
        symbols=symbols or None,
        output=args.output,
        domain=args.domain,
        strategy=args.strategy,
        return_result=True,
    )
    if args.json:
        payload = {
            "status": result.status,
            "generic_formula": str(result.generic_formula),
            "exceptional_formula": str(result.exceptional_formula),
            "generic_case_count": len(result.generic_cases),
            "exceptional_case_count": len(result.exceptional_cases),
        }
        print(json.dumps(payload, indent=2, sort_keys=True, default=_json_default))
    else:
        print("generic:", result.generic_formula)
        print("exceptional:", result.exceptional_formula)


def run_components(args: argparse.Namespace) -> None:
    symbols = _symbol_table(args.variables)
    result = component_instances_text(
        args.formula,
        variables=tuple(symbols.values()) or None,
        symbols=symbols or None,
        domain=args.domain,
        strategy=args.strategy,
        max_components=args.count,
        return_result=True,
    )
    if args.json:
        payload = {
            "status": result.status,
            "instances": [{str(k): str(v) for k, v in inst.items()} for inst in result.instances],
            "approx_instances": [
                {str(k): v for k, v in inst.items()} for inst in result.approx_instances
            ],
            "component_count": len(result),
        }
        print(json.dumps(payload, indent=2, sort_keys=True, default=_json_default))
    else:
        for inst in result.instances:
            print(inst)


def run_instance(args: argparse.Namespace) -> None:
    symbols = _symbol_table(args.variables)
    result = find_instance_text(
        args.formula,
        symbols=symbols or None,
        variables=tuple(symbols.values()) or None,
        domain=args.domain,
        count=args.count,
        strategy=args.strategy,
        random_seed=args.random_seed,
        return_result=True,
    )
    if args.json:
        payload = {
            "status": result.status,
            "method": result.method,
            "domain": result.domain.value,
            "instances": [{str(k): str(v) for k, v in inst.items()} for inst in result.instances],
            "approximate": [{str(k): v for k, v in inst.items()} for inst in result.approximate],
            "diagnostics": {key: str(value) for key, value in result.diagnostics.items()},
        }
        print(json.dumps(payload, indent=2, sort_keys=True, default=_json_default))
    else:
        for inst in result.instances:
            print(inst)


def run_validate(args: argparse.Namespace) -> None:
    cases = built_in_smoke_cases() if args.corpus is None else read_jsonl_cases(args.corpus)
    checkers = [SymPyInequalityChecker()]
    report = run_validation_cases(cases, checkers=tuple(checkers))
    if args.report:
        report.write_json(args.report)
    if args.regressions:
        write_failing_cases(args.regressions, report.results)
    if args.json:
        print(report.to_json())
    else:
        passed = sum(1 for result in report.results if result.passed)
        total = len(report.results)
        print(f"validation: {passed}/{total} cases passed")
        for result in report.results:
            status = "PASS" if result.passed else "FAIL"
            print(f"{status} {result.case.name}: {result.solver_status}")


def _add_common_formula_args(cmd: argparse.ArgumentParser) -> None:
    cmd.add_argument("formula")
    cmd.add_argument("--variables", nargs="*", default=[])
    cmd.add_argument("--domain", default="reals")
    cmd.add_argument("--strategy", default="auto")
    cmd.add_argument("--json", action="store_true")


def main() -> None:
    parser = argparse.ArgumentParser(description="Semialgebraic CAD/QE tools")
    subparsers = parser.add_subparsers(dest="command")

    reduce_cmd = subparsers.add_parser("reduce", help="Eliminate quantifiers from a formula")
    _add_common_formula_args(reduce_cmd)
    reduce_cmd.set_defaults(func=run_reduce)

    resolve_cmd = subparsers.add_parser("resolve", help="Resolve a closed quantified formula")
    _add_common_formula_args(resolve_cmd)
    resolve_cmd.set_defaults(func=run_resolve)

    cad_cmd = subparsers.add_parser("cad", help="Compute a cylindrical algebraic decomposition")
    _add_common_formula_args(cad_cmd)
    cad_cmd.add_argument(
        "--output", choices=("formula", "cells", "function", "tree"), default="formula"
    )
    cad_cmd.add_argument("--operation", choices=("closure", "interior", "boundary"))
    cad_cmd.set_defaults(func=run_cad)

    generic_cmd = subparsers.add_parser("generic-cad", help="Compute a generic parameter-space CAD")
    _add_common_formula_args(generic_cmd)
    generic_cmd.add_argument("--parameters", nargs="*", default=[])
    generic_cmd.add_argument(
        "--output", choices=("formula", "cases", "cells", "function"), default="formula"
    )
    generic_cmd.set_defaults(func=run_generic_cad)

    components_cmd = subparsers.add_parser("components", help="Return one sample per component")
    _add_common_formula_args(components_cmd)
    components_cmd.add_argument("--count", type=int)
    components_cmd.set_defaults(func=run_components)

    instance_cmd = subparsers.add_parser("instance", help="Find satisfying assignments")
    _add_common_formula_args(instance_cmd)
    instance_cmd.add_argument("--count", type=int, default=1)
    instance_cmd.add_argument("--random-seed", type=int)
    instance_cmd.set_defaults(func=run_instance)

    validate_cmd = subparsers.add_parser("validate", help="Run validation corpus cases")
    validate_cmd.add_argument("--corpus", type=Path)
    validate_cmd.add_argument("--report", type=Path)
    validate_cmd.add_argument("--regressions", type=Path)
    validate_cmd.add_argument("--json", action="store_true")
    validate_cmd.set_defaults(func=run_validate)

    args = parser.parse_args()
    if not hasattr(args, "func"):
        parser.print_help()
        return
    args.func(args)


if __name__ == "__main__":
    main()
