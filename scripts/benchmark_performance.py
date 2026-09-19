#!/usr/bin/env python3
"""Run semialg's dedicated performance-regression probes and emit JSON."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from semialg.benchmarks.performance import probes_as_dict


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    payload = {"schema": 1, "probes": probes_as_dict()}
    text = json.dumps(payload, indent=2, sort_keys=True)
    if args.output:
        args.output.write_text(text + "\n", encoding="utf-8")
    else:
        print(text)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
