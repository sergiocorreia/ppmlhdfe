#!/usr/bin/env python3
"""Run the maintained Schur benchmark cells in suite.csv."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys

import run as engine

HERE = Path(__file__).resolve().parent
FIELDS = {
    "platform",
    "design",
    "observations",
    "controls",
    "processors",
    "first_source",
    "annual_mobility",
    "fit_timeout",
    "budget",
}


def load_cells(path=HERE / "suite.csv"):
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        if set(reader.fieldnames or ()) != FIELDS:
            raise ValueError("Unexpected suite.csv columns")
        cells = list(reader)
    if not cells:
        raise ValueError("Benchmark suite is empty")
    seen = set()
    for cell in cells:
        for key in ("observations", "controls", "processors", "fit_timeout", "budget"):
            cell[key] = int(cell[key])
        key = (cell["platform"], cell["design"], cell["controls"], cell["processors"])
        if key in seen:
            raise ValueError("Duplicate benchmark cell: " + str(key))
        seen.add(key)
        if cell["platform"] not in ("scc", "mac") or cell["design"] not in engine.CASES:
            raise ValueError("Unknown platform or design")
        allowed = (1, 4) if cell["platform"] == "scc" else (4, 8)
        if cell["processors"] not in allowed or cell["controls"] not in (4, 16, 64):
            raise ValueError("Unsupported processor or control count")
        if cell["first_source"] not in ("original", "candidate"):
            raise ValueError("Invalid source order")
        if min(cell["fit_timeout"], cell["budget"]) <= 0:
            raise ValueError("Time limits must be positive")
        rate = {"akm_low": ".02", "akm_medium": ".10", "akm_high": ".30"}.get(
            cell["design"], ""
        )
        if cell["annual_mobility"] != rate or cell["observations"] != (
            200000 if rate else 100000
        ):
            raise ValueError("Fixture size or mobility differs from design")
    # Processor pairs run consecutively, independent of CSV row ordering.
    order = {name: i for i, name in enumerate(engine.CASES)}
    return sorted(
        cells,
        key=lambda c: (
            c["platform"],
            order[c["design"]],
            c["controls"],
            c["processors"],
        ),
    )


def cell_directory(root, cell):
    return root / cell["design"] / f"k{cell['controls']}_p{cell['processors']}"


def reuse_option(parent):
    """A failed fit can still supply an intact fixture, never a reference fit."""
    if parent is None:
        return []
    summary = parent / "summary.json"
    try:
        if summary.exists() and json.loads(summary.read_text()).get("complete"):
            return ["--reference-from", str(parent)]
    except (OSError, ValueError, AttributeError):
        pass  # An unreadable fit result can still have an independently valid fixture.
    manifest = parent / "manifest.json"
    if manifest.exists():
        prior = json.loads(manifest.read_text())
        if prior.get("fixture_sha256"):
            return ["--fixture-from", str(parent)]
    return []


def command(args, cell, directory, reuse):
    result = [
        sys.executable,
        str(HERE / "run.py"),
        "--original",
        str(args.original),
        "--candidate",
        str(args.candidate),
        "--git-repo",
        str(args.git_repo),
        "--candidate-revision",
        args.candidate_revision,
        "--dependency-dir",
        str(args.dependency_dir),
        "--output",
        str(directory),
        "--cases",
        cell["design"],
        "--n",
        str(cell["observations"]),
        "--controls",
        str(cell["controls"]),
        "--processors",
        str(cell["processors"]),
        "--precision",
        "tight",
        "--repetitions",
        "3",
        "--first-source",
        cell["first_source"],
        "--fit-timeout",
        str(cell["fit_timeout"]),
        "--budget",
        str(cell["budget"]),
    ]
    if args.stata:
        result += ["--stata", str(args.stata)]
    return result + reuse


def save_report(root, cells):
    result = {"complete": all(c["status"] == "passed" for c in cells), "cells": cells}
    (root / "summary.json").write_text(json.dumps(result, indent=2) + "\n")
    lines = [
        "# Schur benchmark results",
        "",
        "Seconds are medians of three fits after a warmup. Both sources use the same tight settings.",
        "",
        "| Design | K | P | Original, s | Schur, s | Speedup | Status |",
        "|---|---:|---:|---:|---:|---:|---|",
    ]
    for c in cells:
        m = c.get("measurements")
        values = (
            [
                f"{m['original']['median']:.3f}",
                f"{m['candidate']['median']:.3f}",
                f"{m['speedup']:.3f}×",
            ]
            if c["status"] == "passed"
            else ["—", "—", "—"]
        )
        lines.append(
            "| "
            + " | ".join(
                map(
                    str,
                    [c["design"], c["controls"], c["processors"], *values, c["status"]],
                )
            )
            + " |"
        )
    lines += [
        "",
        "Per-cell folders retain all repetitions, ranges, comparisons and failures.",
        "Only complete comparisons receive a speedup. Platforms are reported separately.",
    ]
    (root / "REPORT.md").write_text("\n".join(lines) + "\n")
    return result["complete"]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--platform", required=True, choices=("scc", "mac"))
    parser.add_argument(
        "--cases", help="Optional comma-separated subset of suite designs"
    )
    for name in (
        "original",
        "candidate",
        "git-repo",
        "dependency-dir",
        "stata",
        "output",
    ):
        parser.add_argument("--" + name, type=Path)
    parser.add_argument("--candidate-revision")
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()
    cells = [c for c in load_cells() if c["platform"] == args.platform]
    if args.cases:
        names = args.cases.split(",")
        if len(set(names)) != len(names) or set(names) - {c["design"] for c in cells}:
            parser.error("--cases must name unique designs in this platform suite")
        cells = [c for c in cells if c["design"] in names]
    if args.dry_run:
        print(json.dumps(cells, indent=2))
        return 0
    for name in (
        "original",
        "candidate",
        "git_repo",
        "dependency_dir",
        "output",
        "candidate_revision",
    ):
        if getattr(args, name) is None:
            parser.error(
                "--" + name.replace("_", "-") + " is required unless --dry-run"
            )
    root = args.output.expanduser().resolve() / args.platform
    if root.exists() and any(root.iterdir()):
        raise RuntimeError("Use a new empty output directory for this platform")
    root.mkdir(parents=True, exist_ok=True)
    (root / "plan.json").write_text(
        json.dumps(
            {
                "cells": cells,
                "suite_sha256": hashlib.sha256(
                    (HERE / "suite.csv").read_bytes()
                ).hexdigest(),
            },
            indent=2,
        )
        + "\n"
    )
    env = dict(os.environ)
    env.update(
        {
            key: "1"
            for key in (
                "OMP_NUM_THREADS",
                "OPENBLAS_NUM_THREADS",
                "MKL_NUM_THREADS",
                "NUMEXPR_NUM_THREADS",
            )
        }
    )
    parents = {}
    results = [dict(cell, status="pending") for cell in cells]
    save_report(root, results)
    for index, cell in enumerate(cells):
        directory = cell_directory(root, cell)
        directory.mkdir(parents=True)
        # Put the launcher log outside the cell: run.py requires an empty output.
        with (directory.parent / (directory.name + ".launcher.log")).open("w") as log:
            parent = parents.get((cell["design"], cell["controls"]))
            reuse_error = None
            try:
                reuse = reuse_option(parent)
            except (OSError, ValueError, AttributeError) as error:
                reuse, reuse_error = [], str(error)
            cmd = command(args, cell, directory, reuse)
            print(
                f"{cell['design']} K{cell['controls']} P{cell['processors']}",
                flush=True,
            )
            process = subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT, env=env)
        record = dict(
            cell,
            returncode=process.returncode,
            command=cmd,
            path=str(directory),
            reuse=reuse[0] if reuse else None,
            reuse_error=reuse_error,
            status="failed",
        )
        summary = directory / "summary.json"
        if summary.exists():
            try:
                data = json.loads(summary.read_text())
                m = data.get("cases", {}).get(cell["design"], {})
                if (
                    process.returncode == 0
                    and data.get("complete")
                    and m.get("validated")
                    and engine.number(m.get("speedup")) is not None
                    and m["speedup"] > 0
                ):
                    for mode in ("original", "candidate"):
                        if (
                            engine.number(m[mode]["median"]) is None
                            or m[mode]["median"] <= 0
                        ):
                            raise ValueError("Missing or invalid median")
                    record.update(status="passed", measurements=m)
            except (OSError, ValueError, TypeError, AttributeError, KeyError) as error:
                record["summary_error"] = str(error)
        results[index] = record
        parents[(cell["design"], cell["controls"])] = directory
        save_report(root, results)
    return 0 if save_report(root, results) else 1


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (OSError, ValueError, RuntimeError, subprocess.SubprocessError) as error:
        print("ERROR: " + str(error), file=sys.stderr)
        raise SystemExit(2)
