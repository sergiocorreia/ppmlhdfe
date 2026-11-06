#!/usr/bin/env python3
"""Bounded independent-source ppmlhdfe benchmark (Python standard library)."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import platform
import shutil
import statistics
import subprocess
import sys
import time

CASES = (
    "random_2", "chain_2", "disconnected_2", "random_3", "communities_3",
    "redundant_3", "imbalanced_4", "balanced_4",
)
TIGHT_CASES = {"chain_2", "random_3", "communities_3", "imbalanced_4"}
DEFAULT_CASES = {"random_2", "random_3", "imbalanced_4"}
STRICT_TIMING_CASES = {"random_3", "communities_3", "imbalanced_4"}
STRICT_OPTIONS = "tol(1e-12) use_exact_partial(1) start_inner_tol(1e-12) itol(1e-12)"
HERE = Path(__file__).resolve().parent


def arguments():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--original", required=True, type=Path, help="Pristine upstream repository/archive root")
    parser.add_argument("--candidate", required=True, type=Path, help="Candidate repository root")
    parser.add_argument("--output", required=True, type=Path, help="New external artifact directory")
    parser.add_argument("--stata", type=Path, default=Path(os.environ.get("STATA", shutil.which("stata-mp") or "/Applications/Stata/StataMP.app/Contents/MacOS/stata-mp")))
    parser.add_argument("--n", type=int, default=30000, help="Observations; multiple of 5000, at least 10000")
    parser.add_argument("--repetitions", type=int, default=3)
    parser.add_argument("--cases", default=",".join(CASES))
    parser.add_argument("--fit-timeout", type=float, default=120)
    parser.add_argument("--budget", type=float, default=600, help="Whole-run seconds, including validation")
    parser.add_argument("--upstream-revision", default="b85665d8f674e93c1e07446a9a2b29be7f797910")
    parser.add_argument("--git-repo", type=Path, help="Optional Git repository containing upstream object; verifies original src file-for-file")
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()
    if args.n < 10000 or args.n % 5000:
        parser.error("--n must be a multiple of 5000 and at least 10000")
    if args.repetitions < 1 or args.fit_timeout <= 0 or args.budget <= 0:
        parser.error("repetitions, fit timeout and total budget must be positive")
    args.cases = args.cases.split(",")
    if len(set(args.cases)) != len(args.cases) or any(case not in CASES for case in args.cases):
        parser.error("--cases must be unique comma-separated names from " + ",".join(CASES))
    for name in ("original", "candidate", "output", "stata"):
        path = getattr(args, name).expanduser().resolve()
        if any(character in str(path) for character in ('"', "\n", "\r")):
            parser.error("Paths containing double quotes or newlines are unsupported")
        setattr(args, name, path)
    for path in (args.original, args.candidate):
        if not (path / "src" / "ppmlhdfe.ado").is_file():
            parser.error(f"Missing ppmlhdfe source: {path}")
    if not args.stata.is_file():
        parser.error(f"Stata executable not found: {args.stata}")
    return args


def hashes(root):
    return {str(path.relative_to(root)): hashlib.sha256(path.read_bytes()).hexdigest()
            for path in sorted(root.rglob("*")) if path.is_file()}


def verify_original(args):
    if args.git_repo is None:
        return {"verified_against_git": False,
                "note": "Source hashes recorded; caller must establish pristine archive provenance."}
    files = subprocess.check_output(
        ["git", "-C", str(args.git_repo), "ls-tree", "-r", "--name-only", args.upstream_revision, "src"],
        text=True).splitlines()
    expected = {}
    for name in files:
        content = subprocess.check_output(["git", "-C", str(args.git_repo), "show", f"{args.upstream_revision}:{name}"])
        expected[name.removeprefix("src/")] = hashlib.sha256(content).hexdigest()
    actual = hashes(args.original / "src")
    if expected != actual:
        differences = sorted(name for name in set(expected) | set(actual) if expected.get(name) != actual.get(name))
        raise RuntimeError("Original source does not match upstream revision: " + ", ".join(differences))
    return {"verified_against_git": True, "revision": args.upstream_revision}


def read_csv(path):
    if not path.exists():
        return []
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def number(value):
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def run_stata(args, script, script_args, directory, deadline, marker):
    directory.mkdir(parents=True, exist_ok=True)
    # A no-argument driver avoids platform-specific batch-log names derived
    # from trailing command arguments and opens a license-free explicit log.
    log = directory / "run.log"
    driver = directory / "driver.do"
    do_command = "do " + " ".join('"' + str(value) + '"' for value in [HERE / script, *script_args])
    driver.write_text('capture log close _all\nlog using "' + str(log) + '", text replace\n'
                      'display "REDBENCH_START"\n' + do_command + '\nlog close _all\n')
    command = [str(args.stata), "-b", "do", str(driver)]
    if time.monotonic() >= deadline:
        return {"status": "budget_exhausted", "command": command}
    started = time.monotonic()
    current_stage = "startup"
    stage_started = started
    stage_path = directory / "stage.txt"
    process = subprocess.Popen(command, cwd=directory, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    timeout = None
    while process.poll() is None:
        now = time.monotonic()
        if stage_path.exists():
            stage = stage_path.read_text().strip()
            if stage and stage != current_stage:
                current_stage, stage_started = stage, now
        if now >= deadline:
            timeout = "budget_exhausted"
        elif now - stage_started > args.fit_timeout:
            timeout = "fit_timeout"
        if timeout:
            process.terminate()
            try:
                process.wait(timeout=5)
            except subprocess.TimeoutExpired:
                process.kill()
                process.wait()
            break
        time.sleep(0.2)
    log_text = log.read_text(errors="replace") if log.exists() else ""
    # Never retain Stata's opening license banner in shareable benchmark evidence.
    start = log_text.find("REDBENCH_START")
    if start >= 0:
        log_text = log_text[start:]
    else:
        log_text = "[No benchmark start marker; pre-script output omitted.]\n"
    log.write_text(log_text)
    for batch_log in directory.glob("*.log"):
        if batch_log != log:
            batch_log.unlink()
    status = timeout or ("ok" if marker in log_text else "script_failed")
    result = {"status": status, "returncode": process.returncode,
              "wall_seconds": time.monotonic() - started, "last_stage": current_stage,
              "command": command, "log": str(log)}
    (directory / "process.json").write_text(json.dumps(result, indent=2) + "\n")
    print(f"{directory.name}: {status} ({result['wall_seconds']:.1f}s)", flush=True)
    return result


def save_json(path, value):
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")


def compare_fits(args, case, reference_mode, modes, deadline):
    reference = args.output / "fits" / f"{case}_{reference_mode}" / "fit_1"
    comparisons = []
    for mode in modes:
        repetitions = 1 if mode == "default" else args.repetitions
        for repetition in range(1, repetitions + 1):
            target = args.output / "fits" / f"{case}_{mode}" / f"fit_{repetition}"
            label = f"{case}_{mode}_{repetition}_vs_{reference_mode}"
            directory = args.output / "comparisons" / label
            comparison = {"case": case, "mode": mode, "rep": repetition, "reference": reference_mode}
            if not all(path.with_suffix(suffix).exists() for path in (reference, target) for suffix in (".ster", ".dta")):
                comparison["status"] = "missing_fit"
            else:
                result = run_stata(args, "compare.do", [reference, target, directory / "comparison.csv", int(mode == "default")], directory, deadline, "REDBENCH_COMPARE_OK")
                comparison["status"] = result["status"]
                rows = read_csv(directory / "comparison.csv")
                if rows:
                    comparison.update({key: number(value) for key, value in rows[0].items()})
                    if comparison["passed"] != 1:
                        comparison["status"] = "accuracy_failed"
            comparisons.append(comparison)
    return comparisons


def main():
    args = arguments()
    plan = {"cases": args.cases, "N": args.n, "repetitions": args.repetitions,
            "source_order": {case: (["original", "candidate"] if index % 2 == 0 else ["candidate", "original"])
                             for index, case in enumerate(args.cases)},
            "tight_cases": sorted(set(args.cases) & TIGHT_CASES),
            "default_cases": sorted(set(args.cases) & DEFAULT_CASES),
            "precision_options": {case: STRICT_OPTIONS if case in STRICT_TIMING_CASES else "tol(1e-10)"
                                  for case in args.cases},
            "budget_seconds": args.budget, "fit_timeout_seconds": args.fit_timeout}
    if args.dry_run:
        print(json.dumps(plan, indent=2))
        return 0
    if args.output.exists() and any(args.output.iterdir()):
        raise RuntimeError("Use a new empty output directory; existing evidence is never overwritten.")
    args.output.mkdir(parents=True, exist_ok=True)
    provenance = verify_original(args)
    manifest = {"created_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
                "platform": platform.platform(), "python": sys.version,
                "original": str(args.original), "candidate": str(args.candidate),
                "stata": str(args.stata), "upstream": provenance, "plan": plan,
                "original_src_sha256": hashes(args.original / "src"),
                "candidate_src_sha256": hashes(args.candidate / "src"),
                "harness_sha256": {name: hashlib.sha256((HERE / name).read_bytes()).hexdigest()
                                   for name in ("run.py", "fixtures.do", "fit.do", "compare.do")}}
    save_json(args.output / "manifest.json", manifest)
    deadline = time.monotonic() + args.budget
    fixtures = args.output / "fixtures"
    fixtures.mkdir()
    generated = run_stata(args, "fixtures.do", [fixtures, args.n], fixtures, deadline, "REDBENCH_FIXTURES_OK")
    if generated["status"] != "ok":
        save_json(args.output / "summary.json", {"status": "fixture_failure", "process": generated})
        return 1
    manifest["fixture_sha256"] = {path.name: hashlib.sha256(path.read_bytes()).hexdigest() for path in fixtures.glob("*.dta")}
    save_json(args.output / "manifest.json", manifest)
    runs = {}
    timings = []
    for case in args.cases:
        for mode in plan["source_order"][case]:
            directory = args.output / "fits" / f"{case}_{mode}"
            source = args.original if mode == "original" else args.candidate
            result = run_stata(args, "fit.do", [source, fixtures / f"{case}.dta", directory, mode, args.repetitions, plan["precision_options"][case]], directory, deadline, "REDBENCH_FITS_OK")
            runs[f"{case}_{mode}"] = result
            rows = read_csv(directory / "timings.csv")
            if any(number(row["rc"]) != 0 or number(row["converged"]) != 1 for row in rows):
                result["status"] = "fit_failed"
            if mode == "candidate" and any(
                (number(row.get("schur_calls")) or 0) + (number(row.get("schur_col_calls")) or 0) == 0
                for row in rows if number(row["rc"]) == 0
            ):
                result["status"] = "schur_inactive"
            timings.extend(rows)
    # Validate ordinary fits before expensive tight references, preserving
    # useful diagnostics even when a later reference consumes the time budget.
    comparisons = []
    for case in args.cases:
        comparisons.extend(compare_fits(args, case, "original", ("original", "candidate"), deadline))
    # Untimed candidate controls require exact equality to pristine output.
    for case in plan["default_cases"]:
        directory = args.output / "fits" / f"{case}_default"
        runs[f"{case}_default"] = run_stata(args, "fit.do", [args.candidate, fixtures / f"{case}.dta", directory, "default", 1, plan["precision_options"][case]], directory, deadline, "REDBENCH_FITS_OK")
        rows = read_csv(directory / "timings.csv")
        if any(number(row["rc"]) != 0 or number(row["converged"]) != 1 for row in rows):
            runs[f"{case}_default"]["status"] = "fit_failed"
        if any((number(row.get("schur_calls")) or 0) != 0 or (number(row.get("schur_col_calls")) or 0) != 0 for row in rows):
            runs[f"{case}_default"]["status"] = "unexpected_schur_attempt"
        timings.extend(rows)
        comparisons.extend(compare_fits(args, case, "original", ("default",), deadline))
    # Final speedup gates still require every planned tight-reference check.
    for case in plan["tight_cases"]:
        directory = args.output / "fits" / f"{case}_tight"
        runs[f"{case}_tight"] = run_stata(args, "fit.do", [args.original, fixtures / f"{case}.dta", directory, "tight", 1, STRICT_OPTIONS], directory, deadline, "REDBENCH_FITS_OK")
        timings.extend(read_csv(directory / "timings.csv"))
        comparisons.extend(compare_fits(args, case, "tight", ("original", "candidate"), deadline))
    summary = {"runs": runs, "comparisons": comparisons, "cases": {}}
    for case in args.cases:
        case_summary = {}
        for mode in ("original", "candidate"):
            rows = [row for row in timings if row["design"] == case and row["mode"] == mode and row["phase"] == "timed"]
            valid = [number(row["seconds"]) for row in rows if number(row["rc"]) == 0 and number(row["converged"]) == 1]
            case_summary[mode] = {"timed_fits": len(rows), "successful_fits": len(valid)}
            if len(valid) == args.repetitions:
                case_summary[mode].update({"median": statistics.median(valid), "minimum": min(valid), "maximum": max(valid)})
        both_complete = all("median" in case_summary[mode] for mode in ("original", "candidate"))
        checked = [row for row in comparisons if row["case"] == case]
        relevant_runs = [run for name, run in runs.items() if name.startswith(case + "_")]
        accurate = checked and all(row.get("status") == "ok" and row.get("passed") == 1 for row in checked)
        accurate = accurate and all(run["status"] == "ok" for run in relevant_runs)
        case_summary["validated"] = bool(accurate)
        if both_complete and accurate:
            case_summary["speedup"] = case_summary["original"]["median"] / case_summary["candidate"]["median"]
        summary["cases"][case] = case_summary
    summary["sources_unchanged"] = (manifest["original_src_sha256"] == hashes(args.original / "src") and manifest["candidate_src_sha256"] == hashes(args.candidate / "src"))
    summary["complete"] = summary["sources_unchanged"] and all("speedup" in case for case in summary["cases"].values())
    save_json(args.output / "summary.json", summary)
    for name, rows in (("timings.csv", timings), ("comparisons.csv", comparisons)):
        fields = list(dict.fromkeys(key for row in rows for key in row))
        with (args.output / name).open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields)
            writer.writeheader()
            writer.writerows(rows)
    report = ["# Independent-source ppmlhdfe benchmark", "", f"N={args.n:,}; one processor; one warmup and {args.repetitions} timed fits per source/case.",
              "", "| Case | Original median [range], s | Candidate median [range], s | Validated speedup |", "|---|---:|---:|---:|"]
    for case, values in summary["cases"].items():
        cells = []
        for mode in ("original", "candidate"):
            row = values[mode]
            cells.append(f"{row['median']:.3f} [{row['minimum']:.3f}, {row['maximum']:.3f}]" if "median" in row else f"Incomplete ({row['successful_fits']}/{args.repetitions})")
        speedup = f"{values['speedup']:.3f}x" if "speedup" in values else "Not validated"
        report.append(f"| {case} | {cells[0]} | {cells[1]} | {speedup} |")
    report += ["", "All timings, failures, precision checks, source/fixture hashes, and process outcomes are retained in CSV/JSON/log files.",
               "Both sources use tol(1e-10), except random_3, communities_3 and imbalanced_4, which use identical fully tight settings: " + STRICT_OPTIONS + ".",
               "Tight pristine references are untimed; every primary fit in each reference case is checked against them.",
               "Untimed candidate use_schur(0) controls require exact equality to pristine results in selected 2FE, 3FE and 4FE cases.",
               "This is a bounded local development benchmark; it does not establish universal speedup or release qualification.",
               "", f"Complete and validated: {summary['complete']}. Sources unchanged during run: {summary['sources_unchanged']}."]
    (args.output / "REPORT.md").write_text("\n".join(report) + "\n")
    print(f"REPORT: {args.output / 'REPORT.md'}", flush=True)
    return 0 if summary["complete"] else 1


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (OSError, RuntimeError, subprocess.SubprocessError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        raise SystemExit(2)
