#!/usr/bin/env python3
"""Bounded independent-source ppmlhdfe benchmark (Python standard library)."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import platform
import shutil
import statistics
import subprocess
import sys
import time

CASES = ("chain_2", "communities_3", "imbalanced_4", "balanced_4",
         "akm_low", "akm_medium", "akm_high")
AKM_CASES = {"akm_low", "akm_medium", "akm_high"}
DEFAULT_CASES = {"imbalanced_4", "akm_low"}
STRICT_OPTIONS = "tol(1e-12) use_exact_partial(1) start_inner_tol(1e-12) itol(1e-12)"
HERE = Path(__file__).resolve().parent
# Fixtures store N/200 firm levels as Stata int (largest nonmissing value 32740).
# Round the storage bound down to the required multiple of 5000 observations.
MAX_N = (32740 * 200 // 5000) * 5000
COMPARISON_FIELDS = {
    "stripe_equal", "omitted_equal", "counts_equal", "sample_equal", "predictions_valid",
    "b_scaled", "b_se", "V_scaled", "V_cor_scaled", "V_diag_relative", "se_relative",
    "dev_scaled", "mu_scaled", "eta_max", "mu_weighted_relative", "exact_results", "passed",
}


def arguments():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--original", required=True, type=Path, help="Pristine upstream repository/archive root")
    parser.add_argument("--candidate", required=True, type=Path, help="Candidate repository root")
    parser.add_argument("--output", required=True, type=Path, help="New external artifact directory")
    parser.add_argument("--stata", type=Path, default=Path(os.environ.get("STATA", shutil.which("stata-mp") or "/Applications/Stata/StataMP.app/Contents/MacOS/stata-mp")))
    parser.add_argument("--n", type=int, default=100000, help=f"Observations; multiple of 5000, from 10000 through {MAX_N}")
    parser.add_argument("--repetitions", type=int, default=3)
    parser.add_argument("--processors", type=int, default=1, help="Required actual Stata processor count")
    parser.add_argument("--controls", type=int, default=4, help="Estimate exactly K controls (0..64)")
    parser.add_argument("--precision", choices=("tight",), default="tight", help="Common fully tight settings for both sources")
    parser.add_argument("--fixture-from", type=Path, help="Reuse a verified saved fixture only; compute a fresh one-core reference")
    parser.add_argument("--reference-from", type=Path, help="Completed run with verified one-core tight references and matching N/K/sources/dependencies")
    parser.add_argument("--dependency-dir", type=Path, help="Isolated pinned ado/Mata dependencies; drivers use only BASE, this directory and the selected source")
    parser.add_argument("--candidate-revision", help="Verify candidate src against this revision in --git-repo")
    parser.add_argument("--first-source", choices=("original", "candidate"), help="Override first source block; otherwise alternate across cases")
    parser.add_argument("--candidate-label", default="candidate", help="Recorded identity; keep committed and optimized variants separate")
    parser.add_argument("--cases", default="balanced_4")
    parser.add_argument("--fit-timeout", type=float, default=120)
    parser.add_argument("--budget", type=float, default=600, help="Elapsed Stata-stage budget in seconds, including validation; setup/finalization add overhead")
    parser.add_argument("--upstream-revision", default="b85665d8f674e93c1e07446a9a2b29be7f797910")
    parser.add_argument("--git-repo", type=Path, help="Optional Git repository containing upstream object; verifies original src file-for-file")
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()
    if not 10000 <= args.n <= MAX_N or args.n % 5000:
        parser.error(f"--n must be a multiple of 5000 between 10000 and {MAX_N}")
    if args.repetitions < 1:
        parser.error("--repetitions must be positive")
    if any(not math.isfinite(value) or value <= 0 for value in (args.fit_timeout, args.budget)):
        parser.error("fit timeout and total budget must be finite and positive")
    if args.processors < 1 or (args.controls is not None and not 0 <= args.controls <= 64):
        parser.error("--processors must be positive; --controls must be between 0 and 64")
    if args.candidate_revision and args.git_repo is None:
        parser.error("--candidate-revision requires --git-repo")
    args.cases = args.cases.split(",")
    if len(set(args.cases)) != len(args.cases) or any(case not in CASES for case in args.cases):
        parser.error("--cases must be unique comma-separated names from " + ",".join(CASES))
    if set(args.cases) & AKM_CASES and args.n != 200000:
        parser.error("AKM cases require --n 200000 (20000 workers, 1000 firms, 10 years)")
    if args.fixture_from and args.reference_from:
        parser.error("Use either --fixture-from or --reference-from, not both")
    for name in ("original", "candidate", "output", "stata"):
        path = getattr(args, name).expanduser().resolve()
        if any(character in str(path) for character in ('"', "\n", "\r")):
            parser.error("Paths containing double quotes or newlines are unsupported")
        setattr(args, name, path)
    if args.dependency_dir is not None:
        args.dependency_dir = args.dependency_dir.expanduser().resolve()
        if not args.dependency_dir.is_dir() or any(character in str(args.dependency_dir) for character in ('"', "\n", "\r")):
            parser.error("--dependency-dir must be an existing directory without quotes/newlines")
    for name in ("reference_from", "fixture_from"):
        if getattr(args, name) is not None:
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


def cpu_model():
    cpuinfo = Path("/proc/cpuinfo")
    if cpuinfo.exists():
        for line in cpuinfo.read_text(errors="replace").splitlines():
            if line.startswith("model name"):
                return line.partition(":")[2].strip()
    return platform.processor() or platform.machine()


def hashes(root):
    return {str(path.relative_to(root)): hashlib.sha256(path.read_bytes()).hexdigest()
            for path in sorted(root.rglob("*")) if path.is_file()}


def verify_source(args, source, revision):
    if args.git_repo is None:
        return {"verified_against_git": False,
                "note": "Source hashes recorded; caller must establish pristine archive provenance."}
    files = subprocess.check_output(
        ["git", "-C", str(args.git_repo), "ls-tree", "-r", "--name-only", revision, "src"],
        text=True).splitlines()
    expected = {}
    for name in files:
        content = subprocess.check_output(["git", "-C", str(args.git_repo), "show", f"{revision}:{name}"])
        expected[name.removeprefix("src/")] = hashlib.sha256(content).hexdigest()
    actual = hashes(source / "src")
    if expected != actual:
        differences = sorted(name for name in set(expected) | set(actual) if expected.get(name) != actual.get(name))
        raise RuntimeError("Source does not match requested revision: " + ", ".join(differences))
    return {"verified_against_git": True, "revision": revision}


def verify_original(args):
    return verify_source(args, args.original, args.upstream_revision)


def read_environment(directory, dependency_dir=None):
    with (directory / "environment.tsv").open(newline="") as handle:
        metadata = {row["key"]: row["value"] for row in csv.DictReader(handle, delimiter="\t")}
    with (directory / "dependencies.tsv").open(newline="") as handle:
        dependencies = {row["name"]: {"path": row["path"],
                        "sha256": hashlib.sha256(Path(row["path"]).read_bytes()).hexdigest()}
                        for row in csv.DictReader(handle, delimiter="\t")}
    if dependency_dir is not None:
        for name, item in dependencies.items():
            try:
                Path(item["path"]).resolve().relative_to(dependency_dir)
            except ValueError:
                raise RuntimeError("Dependency resolved outside pinned directory: " + name)
    return {"runtime": metadata, "dependencies": dependencies}


def dependency_hashes(environment):
    return {name: item["sha256"] for name, item in environment["dependencies"].items()}


def validate_fixture(args, manifest, source):
    """Fixture reuse does not require a successful parent model or reference."""
    prior = json.loads((source / "manifest.json").read_text())
    for key in ("N", "control_count"):
        if prior["plan"].get(key) != manifest["plan"].get(key):
            raise RuntimeError("Fixture study differs in " + key)
    generator = manifest["harness_sha256"]["fixtures.do"]
    if prior["harness_sha256"].get("fixtures.do") != generator:
        raise RuntimeError("Fixture generator differs")
    for case in args.cases:
        if case not in prior["plan"]["cases"]:
            raise RuntimeError("Missing fixture case: " + case)
        fixture = source / "fixtures" / (case + ".dta")
        if hashlib.sha256(fixture.read_bytes()).hexdigest() != prior["fixture_sha256"].get(fixture.name):
            raise RuntimeError("Reference fixture changed: " + case)
    return prior


def validate_reference(args, manifest):
    """Reject stale or incompatible saved reference evidence before any timed fit."""
    source = args.reference_from
    prior = json.loads((source / "manifest.json").read_text())
    result = json.loads((source / "summary.json").read_text())
    if not result.get("complete") or prior["plan"].get("tight_reference_processors") != 1:
        raise RuntimeError("--reference-from must be a complete validated study with one-core tight references")
    for key in ("N", "control_count"):
        if prior["plan"].get(key) != manifest["plan"].get(key):
            raise RuntimeError("Reference study differs in " + key)
    if prior["original_src_sha256"] != manifest["original_src_sha256"]:
        raise RuntimeError("Reference pristine source differs")
    if dependency_hashes(prior["environment"]) != dependency_hashes(manifest["environment"]):
        raise RuntimeError("Reference dependencies differ")
    if prior.get("pinned_dependencies_sha256") != manifest.get("pinned_dependencies_sha256"):
        raise RuntimeError("Reference complete pinned dependencies differ")
    if not manifest.get("stata_sha256") or prior.get("stata_sha256") != manifest["stata_sha256"]:
        raise RuntimeError("Reference Stata executable differs or lacks a recorded hash")
    if prior.get("harness_sha256") != manifest.get("harness_sha256"):
        raise RuntimeError("Reference benchmark harness differs")
    if prior["environment"]["runtime"]["stata_version"] != manifest["environment"]["runtime"]["stata_version"]:
        raise RuntimeError("Reference Stata version differs")
    if prior.get("thread_environment") != manifest.get("thread_environment"):
        raise RuntimeError("Reference auxiliary thread environment differs")
    for case in args.cases:
        if case not in prior["plan"]["tight_cases"]:
            raise RuntimeError("Missing one-core tight reference for " + case)
        if prior["plan"].get("precision_options", {}).get(case) != manifest["plan"]["precision_options"][case]:
            raise RuntimeError("Reference precision settings differ: " + case)
        # The parent's timed fits may use several processors. Only its saved
        # tight reference must be an actual validated one-core fit.
        tight_result = dict(result.get("runs", {}).get(case + "_tight", {}))
        if tight_result.get("status") != "ok":
            raise RuntimeError("Reference tight fit was not validated: " + case)
        tight_rows = read_csv(source / "fits" / (case + "_tight") / "timings.csv")
        check_fit(tight_result, tight_rows, case, "tight", 1, 1,
                  manifest["plan"]["controls_by_case"][case])
        if tight_result["status"] != "ok":
            raise RuntimeError("Reference tight fit has invalid one-core timing evidence: " + case)
        fixture = source / "fixtures" / (case + ".dta")
        if hashlib.sha256(fixture.read_bytes()).hexdigest() != prior["fixture_sha256"][fixture.name]:
            raise RuntimeError("Reference fixture changed: " + case)
        reference = source / "fits" / (case + "_tight") / "fit_1"
        for suffix in (".ster", ".dta"):
            artifact = reference.with_suffix(suffix)
            if hashlib.sha256(artifact.read_bytes()).hexdigest() != prior["reference_sha256"][case][artifact.name]:
                raise RuntimeError("Reference artifact changed: " + case)
    return prior


def read_csv(path):
    if not path.exists():
        return []
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def number(value):
    try:
        result = float(value)
        return result if math.isfinite(result) else None
    except (TypeError, ValueError):
        return None


def timing_error(rows, case, mode, repetitions, processors, controls):
    """Require each planned fit exactly once, with the requested specification."""
    expected = ({("validation", 1)} if mode in ("default", "tight") else
                {("warmup", 0)} | {("timed", rep) for rep in range(1, repetitions + 1)})
    actual = [(row.get("phase"), number(row.get("rep"))) for row in rows]
    if len(actual) != len(expected) or set(actual) != expected:
        return "invalid_timing_rows"
    for row in rows:
        if (row.get("design") != case or row.get("mode") != mode or
                number(row.get("processors")) != processors or
                number(row.get("control_count")) != controls):
            return "invalid_timing_metadata"
        if number(row.get("rc")) != 0 or number(row.get("converged")) != 1:
            return "fit_failed"
        seconds = number(row.get("seconds"))
        if row["phase"] != "validation" and (seconds is None or seconds <= 0):
            return "invalid_timing_seconds"
        if mode in ("candidate", "default"):
            calls = [number(row.get(key)) for key in ("schur_calls", "schur_col_calls")]
            if any(value is None or value < 0 or not value.is_integer() for value in calls):
                return "invalid_schur_diagnostics"
            if mode == "candidate" and sum(calls) == 0:
                return "schur_inactive"
            if mode == "default" and sum(calls) != 0:
                return "unexpected_schur_attempt"
    return None


def check_fit(result, rows, case, mode, repetitions, processors, controls):
    error = timing_error(rows, case, mode, repetitions, processors, controls)
    if error:
        result["timing_error"] = error
        # Keep a timeout or process failure as the primary diagnosis.
        if result["status"] == "ok":
            result["status"] = error


def comparison_error(rows, case, plan):
    """A passing subset, or duplicated checks, cannot validate a whole case."""
    expected = {(mode, rep, reference)
                for reference in (["original", "tight"] if case in plan["tight_cases"] else ["original"])
                for mode in ("original", "candidate")
                for rep in range(1, plan["repetitions"] + 1)}
    if case in plan["default_cases"]:
        expected.add(("default", 1, "original"))
    actual = [(row.get("mode"), row.get("rep"), row.get("reference")) for row in rows]
    if len(actual) != len(expected) or set(actual) != expected or any(row.get("case") != case for row in rows):
        return "incomplete_comparisons"
    if any(row.get("status") != "ok" or row.get("passed") != 1 for row in rows):
        return "failed_comparisons"
    return None


def run_stata(args, script, script_args, directory, deadline, marker):
    directory.mkdir(parents=True, exist_ok=True)
    # A no-argument driver avoids platform-specific batch-log names derived
    # from trailing command arguments and opens a license-free explicit log.
    log = directory / "run.log"
    driver = directory / "driver.do"
    do_command = "do " + " ".join('"' + str(value) + '"' for value in [HERE / script, *script_args])
    isolated_ado = ""
    if getattr(args, "dependency_dir", None) is not None:
        isolated_ado = 'global S_ADO "BASE"\nadopath ++ "' + str(args.dependency_dir) + '"\n'
    driver.write_text('capture log close _all\nlog using "' + str(log) + '", text replace\n'
                      'display "REDBENCH_START"\n' + isolated_ado + do_command + '\nlog close _all\n')
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
    status = timeout or ("ok" if process.returncode == 0 and marker in log_text else "script_failed")
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
                if len(rows) != 1 or set(rows[0]) != COMPARISON_FIELDS:
                    if comparison["status"] == "ok":
                        comparison["status"] = "invalid_comparison_rows"
                else:
                    metrics = {key: number(value) for key, value in rows[0].items()}
                    comparison.update(metrics)
                    if comparison["status"] == "ok" and (any(value is None for value in metrics.values()) or metrics["passed"] != 1):
                        comparison["status"] = "accuracy_failed"
            comparisons.append(comparison)
    return comparisons


def main():
    args = arguments()
    plan = {"cases": args.cases, "N": args.n, "repetitions": args.repetitions,
            "processors": args.processors, "control_count": args.controls,
            "controls_by_case": {case: args.controls for case in args.cases},
            "tight_reference_processors": 1, "candidate_label": args.candidate_label,
            "source_order": {case: (["original", "candidate"] if (args.first_source == "original" or (args.first_source is None and index % 2 == 0)) else ["candidate", "original"])
                             for index, case in enumerate(args.cases)},
            "tight_cases": sorted(args.cases),
            "default_cases": sorted(set(args.cases) & DEFAULT_CASES),
            "precision_options": {case: STRICT_OPTIONS
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
                "hostname": platform.node(), "cpu_model": cpu_model(),
                "thread_environment": {name: os.environ.get(name) for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS")},
                "original": str(args.original), "candidate": str(args.candidate),
                "stata": str(args.stata), "stata_sha256": hashlib.sha256(args.stata.read_bytes()).hexdigest(),
                "upstream": provenance, "plan": plan,
                "dependency_dir": str(args.dependency_dir) if args.dependency_dir else None,
                "pinned_dependencies_sha256": hashes(args.dependency_dir) if args.dependency_dir else None,
                "candidate_provenance": verify_source(args, args.candidate, args.candidate_revision) if args.candidate_revision else {"verified_against_git": False},
                "original_src_sha256": hashes(args.original / "src"),
                "candidate_src_sha256": hashes(args.candidate / "src"),
                "harness_sha256": {name: hashlib.sha256((HERE / name).read_bytes()).hexdigest()
                                   for name in ("run.py", "fixtures.do", "fit.do", "compare.do", "environment.do")}}
    save_json(args.output / "manifest.json", manifest)
    deadline = time.monotonic() + args.budget
    environment_dir = args.output / "environment"
    environment_run = run_stata(args, "environment.do", [environment_dir, args.processors, args.original], environment_dir, deadline, "REDBENCH_ENVIRONMENT_OK")
    if environment_run["status"] != "ok":
        save_json(args.output / "summary.json", {"status": "environment_failure", "process": environment_run})
        return 1
    manifest["environment"] = read_environment(environment_dir, args.dependency_dir)
    if int(float(manifest["environment"]["runtime"]["processors"])) != args.processors:
        raise RuntimeError("Actual processor count differs from requested count")
    fixtures = args.output / "fixtures"
    fixtures.mkdir()
    fixture_source = args.reference_from or getattr(args, "fixture_from", None)
    if fixture_source:
        if args.reference_from:
            validate_reference(args, manifest)
        else:
            validate_fixture(args, manifest, fixture_source)
        for case in args.cases:
            shutil.copy2(fixture_source / "fixtures" / (case + ".dta"), fixtures)
        key = "reference_from" if args.reference_from else "fixture_from"
        manifest[key] = {"path": str(fixture_source), "manifest_sha256": hashlib.sha256((fixture_source / "manifest.json").read_bytes()).hexdigest()}
    else:
        generated = run_stata(args, "fixtures.do", [fixtures, args.n, max(4, args.controls or 4), " ".join(args.cases)], fixtures, deadline, "REDBENCH_FIXTURES_OK")
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
            result = run_stata(args, "fit.do", [source, fixtures / f"{case}.dta", directory, mode, args.repetitions, plan["precision_options"][case], args.processors, args.controls if args.controls is not None else ""], directory, deadline, "REDBENCH_FITS_OK")
            runs[f"{case}_{mode}"] = result
            rows = read_csv(directory / "timings.csv")
            check_fit(result, rows, case, mode, args.repetitions, args.processors, plan["controls_by_case"][case])
            timings.extend(rows)
    # Validate ordinary fits before expensive tight references, preserving
    # useful diagnostics even when a later reference consumes the time budget.
    comparisons = []
    for case in args.cases:
        comparisons.extend(compare_fits(args, case, "original", ("original", "candidate"), deadline))
    # Untimed candidate controls require exact equality to pristine output.
    for case in plan["default_cases"]:
        directory = args.output / "fits" / f"{case}_default"
        runs[f"{case}_default"] = run_stata(args, "fit.do", [args.candidate, fixtures / f"{case}.dta", directory, "default", 1, plan["precision_options"][case], args.processors, args.controls if args.controls is not None else ""], directory, deadline, "REDBENCH_FITS_OK")
        rows = read_csv(directory / "timings.csv")
        check_fit(runs[f"{case}_default"], rows, case, "default", 1, args.processors, plan["controls_by_case"][case])
        timings.extend(rows)
        comparisons.extend(compare_fits(args, case, "original", ("default",), deadline))
    # Final speedup gates still require every planned tight-reference check.
    for case in plan["tight_cases"]:
        directory = args.output / "fits" / f"{case}_tight"
        if args.reference_from:
            previous = args.reference_from / "fits" / f"{case}_tight"
            shutil.copytree(previous, directory)
            runs[f"{case}_tight"] = {"status": "ok", "reused_from": str(previous), "processors": 1}
        else:
            runs[f"{case}_tight"] = run_stata(args, "fit.do", [args.original, fixtures / f"{case}.dta", directory, "tight", 1, STRICT_OPTIONS, 1, args.controls if args.controls is not None else ""], directory, deadline, "REDBENCH_FITS_OK")
        rows = read_csv(directory / "timings.csv")
        check_fit(runs[f"{case}_tight"], rows, case, "tight", 1, 1, plan["controls_by_case"][case])
        timings.extend(rows)
        comparisons.extend(compare_fits(args, case, "tight", ("original", "candidate"), deadline))
    final_environment_dir = args.output / "environment_after"
    environment_after = run_stata(args, "environment.do", [final_environment_dir, args.processors, args.original], final_environment_dir, deadline, "REDBENCH_ENVIRONMENT_OK")
    unchanged_dependencies = environment_after["status"] == "ok" and read_environment(final_environment_dir, args.dependency_dir) == manifest["environment"]
    if args.dependency_dir:
        unchanged_dependencies = unchanged_dependencies and manifest["pinned_dependencies_sha256"] == hashes(args.dependency_dir)
    manifest["reference_sha256"] = {}
    for case in plan["tight_cases"]:
        reference = args.output / "fits" / f"{case}_tight" / "fit_1"
        manifest["reference_sha256"][case] = {reference.with_suffix(suffix).name: hashlib.sha256(reference.with_suffix(suffix).read_bytes()).hexdigest() for suffix in (".ster", ".dta") if reference.with_suffix(suffix).exists()}
    save_json(args.output / "manifest.json", manifest)
    summary = {"runs": runs, "comparisons": comparisons, "cases": {}, "dependencies_unchanged": unchanged_dependencies}
    summary["fixtures_unchanged"] = manifest["fixture_sha256"] == {
        path.name: hashlib.sha256(path.read_bytes()).hexdigest() for path in fixtures.glob("*.dta")}
    summary["stata_unchanged"] = manifest["stata_sha256"] == hashlib.sha256(args.stata.read_bytes()).hexdigest()
    summary["harness_unchanged"] = all(digest == hashlib.sha256((HERE / name).read_bytes()).hexdigest()
                                       for name, digest in manifest["harness_sha256"].items())
    summary["sources_unchanged"] = (manifest["original_src_sha256"] == hashes(args.original / "src") and manifest["candidate_src_sha256"] == hashes(args.candidate / "src"))
    for case in args.cases:
        case_summary = {}
        for mode in ("original", "candidate"):
            rows = [row for row in timings if row.get("design") == case and row.get("mode") == mode and row.get("phase") == "timed"]
            valid = [number(row.get("seconds")) for row in rows
                     if number(row.get("rc")) == 0 and number(row.get("converged")) == 1
                     and number(row.get("seconds")) is not None and number(row.get("seconds")) > 0]
            case_summary[mode] = {"timed_fits": len(rows), "successful_fits": len(valid)}
            if len(valid) == args.repetitions:
                case_summary[mode].update({"median": statistics.median(valid), "minimum": min(valid), "maximum": max(valid)})
        both_complete = all("median" in case_summary[mode] for mode in ("original", "candidate"))
        checked = [row for row in comparisons if row["case"] == case]
        relevant_runs = [run for name, run in runs.items() if name.startswith(case + "_")]
        case_summary["comparison_error"] = comparison_error(checked, case, plan)
        accurate = case_summary["comparison_error"] is None
        accurate = (accurate and all(run["status"] == "ok" for run in relevant_runs)
                    and unchanged_dependencies and summary["sources_unchanged"]
                    and summary["fixtures_unchanged"] and summary["stata_unchanged"] and summary["harness_unchanged"])
        case_summary["validated"] = bool(accurate)
        if both_complete and accurate:
            case_summary["speedup"] = case_summary["original"]["median"] / case_summary["candidate"]["median"]
        summary["cases"][case] = case_summary
    summary["complete"] = summary["sources_unchanged"] and all("speedup" in case for case in summary["cases"].values())
    save_json(args.output / "summary.json", summary)
    for row in timings:
        row.update(hostname=manifest["hostname"], cpu_model=manifest["cpu_model"],
                   requested_processors=1 if row.get("mode") == "tight" else args.processors,
                   precision_options=STRICT_OPTIONS if row.get("mode") == "tight" else plan["precision_options"].get(row.get("design")),
                   candidate_label=args.candidate_label)
    for row in comparisons:
        row.update(processors=args.processors, control_count=plan["controls_by_case"][row["case"]],
                   reference_processors=1 if row["reference"] == "tight" else args.processors)
    for name, rows in (("timings.csv", timings), ("comparisons.csv", comparisons)):
        fields = list(dict.fromkeys(key for row in rows for key in row))
        with (args.output / name).open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields)
            writer.writeheader()
            writer.writerows(rows)
    report = ["# Independent-source ppmlhdfe benchmark", "", f"N={args.n:,}; processors={args.processors}; controls={plan['controls_by_case']}; one warmup and {args.repetitions} timed fits per source/case.",
              "", "| Case | Original median [range], s | Candidate median [range], s | Validated speedup |", "|---|---:|---:|---:|"]
    for case, values in summary["cases"].items():
        cells = []
        for mode in ("original", "candidate"):
            row = values[mode]
            cells.append(f"{row['median']:.3f} [{row['minimum']:.3f}, {row['maximum']:.3f}]" if "median" in row else f"Incomplete ({row['successful_fits']}/{args.repetitions})")
        speedup = f"{values['speedup']:.3f}x" if "speedup" in values else "Not validated"
        report.append(f"| {case} | {cells[0]} | {cells[1]} | {speedup} |")
    report += ["", "All timings, failures, precision checks, source/fixture hashes, and process outcomes are retained in CSV/JSON/log files.",
               "Common precision options by case: " + json.dumps(plan["precision_options"]) + ".",
               "Tight pristine references use one processor and are untimed; every primary fit in each reference case is checked against them.",
               "Untimed candidate use_schur(0) controls require exact equality to pristine results in the imbalanced and low-mobility AKM cases.",
               "This is a bounded development benchmark; it does not establish universal speedup or release qualification.",
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
