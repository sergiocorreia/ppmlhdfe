"""Standard-library adversarial checks for benchmark evidence acceptance."""
import copy
import csv
import hashlib
import importlib.util
import io
import json
from pathlib import Path
import tempfile
import time
from types import SimpleNamespace
import unittest
from unittest.mock import patch

spec = importlib.util.spec_from_file_location("schur_runner", Path(__file__).with_name("run.py"))
runner = importlib.util.module_from_spec(spec)
spec.loader.exec_module(runner)


def timing_rows(case="balanced_4", mode="candidate", repetitions=3, processors=4, controls=16):
    identities = ([("validation", 1)] if mode in ("tight", "default") else
                  [("warmup", 0)] + [("timed", rep) for rep in range(1, repetitions + 1)])
    return [{"design": case, "mode": mode, "phase": phase, "rep": str(rep),
             "rc": "0", "converged": "1", "seconds": "." if phase == "validation" else "1.5",
             "processors": str(processors), "control_count": str(controls),
             "schur_calls": "5" if mode == "candidate" else "0", "schur_col_calls": "0"}
            for phase, rep in identities]


def comparison_metrics():
    return {field: "1" if field in {"stripe_equal", "omitted_equal", "counts_equal",
            "sample_equal", "predictions_valid", "exact_results", "passed"} else "0"
            for field in runner.COMPARISON_FIELDS}


def write_csv(path, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]) if rows else ["passed"])
        writer.writeheader()
        writer.writerows(rows)


class ArgumentValidationTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        root = Path(self.temporary.name)
        source = root / "source"
        (source / "src").mkdir(parents=True)
        (source / "src/ppmlhdfe.ado").write_text("source")
        stata = root / "stata"
        stata.write_bytes(b"not executed")
        self.argv = ["run.py", "--original", str(source), "--candidate", str(source),
                     "--stata", str(stata), "--output", str(root / "results")]

    def parse(self, *extra):
        with patch("sys.argv", self.argv + list(extra)), patch("sys.stderr", new_callable=io.StringIO):
            return runner.arguments()

    def test_storage_bound_accepts_supported_endpoints(self):
        self.assertEqual(runner.MAX_N, 6545000)
        self.assertLessEqual(runner.MAX_N // 200, 32740)
        self.assertGreater((runner.MAX_N + 5000) // 200, 32740)
        for observations in (10000, 30000, 100000, runner.MAX_N):
            with self.subTest(n=observations):
                self.assertEqual(self.parse("--n", str(observations)).n, observations)

    def test_invalid_observation_bounds_rejected(self):
        for observations in (-10000, 0, 5000, 10001, 6548000, 6550000):
            with self.subTest(n=observations), self.assertRaises(SystemExit) as error:
                self.parse("--n", str(observations))
            self.assertEqual(error.exception.code, 2)

    def test_positive_finite_time_limits_accepted(self):
        for option, field in (("--budget", "budget"), ("--fit-timeout", "fit_timeout")):
            for value in ("0.01", "600", "1e4"):
                with self.subTest(option=option, value=value):
                    self.assertEqual(getattr(self.parse(option + "=" + value), field), float(value))

    def test_nonfinite_or_nonpositive_time_limits_rejected(self):
        for option in ("--budget", "--fit-timeout"):
            for value in ("0", "-1", "nan", "inf", "-inf"):
                with self.subTest(option=option, value=value), self.assertRaises(SystemExit) as error:
                    self.parse(option + "=" + value)
                self.assertEqual(error.exception.code, 2)


class ReferenceIntegrityTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        self.args = SimpleNamespace(reference_from=self.root, cases=["balanced_4"])
        self.manifest = {
            "plan": {"N": 10000, "control_count": 16, "processors": 1, "tight_reference_processors": 1,
                     "controls_by_case": {"balanced_4": 16}, "precision_options": {"balanced_4": runner.STRICT_OPTIONS},
                     "tight_cases": ["balanced_4"]},
            "original_src_sha256": {"ppmlhdfe.mata": "source-hash"},
            "stata_sha256": "executable-hash", "harness_sha256": {"fit.do": "fit-hash"},
            "pinned_dependencies_sha256": {"ftools.mata": "dependency-hash", "other.ado": "other-hash"},
            "thread_environment": {"OMP_NUM_THREADS": "1"},
            "environment": {"runtime": {"stata_version": "19"},
                            "dependencies": {"ftools.mata": {"sha256": "dependency-hash"}}},
            "fixture_sha256": {}, "reference_sha256": {"balanced_4": {}}}
        fixture = self.root / "fixtures/balanced_4.dta"
        fixture.parent.mkdir()
        fixture.write_bytes(b"identical fixture bytes")
        self.manifest["fixture_sha256"][fixture.name] = hashlib.sha256(fixture.read_bytes()).hexdigest()
        for suffix in (".ster", ".dta"):
            target = self.root / "fits/balanced_4_tight" / ("fit_1" + suffix)
            target.parent.mkdir(parents=True, exist_ok=True)
            target.write_bytes(("saved reference" + suffix).encode())
            self.manifest["reference_sha256"]["balanced_4"][target.name] = hashlib.sha256(target.read_bytes()).hexdigest()
        self.save(self.manifest)
        self.summary = {"complete": True, "runs": {"balanced_4_tight": {"status": "ok"}}}
        self.save_summary(self.summary)
        write_csv(self.root / "fits/balanced_4_tight/timings.csv", timing_rows(mode="tight", processors=1))

    def save(self, manifest):
        (self.root / "manifest.json").write_text(json.dumps(manifest))

    def save_summary(self, summary):
        (self.root / "summary.json").write_text(json.dumps(summary))

    def test_unchanged_reference_accepted(self):
        self.assertEqual(runner.validate_reference(self.args, self.manifest), self.manifest)

    def test_four_core_parent_with_actual_one_core_reference_accepted(self):
        prior = copy.deepcopy(self.manifest)
        prior["plan"]["processors"] = 4
        self.save(prior)
        self.assertEqual(runner.validate_reference(self.args, self.manifest), prior)

    def test_actual_multicore_reference_rejected_even_with_one_core_plan(self):
        prior = copy.deepcopy(self.manifest)
        prior["plan"]["processors"] = 4
        self.save(prior)
        write_csv(self.root / "fits/balanced_4_tight/timings.csv", timing_rows(mode="tight", processors=4))
        with self.assertRaisesRegex(RuntimeError, "invalid one-core timing evidence"):
            runner.validate_reference(self.args, self.manifest)

    def test_missing_duplicate_or_failed_tight_rows_rejected(self):
        rows = timing_rows(mode="tight", processors=1)
        for altered in ([], rows + rows, [dict(rows[0], rc="430")], [dict(rows[0], converged="0")],
                        [dict(rows[0], mode="original")], [dict(rows[0], control_count="64")]):
            with self.subTest(rows=altered):
                write_csv(self.root / "fits/balanced_4_tight/timings.csv", altered)
                with self.assertRaisesRegex(RuntimeError, "invalid one-core timing evidence"):
                    runner.validate_reference(self.args, self.manifest)
        (self.root / "fits/balanced_4_tight/timings.csv").unlink()
        with self.assertRaisesRegex(RuntimeError, "invalid one-core timing evidence"):
            runner.validate_reference(self.args, self.manifest)

    def test_missing_or_failed_tight_summary_rejected(self):
        for status in (None, "fit_failed", "script_failed", "fit_timeout"):
            with self.subTest(status=status):
                summary = copy.deepcopy(self.summary)
                if status is None:
                    summary["runs"].clear()
                else:
                    summary["runs"]["balanced_4_tight"]["status"] = status
                self.save_summary(summary)
                with self.assertRaisesRegex(RuntimeError, "tight fit was not validated"):
                    runner.validate_reference(self.args, self.manifest)

    def test_changed_artifact_rejected(self):
        (self.root / "fits/balanced_4_tight/fit_1.dta").write_bytes(b"changed")
        with self.assertRaisesRegex(RuntimeError, "artifact changed"):
            runner.validate_reference(self.args, self.manifest)

    def test_changed_fixture_rejected(self):
        (self.root / "fixtures/balanced_4.dta").write_bytes(b"changed")
        with self.assertRaisesRegex(RuntimeError, "fixture changed"):
            runner.validate_reference(self.args, self.manifest)

    def test_incompatible_reference_rejected(self):
        variants = [
            ("source", lambda m: m["original_src_sha256"].update(new="source")),
            ("dependencies", lambda m: m["environment"]["dependencies"]["ftools.mata"].update(sha256="changed")),
            ("complete pinned dependencies", lambda m: m["pinned_dependencies_sha256"].update({"other.ado": "changed"})),
            ("complete pinned dependencies", lambda m: m.pop("pinned_dependencies_sha256")),
            ("executable", lambda m: m.update(stata_sha256="changed")),
            ("executable", lambda m: m.pop("stata_sha256")),
            ("harness", lambda m: m["harness_sha256"].update({"fit.do": "changed"})),
            ("Stata version", lambda m: m["environment"]["runtime"].update(stata_version="18")),
            ("control_count", lambda m: m["plan"].update(control_count=64)),
            ("precision settings", lambda m: m["plan"]["precision_options"].update(balanced_4="tol(1e-8)")),
            ("N", lambda m: m["plan"].update(N=100000)),
            ("one-core", lambda m: m["plan"].update(tight_reference_processors=4)),
            ("one-core", lambda m: m["plan"].pop("tight_reference_processors")),
            ("thread environment", lambda m: m["thread_environment"].update(OMP_NUM_THREADS="4")),
            ("Missing one-core", lambda m: m["plan"].update(tight_cases=[])),
        ]
        for expected, mutate in variants:
            with self.subTest(expected=expected, mutate=mutate):
                prior = copy.deepcopy(self.manifest)
                mutate(prior)
                self.save(prior)
                with self.assertRaisesRegex(RuntimeError, expected):
                    runner.validate_reference(self.args, self.manifest)

    def test_incomplete_reference_rejected(self):
        (self.root / "summary.json").write_text(json.dumps({"complete": False}))
        with self.assertRaisesRegex(RuntimeError, "complete validated"):
            runner.validate_reference(self.args, self.manifest)


class TimingEvidenceTests(unittest.TestCase):
    def check(self, rows, mode="candidate", processors=4):
        return runner.timing_error(rows, "balanced_4", mode, 3, processors, 16)

    def test_complete_timing_block_accepted(self):
        self.assertIsNone(self.check(timing_rows()))
        self.assertIsNone(self.check(timing_rows(mode="original"), mode="original"))

    def test_untimed_controls_accepted(self):
        for mode, processors in (("default", 4), ("tight", 1)):
            with self.subTest(mode=mode):
                self.assertIsNone(self.check(timing_rows(mode=mode, processors=processors), mode, processors))

    def test_exact_repetitions_and_warmup_required(self):
        rows = timing_rows()
        variants = [[], rows[:-1], rows[1:], rows + [rows[-1]], rows[:-1] + [rows[-2]]]
        for altered in variants:
            with self.subTest(rows=altered):
                self.assertEqual(self.check(altered), "invalid_timing_rows")

    def test_exact_specification_required(self):
        for key, value in (("design", "chain_2"), ("mode", "original"), ("processors", "1"), ("control_count", "64")):
            with self.subTest(key=key):
                rows = timing_rows()
                rows[-1][key] = value
                self.assertEqual(self.check(rows), "invalid_timing_metadata")

    def test_rc_and_convergence_required(self):
        for key, value in (("rc", "430"), ("converged", "0"), ("converged", ".")):
            with self.subTest(key=key):
                rows = timing_rows()
                rows[-1][key] = value
                self.assertEqual(self.check(rows), "fit_failed")

    def test_positive_finite_timings_required(self):
        for value in (".", "nan", "inf", "-inf", "0", "-1", "", None):
            with self.subTest(seconds=value):
                rows = timing_rows()
                rows[-1]["seconds"] = value
                self.assertEqual(self.check(rows), "invalid_timing_seconds")

    def test_accelerated_attempt_required_for_every_fit(self):
        rows = timing_rows()
        rows[0]["schur_calls"] = "0"
        self.assertEqual(self.check(rows), "schur_inactive")
        rows[0]["schur_col_calls"] = "1"
        self.assertIsNone(self.check(rows))

    def test_diagnostics_cannot_be_missing_or_invalid(self):
        for value in (".", "nan", "-1", "1.5", None):
            with self.subTest(calls=value):
                rows = timing_rows()
                rows[-1]["schur_calls"] = value
                self.assertEqual(self.check(rows), "invalid_schur_diagnostics")

    def test_default_cannot_attempt_schur(self):
        rows = timing_rows(mode="default")
        rows[0]["schur_col_calls"] = "1"
        self.assertEqual(self.check(rows, mode="default"), "unexpected_schur_attempt")

    def test_missing_tight_control_rejected(self):
        self.assertEqual(self.check([], mode="tight", processors=1), "invalid_timing_rows")

    def test_timeout_diagnosis_is_preserved(self):
        result = {"status": "fit_timeout"}
        runner.check_fit(result, [], "balanced_4", "candidate", 3, 4, 16)
        self.assertEqual(result, {"status": "fit_timeout", "timing_error": "invalid_timing_rows"})


class ComparisonEvidenceTests(unittest.TestCase):
    def setUp(self):
        self.plan = {"repetitions": 3, "tight_cases": ["balanced_4"], "default_cases": ["balanced_4"]}
        self.rows = [{"case": "balanced_4", "mode": mode, "rep": rep, "reference": reference, "status": "ok", "passed": 1}
                     for reference in ("original", "tight") for mode in ("original", "candidate") for rep in range(1, 4)]
        self.rows.append({"case": "balanced_4", "mode": "default", "rep": 1, "reference": "original", "status": "ok", "passed": 1})

    def test_exact_comparison_set_required(self):
        self.assertIsNone(runner.comparison_error(self.rows, "balanced_4", self.plan))
        for altered in ([], self.rows[:-1], self.rows + [self.rows[-1]], self.rows[:-1] + [self.rows[-2]]):
            with self.subTest(count=len(altered)):
                self.assertEqual(runner.comparison_error(altered, "balanced_4", self.plan), "incomplete_comparisons")

    def test_failed_comparison_rejected(self):
        for key, value in (("passed", 0), ("status", "script_failed")):
            with self.subTest(key=key):
                rows = copy.deepcopy(self.rows)
                rows[-1][key] = value
                self.assertEqual(runner.comparison_error(rows, "balanced_4", self.plan), "failed_comparisons")

    def test_comparison_csv_requires_exactly_one_complete_finite_row(self):
        with tempfile.TemporaryDirectory() as temporary:
            output = Path(temporary)
            for mode in ("original", "candidate"):
                folder = output / "fits" / ("balanced_4_" + mode)
                folder.mkdir(parents=True)
                for suffix in (".ster", ".dta"):
                    (folder / ("fit_1" + suffix)).write_bytes(b"artifact")
            args = SimpleNamespace(output=output, repetitions=1)
            metrics = comparison_metrics()
            variants = [("ok", [metrics]), ("invalid_comparison_rows", []),
                        ("invalid_comparison_rows", [metrics, metrics]),
                        ("invalid_comparison_rows", [{"passed": "1"}]),
                        ("accuracy_failed", [dict(metrics, mu_scaled="nan")]),
                        ("accuracy_failed", [dict(metrics, passed="0")])]
            for expected, rows in variants:
                with self.subTest(expected=expected, rows=rows):
                    def fake_run(*unused):
                        write_csv(output / "comparisons/balanced_4_candidate_1_vs_original/comparison.csv", rows)
                        return {"status": "ok"}
                    with patch.object(runner, "run_stata", side_effect=fake_run):
                        results = runner.compare_fits(args, "balanced_4", "original", ("candidate",), 999)
                    self.assertEqual(results[0]["status"], expected)
            def failed_run(*unused):
                write_csv(output / "comparisons/balanced_4_candidate_1_vs_original/comparison.csv", [metrics])
                return {"status": "script_failed"}
            with patch.object(runner, "run_stata", side_effect=failed_run):
                results = runner.compare_fits(args, "balanced_4", "original", ("candidate",), 999)
            self.assertEqual(results[0]["status"], "script_failed")


class ProcessOutcomeTests(unittest.TestCase):
    def test_success_requires_zero_exit_and_marker(self):
        for returncode, marker, expected in ((0, True, "ok"), (7, True, "script_failed"),
                                             (0, False, "script_failed"), (7, False, "script_failed")):
            with self.subTest(returncode=returncode, marker=marker), tempfile.TemporaryDirectory() as temporary:
                directory = Path(temporary)
                (directory / "run.log").write_text("license banner\nREDBENCH_START\n" + ("SUCCESS" if marker else "failure"))
                process = SimpleNamespace(returncode=returncode, poll=lambda: returncode)
                args = SimpleNamespace(stata=Path("/fake/stata"), fit_timeout=1)
                with patch.object(runner.subprocess, "Popen", return_value=process), patch("sys.stdout", new_callable=io.StringIO):
                    result = runner.run_stata(args, "fit.do", [], directory, time.monotonic() + 60, "SUCCESS")
                self.assertEqual(result["status"], expected)
                self.assertNotIn("license banner", (directory / "run.log").read_text())


class RunnerIntegrationTests(unittest.TestCase):
    def test_only_complete_validated_evidence_gets_speedup(self):
        for damage in (None, "duplicate_timing", "missing_tight_row", "missing_comparison", "failed_process", "changed_fixture", "changed_executable"):
            with self.subTest(damage=damage), tempfile.TemporaryDirectory() as temporary:
                root = Path(temporary)
                original, candidate = root / "original", root / "candidate"
                for source in (original, candidate):
                    (source / "src").mkdir(parents=True)
                    (source / "src/ppmlhdfe.ado").write_text("source")
                stata = root / "stata"
                stata.write_bytes(b"executable")
                args = SimpleNamespace(original=original, candidate=candidate, output=root / "results", stata=stata,
                    n=10000, repetitions=3, processors=4, controls=16, precision="tight", reference_from=None,
                    dependency_dir=None, candidate_revision=None, candidate_label="candidate", first_source=None,
                    cases=["balanced_4"], fit_timeout=60, budget=600, dry_run=False)
                def fake_run(args, script, script_args, directory, deadline, marker):
                    directory.mkdir(parents=True, exist_ok=True)
                    if script == "environment.do" and directory.name == "environment_after":
                        if damage == "changed_fixture":
                            (args.output / "fixtures/balanced_4.dta").write_bytes(b"changed fixture")
                        if damage == "changed_executable":
                            args.stata.write_bytes(b"changed executable")
                    if script == "fixtures.do":
                        (directory / "balanced_4.dta").write_bytes(b"fixture")
                    elif script == "fit.do":
                        mode = script_args[3]
                        rows = timing_rows(mode=mode, processors=1 if mode == "tight" else 4)
                        if damage == "duplicate_timing" and mode == "candidate":
                            rows[-1] = rows[-2]
                        if damage == "missing_tight_row" and mode == "tight":
                            rows = []
                        write_csv(directory / "timings.csv", rows)
                        for rep in range(1, (1 if mode == "tight" else 3) + 1):
                            for suffix in (".ster", ".dta"):
                                (directory / (f"fit_{rep}" + suffix)).write_bytes(b"fit")
                        if damage == "failed_process" and mode == "candidate":
                            return {"status": "script_failed"}
                    elif script == "compare.do":
                        rows = [] if damage == "missing_comparison" else [comparison_metrics()]
                        write_csv(directory / "comparison.csv", rows)
                    return {"status": "ok"}
                environment = {"runtime": {"processors": "4", "stata_version": "19"}, "dependencies": {}}
                with patch.object(runner, "arguments", return_value=args), patch.object(runner, "verify_original", return_value={}), \
                        patch.object(runner, "read_environment", return_value=environment), \
                        patch.object(runner, "run_stata", side_effect=fake_run), patch("sys.stdout", new_callable=io.StringIO):
                    status = runner.main()
                summary = json.loads((args.output / "summary.json").read_text())
                self.assertEqual(status, 0 if damage is None else 1)
                self.assertEqual(summary["complete"], damage is None)
                self.assertEqual("speedup" in summary["cases"]["balanced_4"], damage is None)


if __name__ == "__main__":
    unittest.main()
