"""Tests for the maintained cell list, execution plan and fixture reuse."""

import copy
import csv
import hashlib
import io
import json
from pathlib import Path
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch

import run as engine
import suite


class SuiteTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)

    def test_catalog_is_requested_successful_subset_plus_akm(self):
        cells = suite.load_cells()
        actual = {
            (c["platform"], c["design"], c["controls"], c["processors"]) for c in cells
        }
        expected = {
            ("scc", d, k, p)
            for d in ("chain_2", "communities_3", "imbalanced_4")
            for k in (4, 16)
            for p in (1, 4)
        }
        expected |= {("scc", "balanced_4", k, p) for k in (4, 16, 64) for p in (1, 4)}
        expected |= {
            ("mac", d, 64, p) for d in ("imbalanced_4", "balanced_4") for p in (4, 8)
        }
        expected |= {
            ("scc", d, 4, p)
            for d in ("akm_low", "akm_medium", "akm_high")
            for p in (1, 4)
        }
        self.assertEqual(actual, expected)
        self.assertEqual(len(cells), 28)

    def write_cells(self, cells):
        p = self.root / "suite.csv"
        with p.open("w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=sorted(suite.FIELDS))
            w.writeheader()
            w.writerows(cells)
        return p

    def test_reordering_catalog_preserves_plan_and_source_order(self):
        cells = suite.load_cells()
        self.assertEqual(
            suite.load_cells(self.write_cells(list(reversed(cells)))), cells
        )

    def test_duplicate_or_invalid_cells_rejected(self):
        cells = suite.load_cells()
        with self.assertRaises(ValueError):
            suite.load_cells(self.write_cells(cells + [cells[0]]))
        for key, value in [
            ("processors", 8),
            ("observations", 30000),
            ("annual_mobility", ".40"),
            ("first_source", "unknown"),
            ("budget", 0),
        ]:
            bad = copy.deepcopy(cells)
            cell = next(c for c in bad if c["platform"] == "scc")
            cell[key] = value
            with self.subTest(key=key), self.assertRaises(ValueError):
                suite.load_cells(self.write_cells(bad))

    def test_complete_parent_uses_reference(self):
        (self.root / "summary.json").write_text('{"complete":true}')
        self.assertEqual(
            suite.reuse_option(self.root), ["--reference-from", str(self.root)]
        )

    def test_failed_parent_uses_fixture_and_requests_fresh_reference(self):
        (self.root / "summary.json").write_text('{"complete":false}')
        (self.root / "manifest.json").write_text(
            '{"fixture_sha256":{"akm_low.dta":"hash"}}'
        )
        reuse = suite.reuse_option(self.root)
        self.assertEqual(reuse, ["--fixture-from", str(self.root)])
        cell = next(
            c
            for c in suite.load_cells()
            if c["design"] == "akm_low" and c["processors"] == 4
        )
        args = SimpleNamespace(
            original=Path("o"),
            candidate=Path("c"),
            git_repo=Path("g"),
            candidate_revision="rev",
            dependency_dir=Path("d"),
            stata=None,
        )
        cmd = suite.command(args, cell, self.root / "output", reuse)
        self.assertIn("--fixture-from", cmd)
        self.assertNotIn("--reference-from", cmd)
        self.assertEqual(cmd[cmd.index("--processors") + 1], "4")

    def test_early_parent_failure_does_not_suppress_next_cell(self):
        self.assertEqual(suite.reuse_option(self.root), [])
        self.assertEqual(suite.reuse_option(None), [])

    def test_pending_and_failed_cells_cannot_complete_report(self):
        cells = [dict(c, status="pending") for c in suite.load_cells()[:2]]
        self.assertFalse(suite.save_report(self.root, cells))
        cells[0].update(status="failed")
        self.assertFalse(suite.save_report(self.root, cells))
        d = json.loads((self.root / "summary.json").read_text())
        self.assertEqual(len(d["cells"]), 2)
        self.assertNotIn("×", (self.root / "REPORT.md").read_text())

    def test_bad_child_summary_does_not_skip_next_processor(self):
        cells = [
            c
            for c in suite.load_cells()
            if c["platform"] == "scc"
            and c["design"] == "balanced_4"
            and c["controls"] == 4
        ]
        output = self.root / "output"
        argv = [
            "suite.py",
            "--platform",
            "scc",
            "--original",
            str(self.root),
            "--candidate",
            str(self.root),
            "--git-repo",
            str(self.root),
            "--candidate-revision",
            "test",
            "--dependency-dir",
            str(self.root),
            "--output",
            str(output),
        ]
        commands = []

        def child(cmd, **kwargs):
            commands.append(cmd)
            dest = Path(cmd[cmd.index("--output") + 1])
            if len(commands) == 1:
                (dest / "summary.json").write_text("{truncated")
                (dest / "manifest.json").write_text("{also truncated")
            else:
                m = {
                    "validated": True,
                    "speedup": 2,
                    "original": {"median": 2},
                    "candidate": {"median": 1},
                }
                (dest / "summary.json").write_text(
                    json.dumps({"complete": True, "cases": {"balanced_4": m}})
                )
            return SimpleNamespace(returncode=0)

        with patch("sys.argv", argv), patch.object(
            suite, "load_cells", return_value=cells
        ), patch.object(suite.subprocess, "run", side_effect=child), patch(
            "sys.stdout", new_callable=io.StringIO
        ):
            self.assertEqual(suite.main(), 1)
        self.assertEqual(len(commands), 2)
        self.assertNotIn("--reference-from", commands[1])
        self.assertNotIn("--fixture-from", commands[1])
        records = json.loads((output / "scc/summary.json").read_text())["cells"]
        self.assertEqual([r["status"] for r in records], ["failed", "passed"])
        self.assertIn("summary_error", records[0])
        self.assertTrue(records[1]["reuse_error"])

    def test_corrupt_summary_still_allows_independent_fixture_check(self):
        (self.root / "summary.json").write_text("{truncated")
        (self.root / "manifest.json").write_text(
            '{"fixture_sha256":{"akm_low.dta":"hash"}}'
        )
        self.assertEqual(
            suite.reuse_option(self.root), ["--fixture-from", str(self.root)]
        )


class FixtureReuseTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        (self.root / "fixtures").mkdir()
        self.fixture = self.root / "fixtures/akm_low.dta"
        self.fixture.write_bytes(b"fixture bytes")
        self.args = SimpleNamespace(cases=["akm_low"])
        self.manifest = {
            "plan": {"cases": ["akm_low"], "N": 200000, "control_count": 4},
            "harness_sha256": {"fixtures.do": "generator"},
            "fixture_sha256": {
                "akm_low.dta": hashlib.sha256(self.fixture.read_bytes()).hexdigest()
            },
        }
        self.save(self.manifest)

    def save(self, data):
        (self.root / "manifest.json").write_text(json.dumps(data))

    def test_fixture_only_reuse_does_not_require_successful_fit(self):
        (self.root / "summary.json").write_text('{"complete":false}')
        self.assertEqual(
            engine.validate_fixture(self.args, self.manifest, self.root), self.manifest
        )

    def test_changed_fixture_rejected(self):
        self.fixture.write_bytes(b"different bytes")
        with self.assertRaises(RuntimeError):
            engine.validate_fixture(self.args, self.manifest, self.root)

    def test_wrong_generator_dimensions_controls_or_case_rejected(self):
        for key, value in [
            ("N", 100000),
            ("control_count", 16),
            ("cases", ["akm_high"]),
        ]:
            d = copy.deepcopy(self.manifest)
            d["plan"][key] = value
            self.save(d)
            with self.subTest(key=key), self.assertRaises(RuntimeError):
                engine.validate_fixture(self.args, self.manifest, self.root)
        d = copy.deepcopy(self.manifest)
        d["harness_sha256"]["fixtures.do"] = "changed"
        self.save(d)
        with self.assertRaises(RuntimeError):
            engine.validate_fixture(self.args, self.manifest, self.root)


if __name__ == "__main__":
    unittest.main()
