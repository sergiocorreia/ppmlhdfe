# Independent-source multiway ppmlhdfe benchmark

This bounded development benchmark compares a pristine upstream archive and a
candidate source tree in separate fresh Stata processes. It does not install or
overwrite any ado package. Python uses only the standard library (3.9+).
The runner uses the macOS/Linux Stata batch interface. The do-files use
Stata 13 version control for reproducible random-number generation and syntax;
this is a benchmark requirement, not a change to the package's supported versions.

```sh
python3 benchmark/schur/run.py \
  --original /path/to/pristine-ppmlhdfe \
  --candidate /path/to/candidate-ppmlhdfe \
  --git-repo /path/to/upstream-git-checkout \
  --output /path/outside/checkout/benchmark-results-30000 \
  --n 30000 --repetitions 3 --fit-timeout 120 --budget 600
```

Use a new empty output directory. `--git-repo` verifies every pristine `src/`
file against upstream `b85665d8f674e93c1e07446a9a2b29be7f797910`. Omitting it
records source hashes but cannot establish pristine-source provenance. The
candidate and original hashes are rechecked after the run. `--dry-run` prints
the plan without launching Stata. `--stata` or `STATA` overrides the executable.

The eight default fixtures use N=30,000; `--n 100000` prepares the corresponding
larger targeted study. N must be a multiple of 5,000 and at least 10,000. There
are N/10 workers, N/200 firms, ten observations per worker and 20% movers.
Movers spend five observations at each endpoint. The balanced four-way design
instead uses four independent factors with N/200 levels each.
For N=100,000, allow `--budget 1200 --fit-timeout 180`: the tight pristine
chain reference can take longer than the ordinary timed fits.

| Case | Fixed effects | Outcome and other features |
|---|---|---|
| random_2 | Worker, firm | Positive continuous; connected random mobility with ring backbone |
| chain_2 | Worker, firm | Positive continuous; adjacent-firm moves |
| disconnected_2 | Worker, firm | Positive continuous; disconnected 25-firm communities |
| random_3 | Worker, firm, time | Poisson; two deliberately all-zero worker groups |
| communities_3 | Worker, firm, time | Positive continuous; sparsely bridged communities; lognormal probability weights |
| redundant_3 | Worker, firm, region | Positive continuous; region nested in firm creates redundant retained effects |
| imbalanced_4 | Worker, firm, time, product | Poisson; 40 crossed products; frequency weights and offset |
| balanced_4 | Four independent categorical effects | Poisson; no estimated covariates, exposing setup overhead |

The first seven fit four Gaussian controls. The cases are deliberately varied
robustness and performance probes; comparisons across rows do not isolate the
causal effects of topology, weights, outcome type or FE count. Original and
candidate use the identical saved fixture within each row.

Each source/case process uses one processor, one warmup, then three timed
commands, all saving the FE sum. Both sources use `tol(1e-10)` except
random_3, communities_3 and imbalanced_4, which use
`tol(1e-12) use_exact_partial(1) start_inner_tol(1e-12) itol(1e-12)`.
These common stricter settings were chosen after development runs: the
original adaptive path itself exceeded the prespecified 1e-4 maximum scaled
fitted-mean gate against its tight reference in those cases (1.54e-3, 1.52e-4
and 8.56e-4, respectively). In random_3 the accelerated adaptive path also
exceeded that gate (1.54e-3), while the two fully tight paths agreed within
5.14e-9. This exposes finite-stopping sensitivity in saved FE predictions even
when coefficients and deviance agree closely. Failed adaptive results remain
in the development evidence; the gate was not relaxed. The manifest
and each fit log record the exact options; timings between different settings
must not be pooled. Timings exclude compilation,
data generation, prediction and validation. Source block order alternates
across cases: original/candidate, then candidate/original. Repetitions within
a block are consecutive, not individually interleaved. No other performance
workload should run concurrently. Minimum, median and maximum are retained;
the speedup is the ratio of medians, not a confidence interval.

Each fit stores its `.ster` result and row IDs, sample indicators and fitted
means in a `.dta` artifact. A third fresh Stata process compares every timed
fit to the first original fit. Ordinary comparisons and default controls run
before expensive tight references, preserving partial validation if a reference
times out. Final speedup gates still require all planned reference comparisons.
For chain_2, random_3, communities_3 and imbalanced_4,
both paths are also compared to an untimed pristine reference using
`tol(1e-12) use_exact_partial(1) start_inner_tol(1e-12) itol(1e-12)`.
Untimed candidate `use_schur(0)` controls in random_2, random_3 and imbalanced_4
require exact coefficient, covariance, deviance and prediction equality to
pristine results, together with the exact sample/omission/count checks.
Every accelerated fit must report a Schur attempt; default controls must report
none. Initial-collinearity and IRLS diagnostics are recorded separately.

Acceptance requires exact coefficient stripes/omission masks, sample,
separation/singleton counts and finite positive predictions. Prespecified
continuous thresholds are scaled coefficient/covariance error <1e-7, deviance
error <1e-9, maximum `abs(mu-mu_ref)/(1+abs(mu_ref))` <1e-4, coefficient error
in reference-SE units <1e-4, covariance error scaled by the corresponding
reference standard deviations <1e-4, diagonal variance relative error <1e-4,
and standard-error relative error <1e-4. Weighted relative L2 fitted-mean error and maximum log-mean difference
are supplementary diagnostics. These gates are fixed in `compare.do`.

The runner enforces a wall-time limit on each fit through a progress marker,
plus a total study budget. It preserves completed timing rows when a later
fit times out. Failures, incomplete runs and numerical disagreements remain in
the output; they never receive a validated speedup. Stata license banners are
removed from archived logs. `manifest.json` records source, fixture and harness
hashes; the logs record Stata and dependency versions/resolved source paths.

This harness is bounded synthetic development evidence, not an exhaustive
matrix, empirical validation or release qualification. Small dummy-variable
projection oracles, deliberate solver-failure tests and postestimation tests
belong in the separate correctness suite.
