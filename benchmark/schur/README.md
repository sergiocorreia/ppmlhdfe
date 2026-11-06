# Independent-source multiway ppmlhdfe benchmark

This bounded development benchmark compares a pristine upstream archive and a
candidate source tree in separate fresh Stata processes. It does not install or
overwrite any ado package. Python uses only the standard library (3.9+).
The runner uses the macOS/Linux Stata/MP batch interface. The do-files use
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
larger targeted study. N must be a multiple of 5,000 between 10,000 and
6,545,000. The upper bound keeps the N/200 firm levels within Stata's largest
nonmissing `int` value, 32,740, without changing fixture storage or draws. It is
a storage limit, not a claim that the largest fixture fits the default time or
memory budget. Both time limits must be finite and positive. There are N/10 workers, N/200 firms, ten observations per worker and 20% movers.
Movers spend five observations at each endpoint. The balanced four-way design
instead uses four independent factors with N/200 levels each.
For the default control counts at N=100,000, allow `--budget 1200 --fit-timeout 180`: the tight pristine
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
| balanced_4 | Four independent categorical effects | Poisson; no estimated covariates by default, exposing setup overhead |

By default, the first seven fit four Gaussian controls. The cases are deliberately varied
robustness and performance probes; comparisons across rows do not isolate the
causal effects of topology, weights, outcome type or FE count. Original and
candidate use the identical saved fixture within each row.

By default, each source/case process uses one processor, one warmup, then three
timed commands, all saving the FE sum. Both sources use `tol(1e-10)` except
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
plus an elapsed budget shared by the Stata stages, including validation. The
budget starts after initial source/provenance checks; synchronous hashing, copying
and final report generation can add wall time. It is not a hard whole-process
wall-time limit. The runner preserves completed timing rows when a later
fit times out. Failures, incomplete runs and numerical disagreements remain in
the output; they never receive a validated speedup. Stata license banners are
removed from archived logs. `manifest.json` records source, fixture and harness
hashes; the logs record Stata and dependency versions/resolved source paths.

This harness is bounded synthetic development evidence, not an exhaustive
matrix, empirical validation or release qualification. Small dummy-variable
projection oracles, deliberate solver-failure tests and postestimation tests
belong in the separate correctness suite.

## Processor and covariate scaling

The defaults use one Stata processor and four controls, except balanced_4,
which estimates no controls. `--processors P` requires exactly P actual Stata
processors, and `--controls K` estimates exactly K controls in every selected
design, including balanced_4. The supported control range is 0 through 64.
The preflight and each fit reject an unavailable processor count.

The additional x5 through x64 are independent zero-effect Gaussian controls,
drawn after every original outcome, weight and FE draw under seed 183271. The
original stream uses seed 98241. Increasing K therefore preserves the original
variables. The focused Stata check accepts an original and an extended fixture:

```stata
do benchmark/schur/verify_fixture_prefix.do original.dta extended-k64.dta
```

Use `--precision tight` for comparisons across processors or K. This applies
`tol(1e-12) use_exact_partial(1) start_inner_tol(1e-12) itol(1e-12)` to both
sources. Every scaling cell additionally compares all timed fits to an untimed
pristine reference fitted with one processor and those same options. Numerical
gates remain those in `compare.do`.

A completed run may supply its exact saved fixture and one-core tight reference
to another processor count via `--reference-from`. The parent's timed fits may
use several processors. Reuse requires a complete validated parent, matching
N/K/options, source and environment hashes, successful tight-reference status,
and a complete tight-reference timing record reporting actual P=1. Fixture
and saved-result bytes are rehashed before reuse.

For example, replace the capitalized paths and `CANDIDATE_SHA` with the source
commit being tested. `GIT_REPO` must contain both that commit and the pinned
upstream revision. Use new empty result directories and the same frozen harness
for both invocations.

```sh
python3 benchmark/schur/run.py --original ORIGINAL --candidate CANDIDATE \
  --git-repo GIT_REPO --candidate-revision CANDIDATE_SHA \
  --dependency-dir DEPENDENCIES --output RESULTS/balanced-k64-p4 \
  --cases balanced_4 --n 100000 --controls 64 --processors 4 --precision tight \
  --repetitions 3 --first-source original --fit-timeout 600 --budget 3600
python3 benchmark/schur/run.py --original ORIGINAL --candidate CANDIDATE \
  --git-repo GIT_REPO --candidate-revision CANDIDATE_SHA \
  --dependency-dir DEPENDENCIES --output RESULTS/balanced-k64-p8 \
  --cases balanced_4 --n 100000 --controls 64 --processors 8 --precision tight \
  --repetitions 3 --first-source candidate --fit-timeout 600 --budget 3600 \
  --reference-from RESULTS/balanced-k64-p4
```

`--dependency-dir` selects an isolated flat installation of the tested public
ado/Mata dependencies. Drivers use BASE, that directory, and the selected
ppmlhdfe source; installed SITE/PERSONAL/PLUS packages are excluded. The preflight
rejects a resolved dependency outside that directory. Prepare dependencies before
running the benchmark:

1. Choose version-pinned source URLs for ftools and reghdfe and record their
   revisions or release identifiers. Install their complete packages into an
   empty disposable PLUS directory using `sysdir set PLUS "SCRATCH_PLUS"`, then
   `net install ftools, from("FTOOLS_PINNED_URL")` and
   `net install reghdfe, from("REGHDFE_PINNED_URL")`. These are preparation
   commands; the benchmark runner does not install packages.
2. Copy the packages' public ado/Mata files and helpers into a separate flat
   directory. Include the `ms_*.ado` helpers and `reghdfe5_parse.ado`, as well as
   `reghdfe5.mata`; copying only the top-level ado and Mata files is insufficient.
   Include platform-compatible public plugin or compiled Mata files when the
   selected packages require them. Reject duplicate basenames with different
   contents. Keep Stata executables, licenses and unrelated packages out.
3. Verify the installation with the intended Stata version and a small isolated
   fit before freezing it. Any package compilation belongs in preparation.
   Preserve the source URLs, versions and complete directory with the result
   archive; the run manifest records every dependency file's hash.

For a quick check of the prepared installation, use the runner with
`--n 10000 --cases balanced_4 --controls 4 --repetitions 1 --precision tight`
and a fresh output directory. The preflight also checks required dependency
resolution. A failed check must be resolved before collecting timings.

Stata executable, source, fixture, harness and dependency changes invalidate a
run. System dynamic libraries are not included in these hashes; use a fixed
Stata installation.

`--candidate-label` records a readable source identity, and `--first-source`
allows source-block order to be balanced across cells. Keep results with
different N, K, P, precision or source/dependency hashes separate. Report exact
hardware and Stata versions, medians and ranges, and failed cells alongside
successful speedups. Three repetitions are descriptive development evidence.

Run `python3 benchmark/schur/test_runner.py` for license-free evidence-integrity
tests. They reject incomplete or duplicated fit/comparison records, invalid
metadata, nonfinite output, failed processes, and incompatible reused evidence.
A successful process requires both zero exit status and a completion marker.
These tests complement the Stata fixture, processor-count and numerical checks.
