# Schur benchmarks

Run `suite.py` to compare pristine ppmlhdfe with the Schur implementation.
`suite.csv` is the single maintained list of benchmark cells. It contains the
22 configurations selected from the completed scaling study, plus six new
AKM-style mobility comparisons. Source order is explicit in each row.

| Platform | Design | Observations | Controls | Processors |
|---|---|---:|---:|---:|
| SCC | chain_2, communities_3, imbalanced_4 | 100,000 | 4, 16 | 1, 4 |
| SCC | balanced_4 | 100,000 | 4, 16, 64 | 1, 4 |
| Mac | imbalanced_4, balanced_4 | 100,000 | 64 | 4, 8 |
| SCC | akm_low, akm_medium, akm_high | 200,000 | 4 | 1, 4 |

The first 22 cells were chosen because they completed the earlier study. That
selection makes this a practical repeatable benchmark, rather than a survey of
all designs. Earlier failed experiments remain in their original results archive;
they are not additional maintained suites. New failures are retained in the report.
The six AKM cells study mobility at a fixed panel size.

## Run the suite

Use Python 3.9+ and Stata/MP. The runner uses only the Python standard library.
First inspect the cells without launching Stata:

```sh
python3 benchmark/schur/suite.py --platform scc --dry-run
python3 benchmark/schur/suite.py --platform mac --dry-run
```

Then run on the selected machine, using new output directories outside the checkout:

```sh
python3 benchmark/schur/suite.py --platform scc \
  --original /path/to/pristine-ppmlhdfe --candidate /path/to/this-checkout \
  --git-repo /path/to/this-checkout --candidate-revision CANDIDATE_SHA \
  --dependency-dir /path/to/pinned-dependencies \
  --stata /path/to/stata-mp --output /path/to/results
```

Use `--platform mac` for the four Mac cells. To run only the new mobility cases,
add `--cases akm_low,akm_medium,akm_high` to the SCC command. On SCC, invoke this
command within a compute allocation with at least four slots; `suite.py` does not
submit jobs. Keep processor pairs on the same host and avoid concurrent benchmarks.
The original tree is checked against upstream
`b85665d8f674e93c1e07446a9a2b29be7f797910`; the Git repository must contain that
revision and the requested candidate revision. Actual processor counts are checked,
so a four-core Stata license cannot silently run an eight-core cell.

Results appear under `results/scc/` or `results/mac/`: `REPORT.md` gives the
comparison, `summary.json` retains every planned row, and each cell has its timing,
numerical-comparison and saved-fit files. Platforms are reported separately.
`run.py` is the single-cell engine used by the suite; it also supports short focused
checks with fewer repetitions. There is no separate cache/prototype benchmark.

## Designs

The chain has adjacent-firm worker moves. The community design adds year effects,
sparse links between groups of firms and probability weights. The imbalanced
four-way design has worker, firm, year and product effects, frequency weights and
an offset. The balanced design has four independent categorical factors. These
100,000-row fixtures retain their original random-number streams; adding controls
uses a separate stream after all original draws. Extra controls have zero effects.

The AKM-style panel has **20,000 workers, 1,000 firms and ten years**, with 20
workers initially assigned to each firm. The absorbed effects are worker, firm
and year. In each of the nine annual transitions, a worker changes firms with
probability **2%, 10% or 30%**. A mover draws uniformly from the other 999 firms;
workers may move repeatedly. These are annual probabilities, not fractions of
workers who ever move. Initial employers, move/destination draws, covariates and
outcome shocks are shared across the three rates.

The positive wage-like outcome has conditional mean
`exp(worker effect + firm effect + year effect + X beta)`, with mean-one
lognormal noise. This exercises PPML on an AKM-style employer/employee panel;
it is not a benchmark of the linear log-wage AKM estimator. Fixture metadata and
logs record actual moves, workers who ever move, and firm-graph components,
including isolated firms. Low mobility need not yield a connected graph.

## Timing and numerical comparison

Every suite cell runs one warmup and three subsequent fits in a fresh process
for each source. Warmup times are recorded but excluded from reported medians.
The timer covers the ppmlhdfe call, including setup and both projection stages;
fixture generation, prediction and comparisons are outside it. Minimum, median
and maximum times are retained. Speedup is original median / candidate median.
Source blocks alternate according to `suite.csv`, not the order of CSV rows.

Both sources use:

```stata
tol(1e-12) use_exact_partial(1) start_inner_tol(1e-12) itol(1e-12)
```

These common tight settings also address fitted-mean sensitivity in original
ppmlhdfe. Every timed fit is compared with original at the same processor count
and with a pristine one-core reference. Imbalanced and low-mobility AKM cells
also check exact agreement with `use_schur(0)`.

A successful lower-processor cell supplies the exact saved fixture and one-core
reference to its paired cell. If its fit fails, an intact fixture can still be
copied after checking its hash and generation settings; the next cell computes
a fresh reference. Reference reuse itself still requires a fully validated parent.
A failure therefore does not automatically skip the other processor count.
Independent reference fits do not establish a shared-reference cross-core comparison.

Comparison requires matching sample, omitted variables, coefficient names and
separation/singleton counts, and finite positive predictions. The fixed numerical
gates are in `compare.do`: scaled coefficient/covariance error below 1e-7,
deviance error below 1e-9, maximum scaled fitted-mean error below 1e-4, and errors
in coefficients/SE units, covariance/reference-SE units, variances and standard
errors below 1e-4. Default-off controls require exact results. A ratio is reported
only after all required fits, comparisons and unchanged-input checks pass.

Each row sets a per-stage timeout and a shared elapsed budget for Stata stages.
The AKM cells use 180s and 900s; retained SCC cells use 1,800s and 7,200s;
Mac cells use 1,200s and 3,600s. Initial checks and synchronous copying/hashing/report
work add overhead. Limits are not hard whole-process deadlines. Completed rows
survive a later timeout, and failures remain visible without a validated ratio.

## Dependencies and tests

Prepare an isolated flat installation of version-pinned public ftools and reghdfe
packages. Install their complete packages into a disposable Stata PLUS directory,
then copy the public ado/Mata files and helpers to the benchmark dependency folder.
Include `ms_*.ado`, `reghdfe5_parse.ado` and `reghdfe5.mata`; copying only the main
ado files is insufficient. Include required platform-compatible compiled helpers.
Reject conflicting duplicate basenames. Freeze the directory after a small fit.
Drivers use BASE, this folder and the selected ppmlhdfe source, excluding installed
SITE/PERSONAL/PLUS packages. Do not include Stata executables or license files.
The manifest records source, harness, executable, fixture and dependency hashes.

```sh
python3 benchmark/schur/test_runner.py
python3 benchmark/schur/test_suite.py
```

The first tests check result validation, source/reference compatibility and timeout
behavior. The second tests check the maintained cell list, stable source ordering,
failed-parent behavior, and independent fixture reuse. No Stata is required.

In Stata, run the focused fixture checks with external output paths:

```stata
do benchmark/schur/verify_akm.do benchmark/schur /path/to/akm-fixture-test
do benchmark/schur/verify_fixture_prefix.do original.dta extended-k64.dta
```

The AKM check verifies panel dimensions, genuine job changes, realized mobility,
metadata and shared random draws. Scientific estimator tests remain under `test/`.
