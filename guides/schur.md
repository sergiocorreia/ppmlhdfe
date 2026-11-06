# Optional categorical Schur projection

```stata
ppmlhdfe y x1 x2, absorb(worker firm year) use_schur(1)
```

`use_schur(1)` accelerates the projection used by the initial collinearity check
and IRLS. The default is `use_schur(0)`. This is an opt-in numerical strategy;
it does not change the model, separation algorithms, rank tests, or convergence
options. Speed depends on the design. Setup and extra checks can make easy or
balanced designs slower.

## Supported models

The projector accepts two or more categorical effects, including categorical
interactions. It eliminates the effect with the most levels in the current
estimation sample; ties go to the first effect in `absorb()`. All remaining
effects are solved jointly. More than two effects do not require separate
pairwise approximations.

One-effect models, heterogeneous slopes, pruning, and the internal RRE or
alpha-storage modes use ordinary absorption. Eligible models can also fall back
on individual projections. Weight-dependent quantities are recomputed on every
call, and grouping is rebuilt when IRLS starts again after separation changes
the sample. Initial guesses and separation-specific collinearity checks retain
the original implementation. Predictions and saved effects use existing code.

## Reduced system

Let `A` be the eliminated dummy block, `B` the concatenation of all other dummy
blocks, and `W` the diagonal matrix of positive observation weights. For any
response/covariate matrix `Y`, define

```text
M_A Y = Y - A (A' W A)^(-1) A' W Y
S     = B' W M_A B
rhs   = B' W M_A Y
S gamma = rhs
residual = M_A (Y - B gamma).
```

The implementation scales retained coefficients by the square roots of their
weighted level counts and uses diagonally preconditioned conjugate gradients.
Its Jacobi diagonal sums `w_ab * (1 - w_ab / w_a)` within each retained level
and divides by that level's total weight `w_b`;
this avoids subtracting two large, nearly equal totals. Zero diagonal entries
are accepted only when their right-hand side is numerically zero.

Neither a dense dummy matrix nor a dense Schur matrix is constructed. With two
effects, multiplication uses coalesced cells and the normalized `I - C'C`
operator. With three or more effects, each multiplication predicts the sum of
retained effects, demeans within `A`, and aggregates to every retained block.
The system may be singular because of component constants or redundant effects;
only the projected residual is needed, not uniquely identified FE coefficients.

Eliminating the largest block minimizes the number of retained coefficients.
It does not necessarily minimize runtime: conditioning, observation scans,
weight variation, and setup costs also matter. Multiple right-hand sides share
the operator but retain independent CG step lengths and convergence flags.

## Numerical safeguards

Weights and input columns are scaled internally to limit overflow. A kernel
failure leaves its input unchanged. Acceptance requires finite output, fresh
reduced-system residuals, weighted group moments for *all* effects, and weighted
norm contraction. Exhausted iterations, invalid weights, or failed checks return
to ordinary absorption on the original input.

The ordinary reghdfe projection always polishes a successful candidate. Its
collinearity threshold depends on each **original input column's** unweighted
norm. The wrapper preserves those norms, the original means, and the `kept`
flags in both the initial check and IRLS. It falls back when a residual/original
squared-norm ratio lies within a factor of 100 of
reghdfe's cutoff `min(1e-6, tolerance/10)`.

A condition check on the small normalized residual Gram matrix also falls back
for zero columns or condition numbers at least `1e10`, protecting the separate
`invsym()` rank decision for nearly dependent regressors. The zero-column check
also catches false zeroing during polishing: a weighted projection may increase
the unweighted norm, making the intermediate norm an unsuitable omission
reference. Zero or unrepresentable original squared norms also trigger fallback.
These are conservative
fallback guards; they do not replace the original omission thresholds. This
check costs O(N K^2 + K^3) for K projected columns and can limit gains with many
controls. Floating-point results need not be bitwise identical when the option
is enabled. Changing the absorption trajectory can change finite-tolerance
IRLS stopping and saved-FE predictions, even when slopes and standard errors
agree closely. For demanding prediction comparisons, use the same tighter
settings in both paths; the benchmark documents those settings and retains
failed ordinary-setting checks. Extremely ill-conditioned models retain the
limitations of the reference absorption routine.

The outer IRLS stopping rule is unchanged. It includes an absolute floor in
the deviance denominator, so scaling all probability weights by a very small
constant can change finite-tolerance stopping in either implementation. The
tiny-weight estimator test checks agreement between the two paths under stricter
settings; it does not establish full-estimator invariance to weight scaling.
Both paths can still differ from an accurately fitted normal-scale reference.
That inherited convergence limitation is outside this projection change. The
projection tests separately verify common-scale invariance against an independent
dense-dummy oracle, and distinguish small group moments from actual residual
accuracy under extreme relative weights.

## Diagnostics

| Result | Meaning |
|---|---|
| `e(schur_calls)` | Eligible projection attempts across IRLS invocations. |
| `e(schur_iterations)` | Total reduced CG iterations, including rejected attempts. |
| `e(schur_fallbacks)` | IRLS attempts rejected in favor of ordinary absorption. |
| `e(schur_col_calls)` | Eligible initial collinearity projection attempts (0 or 1). |
| `e(schur_col_iters)` | Reduced CG iterations in that initial projection. |
| `e(schur_col_fails)` | Rejected initial attempts (0 or 1). |

Unsupported models and default-off fits report zero attempts. `e(ic2)` continues
to count reghdfe iterations; it does not include reduced CG work. There are no
production timing globals or separate switches for the two call sites.

## Tests and benchmarks

From the checkout root, with ppmlhdfe's dependencies installed:

```stata
do test/schur.do .
do test/schur_collinear.do .
do test/schur_projection.do .
do test/schur_weighted_growth.do .
```

The tests cover end-to-end estimation, exact and near collinearity, independent
explicit-dummy projection oracles, redundant/disconnected effects, row and FE
ordering, changing weights, extreme scales, transactional failure, separation,
saved effects, and unsupported models. They are also registered in
`test/testall.do`. The existing full suite reinstalls ppmlhdfe; use a disposable
Stata PLUS directory if the installed package must be preserved.

See [benchmark/schur/README.md](../benchmark/schur/README.md) for a reproducible comparison
against a **pristine upstream source tree** in fresh Stata processes. It reports
both gains and regressions, all repetitions, sample/omission agreement,
coefficient and covariance errors, fitted-mean errors, tighter reference fits,
source hashes, and unsuccessful runs. Timings exclude compilation through an
untimed warmup and include estimator setup and both projection stages.
