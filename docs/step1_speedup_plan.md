# Step 1 speed-up plan

Step 1 (`step1_fitNULLGLMM_qtl.R`) takes most of a CASTIE run, and much of that time is not model
fitting. This page records what was measured, why it is slow, and the three changes planned, in order.

## Benchmark

The data is a caQTL run on SHARE-seq ATAC peaks: 69,947 cells and 260 donors, a Poisson count trait,
38 covariates and 2 binary contexts with `--perContextTau=TRUE`. There is no GRM: donor intercept and
context slopes come through `I_longl`. The benchmark set is the 297 chr22 peaks of Stage A. Runs used
1 thread per peak on a 14-vCPU node (Xeon E5-2620 v4), fork commit `39fbbe7`, with `--tol=0.00001
--maxiter=500`. Numbers are taken from 2026-11-24, after 277 of the 297 peaks had finished.

| | step 1 (s) |
|---|---|
| median | 81 |
| mean | 165 |
| 90th percentile | 441 |
| max | 812 |

- **Unconverged fits:** 22 of 277 fits (7.9%) end with `Status: FAILURE` after `maxiter = 500`. They
  average 681 s each and take 33% of all step 1 time.
- **Step 2:** about 15 s per peak.

Where a typical converged peak's time goes (`[TIMING]` lines, 67 s of step 1 plus 17 s of step 2):

| Phase | s |
|---|---|
| R start, packages, container | ~13 |
| phenotype read and data prep (merge, model matrix, eMat) | 11.6 |
| initial GLM | 1.5 |
| AI-REML (8 iterations; each one Get_Coef 0.1-0.2 s + AI step ~1.2 s) | 13.0 |
| post-processing and save of the `.rda` (51 MB gzip) | 8.2 |
| variance ratio (80 random markers) and a second save | 19.4 |
| step 2: R start, model reload, 753 variants | 16.7 |

## Why it is slow

1. **`tol` does two jobs.**
   - In the R loop (`glmmkin.ai_PCG_Rcpp_multiV_eMat`) it is the convergence test,
     `max(|tau - tau0| / (|tau| + |tau0| + tol)) < tol`.
   - Passed to `fitglmmaiRPCG_multiV_eMat` (`src/Main.cpp`, the `tauVec.elem(find(tauVec < tol ...)).zeros()`
     lines), it is also the cutoff below which a variance component is set to 0. The R loop then stops as
     soon as all components are 0.
   - The variances here are ~1e-3, so `tol` must stay far below that: with `--tol=0.02` every
     component is zeroed after 1-2 iterations. On one test peak the lead GxC p-value went from
     10^-2.46 to 10^-3.44, an anti-conservative test.
   - With `tol = 1e-5`, tau settles within ~20 iterations and then jitters at ~1e-3 relative, so the
     test is never met and the fit runs to `maxiter`. Example: `chr22_32474489_32474989`, tau
     (0.010568, 0.001929, 0.000531) to 4 digits from iteration ~20 to 500.
2. **The traces are estimated stochastically.** `GetTrace_multiV_eMat` (`src/Main.cpp`) estimates
   `tr(P K_k)` with 30 Rademacher vectors, adds 10 at a time until the CV is below `traceCVcutoff = 0.0025`
   (up to 80 here), and costs one Sigma^-1 solve plus kernel products per vector.
   - The seed is reset every call (`set_seed(200)`). Still, the number of vectors changes between
     iterations, and the arithmetic is single precision, so the AI update is not a smooth function of
     tau. This is the source of the jitter in point 1.
   - It is also most of the ~1.2 s of each AI step.
3. **Per-peak fixed cost.** Each peak:
   - starts R twice (steps 1 and 2);
   - rebuilds the same model matrix and contexts;
   - saves the 51 MB model twice;
   - reloads it in step 2;
   - estimates a variance ratio from 80 markers (0.98 on the typical peak).

   This is ~55 s of the ~84 s per typical peak.

## Plan
The sections below are the plan as written before the work; the outcomes are under Status.

### 1. Separate the zeroing cutoff from `tol` (merged)
- New option `--tauZeroCutoff` (branch `feat/tau-zero-cutoff`, `8beb378`). Variance components below it
  are set to 0; NA (the default) uses `tol`, so existing runs are unchanged.
- Only `fitNULLGLMM_multiV` and the two Poisson/binomial loops (`glmmkin.ai_PCG_Rcpp_multiV`,
  `_eMat`) change. `count_nb` is not touched.
- Use for small variances: `--tol=1e-3` (or `1e-4`) `--tauZeroCutoff=1e-5`.
- **Validation (done, see Status):**
  - 32 chr22 peaks, made up of the 18 unconverged at the time plus 14 converged spanning fit times and
    zero/nonzero tau.
  - Each is fitted at `--tol` 1e-3 and 1e-4 with `--tauZeroCutoff=1e-5` and compared against the
    baseline for status, iterations, theta, and step 2 `p.value` / `pval_ge_CCT`.
  - A back-compat check runs the new build with the original flags on 3 peaks.
  - Accept if all fits converge and, on baseline-converged peaks, |delta -log10 p| < 0.05 for main and GxC.
- **Expected:** removes the `maxiter` tail (33% of step 1 time on chr22) and the slow-converging
  fits (e.g. 242 iterations).

### 2. Exact traces for donor-block kernels (merged)
- Here `K_k = Z_k Z_k^T`, with `Z_1 = I_longl` (donor intercept) and `Z_{1+c} = diag(e_c) I_longl`
  (`eKernel`), so Sigma = diag(tau_0 / w) + sum_k tau_k K_k is block diagonal over donors.
- For donor d, let `U_d = [1, e_1, ..., e_C]` on its cells, `G = diag(tau_1, ..., tau_{1+C})` and
  `A_d = U_d^T W_d U_d` (W_d = diag(w / tau_0) on the donor's cells). Then
  `U_d^T Sigma_d^-1 U_d = A_d (I + G A_d)^-1`, which holds also when some `tau_k = 0`.
- With `P = Sigma^-1 - Sigma^-1 X cov X^T Sigma^-1`:

      tr(P K_k) = sum_d [A_d (I + G A_d)^-1]_kk  -  tr(M_k cov M_k^T),   M_k = Z_k^T Sigma^-1 X

  `M_k` is donor sums of `u_k * Sigma_iX` (already computed).
- Cost per iteration is O(n (1+C)^2 + n p + D p^2): milliseconds instead of up to 80 solves, and
  deterministic. This is the same idea as the step 2 donor engine (`docs/gxe_donor_aggregation.md`).
- **Steps:**
  1. Prototype in R on saved Stage A models. Exact vs `GetTrace_multiV_eMat` should agree within its CV.
  2. Implement in C++ behind a flag, used only when there is no GRM and the kernels are `I_longl`/`eKernel`.
  3. A/B test on the 32-peak set and then all 297 chr22 peaks: theta, step 2 p-values, time.
  4. With a deterministic AI update, check whether the tight `tol` converges on its own.
- **Expected:** AI-REML 13 s to ~2 s on a typical peak, no jitter, and fewer iterations.

### 3. Fit many phenotypes per R session (merged)
- One R process per task would build the design, contexts and donor structures once, then loop over
  peaks. Step 2 would run in the same session from the in-memory model, without saving it, and only
  one save would happen if the model is kept.
- This is the largest remaining cost (~55 s of ~84 s per typical peak), but also the largest change.
- **Cheaper pieces to measure first:**
  - fewer variance-ratio markers (`--numRandomMarkerforVarianceRatio`, `--ratioCVcutoff`);
  - saving the model once instead of twice.

## Status
All three changes are merged into `main` (PR #10, merge b26838e; image built from 02b124f).
- [x] Measurements and causes (this page)
- [x] 1: `--tauZeroCutoff` (8beb378). The caQTL runs use `--tol=1e-3 --tauZeroCutoff=1e-5`.
- [x] 2: `--exactTrace` (22b07a0), plus `--activeSetAI` and `--tolEta` (6c86447), which fit the null model
  to the constrained REML optimum.
- [x] 3: `extdata/step12_batch_qtl.R` (e731da6) runs steps 1 and 2 of many peaks in one R session.
  `feat/fast-step1` adds a donor-level variance ratio, in-memory model reuse and cheaper per-peak setup.

Validation for the merge:
- the tutorial (Steps 1-4) is byte-identical to fork commit 39fbbe7 with default flags;
- the 297 chr22 peaks have 0 failed, 0 unconverged and 0 step 1 fallbacks;
- the genome-wide run (14,616 peaks) has 0 unconverged and 0 fallbacks, and its 26 failures are peaks with no
  markers in the cis window.

Speed: about 4.2 s per peak for steps 1 and 2 on a quiet node, against 84 s for the typical peak above.
The genome-wide run took 4 h 41 min at 12 concurrent tasks.
