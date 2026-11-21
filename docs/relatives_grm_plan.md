# Related donors: plan

Goal: get calibrated main and GxC p-values when some donors are related
(families, sib pairs, cryptic relatives), with the dynamic-context model
(`--dynamicCovarColList`) and the donor GxC engine.

This is new modelling, not a bug fix. Today CASTIE treats donors as
independent in every fast path, and the analysis recipe turns the GRM off.

## Current state (main 53d34da)

- **Step 1 solver.** In `src/Main.cpp` `getPCG1ofSigmaAndVector_multiV`:
  - Case 1 (the closed form per donor block) and Case 2 (Woodbury with the
    context matrix) assume independent donors. Neither reads `g_isGRM`.
  - Only the Case 3 PCG fallback (`getCrossprod_multiV_eMat`) uses a GRM. There
    it is the dense genotype GRM, and it *replaces* the iid donor intercept.
  - The solves and the AI score/trace can therefore use different covariance
    models. This is the likely cause of the `missing value` and `dot()`
    failures with `--useGRMtoFitNULL=TRUE` on the tutorial (not proven).
- **Sparse GRM.** `gen_sp_Sigma_multiV` (the `g_isStoreSigma` path) builds
  Sigma without the context-slope term. That is the wrong model whenever
  dynamic contexts are set.
- **Step 2.** `prepare_gxe_condG_donor` and `getMarkerPval_gxe_donor` in
  `src/SAIGE_test.cpp` use exact donor-level sufficient statistics and assume
  independent donors. The main test uses the Step 1 variance ratio.
- **`createSparseGRM.R` segfaults.** `setminMAFforGRM` writes through
  `ptr_gNULLGENOobj` before `setgenoNULL()` has created it.
  - Fix: add `if (isSetGeno) setgenoNULL()` at the top of `createSparseGRM()`.
  - The guard matters because the internal caller passes `isSetGeno = FALSE`
    with genotypes already loaded.

## Model

Sigma = τ0 W⁻¹ + Z (τ1 I + τg K) Zᵀ + context slopes (iid per donor, as now).

- `K` is the sparse GRM, thresholded at the relatedness cutoff, so it is
  block-diagonal by family (connected components). Unrelated donors are
  1×1 blocks.
- The iid term τ1 stays: it carries non-genetic donor effects.
- τg is identified only through related pairs. With none above the cutoff, `K`
  is the identity and the model must reduce to today's, with τg not fitted.
- Version 1 keeps the slopes iid. Kinship-structured slopes (heritable context
  response) come later, and only if the sib simulation shows that the GxC
  tests need them.

The key property is that Sigma stays block-diagonal by **family**, with small
blocks. Every place that loops over donor blocks today can loop over family
blocks instead, with a small dense solve per family.

## Steps

0. **Measure (in progress).** Calibration of the current no-GRM recipe on
   simulated sib pairs. Setup:
   - `castie-calibration/sim.py --sibs`: 75 full-sib pairs from gene-dropping,
     donor intercept and slopes correlated 0.5 within pairs;
   - the null scenario N1, run with the verdict rule from the calibration
     README.

   If main and GxC are both OK, document "prune or ignore" and stop here.

   **Result (2026-11-20).** Relatives inflate the current recipe. Type-I
   ratio at α = 0.01 (observed / nominal, 95% peak-bootstrap interval), 200
   null peaks each:

   | Test | Unrelated (150 donors) | 75 sib pairs |
   |---|---|---|
   | main | 0.93 [0.87–0.98] | 1.28 [1.19–1.38], up to 1.83 in the densest tier |
   | GxC ctx1 | 1.09 [1.01–1.17] | 1.41 [1.32–1.51] |
   | GxC stim | 0.36 [0.31–0.42] | 0.44 [0.38–0.50] |

   λ for main goes from 1.00 to 1.10. The GxC ctx1 inflation comes from the
   sib-correlated slopes, so kinship-structured slopes are probably needed,
   not optional.
1. **`createSparseGRM.R` fix.** Check on the simulated families that the sib
   entries are about 0.5 and the rest are about 0.
2. **Family index.** From the sparse GRM sample IDs, map donors to families.
   Order cells by family, then donor. This generalises `g_I_start_indices`.
3. **Step 1 covariance.**
   - Case 1 and Case 2 solves per family block.
   - AI score and trace for the new τg.
   - Variance-ratio estimation with the same Sigma.
   - Retire or fix the Case 3 GRM branch, so that only one covariance model
     exists.
4. **Step 2.**
   - The donor GxC engine caches family-block statistics instead of donor
     scalars.
   - Decide whether the main test keeps the variance ratio or uses the same
     exact family-block computation. Exact is preferred: it avoids a
     variance-ratio error that depends on each variant's sib genotype
     correlation.
5. **CLI.** Use the existing `--sparseGRMFile`, `--sparseGRMSampleIDFile` and
   `--relatednessCutoff`. Add no new flags unless one is needed.

## Validation

- **Regression.** Without a sparse GRM, outputs must be byte-identical to
  main on the tutorial and on the calibration data. The same holds with a GRM
  that has no related pairs.
- **Calibration with sibs.** Null scenarios with and without the GRM. The GRM
  run must be OK under the README verdict rule for main, `ge_*` and CCT.
- **Power.** The N3 and POS analogues with sibs. Power must not be lower than
  in the no-GRM run on the same data.
- **θ recovery.** τg should approach the simulated genetic intercept
  variance.

## Interim option

Drop one donor of each related pair before analysis (for example
`plink2 --king-cutoff 0.0884` for second degree). This needs no code change
but loses donors. Use it until step 4 passes validation.

## Open questions

- Dataset: number of donors, family structure (nuclear families, sib pairs,
  cryptic), genotyping (array or WGS), cell type.
- Is the context response heritable enough to need kinship-structured slopes?
