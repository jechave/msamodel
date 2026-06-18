# Checklist — current version in flight

Checklist for the version **currently being worked on**. Rewritten from that
version's detailed plan when the version starts; dormant between versions. The
durable roadmap is `dev/plan.md`; the append-only history is `dev/LOG.md`.

Tick items as substeps finish, and add a one-line dated entry to `dev/LOG.md`.

---

## Current version: v0.3a — mode-form structural divergence (`dr2_njm`)

Started 2026-06-18. First/smallest slice of v0.3. Adds the mode response axis
using `dr2` (no new quantity, no motion); predict-only. Detailed plan: the
approved v0.3a plan (raw=all-effects SPM, per-purpose preprocess, slow loop).
Design decisions recorded in the `dev/plan.md` v0.3a bullet + memory
`project-msamodel-v03-mode-axis` / `project-penm-mutscan-tiers`.

- [x] Step 0 — plan docs (plan.md v0.3a bullet, this checklist, LOG entry)
- [x] Step 1 — `generate_spm_data.R`: add `dr2n` + `mode` SPM list-columns;
      update `@return`
- [x] Step 2 — `msa_bayesian_data_preparation.R`: add `preprocess_spm_mode()`
      (builds `dr2nmat`, no `pdb_site` map)
- [x] Step 3 — `msa_model_evaluation.R`: add `calculate_dr2n_msa()` (predict-only,
      no loglik counterpart)
- [x] Step 4 — `msamodel-package.R` globalVariables + `document()` (14 → 16
      exports; `get_mode` stays `:::`, not imported)
- [x] Step 5 — regenerate `znb_spm` fixture; relax the one-time `tmp_src`
      validation to shared columns
- [x] Step 6 — tests: `mode == penm:::get_mode` assertion in test-spm-generate;
      new `test-msa-mode.R` (shape + Parseval site/mode cross-check)
- [x] Step 7 — new `dr2n-analysis` vignette (`.Rmd.orig` + knit); generalize
      `.githooks/pre-commit` to all `*.Rmd.orig` pairs
- [x] Step 8 — DESCRIPTION dev version bump + `NEWS.md` bullet
- [x] Final — `test()` (59 pass) + `check()` at baseline (0E/1W/2N) + spot-check
      prints `v0.3a OK`; site path loglik `-184.3241923285` unchanged

**v0.3a COMPLETE 2026-06-18** (not yet committed). `check()` 0E/1W/2N = baseline;
the data-size NOTE grew (znb_spm 12.5→23.6 MB from the wide `dr2n` column) but is
the same pre-existing finding, not new. NOT committed — awaiting user review.
