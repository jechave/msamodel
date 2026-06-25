# Checklist — work item in progress

Checklist for the **work item currently being executed** (a slice of the current
release cycle — NOT a package version; the version is plain semver in
`DESCRIPTION`). Rewritten from that item's detailed plan when it starts; dormant
between items. The durable roadmap is `dev/plan.md`; the append-only history is
`dev/LOG.md`.

Tick items as substeps finish, and add a one-line dated entry to `dev/LOG.md`.

---

## DORMANT — AGQ inference loop 1 done (next: vignette redesign)

Two items shipped this cycle. Detail in `dev/LOG.md`; decisions in `dev/plan.md`
"Inference rework".

### par_fit/b leak cleanup — DONE (committed `b490943`)
Removed the dead `par_fit = c(a1, b)` field from both ML fitters so no internal
optimizer coordinate is returned; `cov` kept (log2(a2+1)-scale, documented).

### AGQ inference loop 1 — DONE
- [x] Provenance check (find-source.sh: no `tmp_src` source ⇒ new code)
- [x] `fit_lrmsd_i_msa_agq()` — adaptive Gauss–Hermite quadrature, Laplace-referenced
      to the ML fit; `"msa_agq"` object {a1,a2,sd,ci,nodes(a1,a2,log_weight),n_nodes,
      laplace(mu,cov),log_evidence}; deterministic; default `n_nodes=7`
- [x] `predict_lrmsd_i_agq()` — node-propagated credible band (all-quadrature, no
      Gaussian assumption, no seed); coarseness documented, NOT warned (option 1)
- [x] Sanity vs 61-grid ground truth (moments match; converge 3→5→7)
- [x] Roxygen + `document()`; scale convention (compute on t=log2(a2+1), report
      natural a2; cov on t) documented
- [x] `test-fit-agq.R` (40 assertions incl. change-of-measure check, frozen values)
- [x] intro vignette §5.2: AGQ added to the 3-method comparison + a banded-profile plot
- [x] `test()` / `check()` at baseline
- [x] Commit gate: user reviewed HTML, approved
- [x] Satellite pass; LOG entry; memory (scale convention + next-session)

**NEXT loop (separate, plan-mode item): vignette redesign / split.** The intro
vignette has overgrown into a site-branch deep-dive; the prediction sections are
crowded now that there are 3 methods. Rethink what vignettes exist (true intro vs
fitting/methods vs site/mode), split, and restructure. Plan at execution time.

---

## History — completed work (newest first)

### Slice 3 — DONE (3a fit-side naming + 3b mode arm)

### 3a — fit-side naming (function renames + lrmsd column rename)
- [x] Function renames: `fit_msa_ml`→`fit_lrmsd_i_msa_ml`,
  `run_mcmc_msa`→`fit_lrmsd_i_msa_mcmc`, `calculate_loglik_msa`→
  `calculate_loglik_lrmsd_i_msa` (R/, 6 test files, intro vignette). Arg
  `observed_data` kept; `run_msa_bayesian_analysis` name + arg + return-key kept.
- [x] Column rename `lrmsd_obs`→`lrmsd_i_obs`, derived `nlrmsd_obs`→`nlrmsd_i_obs`
  (R/, tests, data-doc, globals, vignette); `znb_profile.rda` regenerated (only that
  `.rda` changed).
- [x] `devtools::document()` (new `.Rd`, old 3 deleted, NAMESPACE updated).
- [x] `devtools::test()` — frozen values unchanged (pure rename verified).
- [x] Removed the scientific-consistency test `ML estimate sits at the seeded MCMC
  posterior mode` (test-fit-ml.R, 4000/1000 MCMC) — not a regression guard; the
  frozen drift-guards cover "output unchanged". Suite 116→110/0F; ~51s faster
  (161.7s→110.7s).
- [x] `devtools::check()` — v0.1 baseline (0E/1W/2N), re-run after the test removal.
- [x] Knit intro vignette (`setwd("vignettes")`) + preview → USER HTML approved → committed `e00f0ef`, pushed `0d0fc26..e00f0ef`.

### 3b — mode fit arm
- [x] Provenance re-check (find-source.sh, name + formula) — no `tmp_src` source ⇒ new code.
- [x] `calculate_loglik_lrmsd_n_msa` + `fit_lrmsd_n_msa_ml` (arg `observed_data`, column `lrmsd_n_obs`, fail-loud n-coverage).
- [x] Synthetic `znb_profile_n` fixture (truth = site ML fit; seed 2025, sd 0.30); data-doc `@source` synthetic; NSE globals. Only that `.rda` added.
- [x] `test-fit-ml-mode.R` (29 assertions: frozen + fixture determinism 1e-12). Recovers truth (0.458,42.30)→(0.449,40.82).
- [x] dr2n-analysis vignette fit section + SYNTHETIC callout → USER HTML approved.
- [x] document/test (139/0F) / check (0E/1W/2N); satellite (LOG/plan/memory).

**Recently completed (0.3.0 fit work, 2026-06-24 — newest first; full detail in
`dev/LOG.md`):**
- **Slice 2 — ML fit arm + MCMC-vs-ML vignette section** (`aa868a3`). New
  `fit_msa_ml()` (`R/fitting.R`): a maximum-likelihood **point** estimator of
  `(a1,a2)` over the slice-1-corrected `calculate_loglik_msa`, parallel to the MCMC.
  L-BFGS-B on `(a1, log2(a2+1))`, grid-max start, `optimHess` covariance,
  delta-method `se_a2`, fail-loud on a singular Hessian. 33-assertion drift-guard
  `test-fit-ml.R`. intro-vignette §5.2 compares MCMC fit / ML fit / observed —
  **both arms predicted at a single `(a1,a2)`** (posterior mean vs point estimate)
  so the only difference is the estimate, plus a direct MCMC-vs-ML scatter
  (R²=0.9999). Provenance: no `tmp_src` source ⇒ new code. Suite 116/0F.
- **Slice 1 — σ correctness fix** (`b7251be`). `calculate_loglik_msa` used
  `sd(residuals)` (divisor `n−1`); fixed to the profiled-Gaussian MLE
  `sqrt(mean(residuals^2))` (divisor `n`). Verified: argmax `(a1,a2)` unchanged (=
  least-squares), logLik shift a constant +0.00111, seeded MCMC posterior
  bit-identical. 3 frozen-value tests updated.

**Recently completed (0.3.0 predict cycle, 2026-06-19 — newest first; full detail in
`dev/LOG.md`):**
- **Mode nested-models + dr2n vignette parity** (`b8c0424`). New
  `calculate_lrmsd_n_nested_models()` (mode counterpart of the site one; `n` +
  `lrmsd_n_mm/ms/ma/msa`, no `pdb_site`). `dr2n-analysis` vignette restructured to
  mirror the intro (profile → nested+decomposition → sweeps), y = `lrmsd_n`, all
  modes, Consistency section removed, mode decomposition uses `free_y`.
- **Pure-vector decomposition + nested-models fn + `lrmsd_i_*` rename** (`65e7b9f`,
  `eeabb17`). `calculate_msa_decomposition(mm,ms,ma,msa)` is now a pure 4-vector
  function (was tibble-with-fixed-columns); new
  `calculate_lrmsd_i_nested_models()`; nested-model cols renamed
  `lrmsd_*`→`lrmsd_i_*`. Intro vignette got a new §3 (nested models +
  decomposition at fixed (a1,a2)) with paper component colours.
- **Decomposition formula bugfix** (`2d12997`, `c7c0403`). `phi_*` was the wrong
  (Shapley) formula; corrected to the sequential M0→MM→MS→MSA form. See
  [[decomposition-not-shapley]].
- **`dr2*` naming convention — DONE 2026-06-18** (`dr2_<indices>`; details in
  [[dr2-naming-convention]] and earlier LOG entries).

Suite 83/0F; `check()` at v0.1 baseline (0E/1W/2N); `DESCRIPTION` `0.3.0.9000`.

**Next work item — defined step-by-step at execution time:** either the remaining
0.3.0 work (observed `dr2_i`/`dr2_n` profiles from homologous structures +
alignment; fit `dr2_n`; joint fit) OR the **0.4 motion arm** (`dh_ijm` → `dh_njm` +
`nh_njm`). The site/mode predict side is now complete and symmetric (forward map,
nested models, decomposition for both `i` and `n`). See `dev/plan.md` + the
`project_next_session` memory. **When the next item starts:** enter plan mode, read
the touched code, write the detailed plan, then rewrite this checklist from it.
