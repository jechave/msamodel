# Checklist — work item in progress

Checklist for the **work item currently being executed** (a slice of the current
release cycle — NOT a package version; the version is plain semver in
`DESCRIPTION`). Rewritten from that item's detailed plan when it starts; dormant
between items. The durable roadmap is `dev/plan.md`; the append-only history is
`dev/LOG.md`.

Tick items as substeps finish, and add a one-line dated entry to `dev/LOG.md`.

---

## IN FLIGHT — Slice 2: ML `lrmsd_i` fit arm + 3-way comparison vignette section

Second of three 0.3.0 fit slices (see `dev/plan.md` §v0.3 (B)). Add `fit_msa_ml()`:
a maximum-likelihood **point** estimator of `(a1, a2)` over the (slice-1-corrected)
shared likelihood `calculate_loglik_msa`, parallel to the MCMC — for speed and
toward `lm`/`gam` conventions. Agile scope: estimator only, NO S3 / `method=` flag,
a **separate** fn from MCMC. Plus an intro-vignette §5.2 comparing MCMC fit, ML fit,
and the observed profile.

Provenance checked 2026-06-24: NO `tmp_src` source for ML-fitting `(a1,a2)` (source
project fit only by MCMC / grid; the lone `optim` is the unmigrated rate arm) ⇒ new
code, not a migration.

- [x] find-source.sh provenance check (`optim`/`optimHess`/objective) — recorded above
- [x] Add `optim optimHess` to `@importFrom stats`; `devtools::document()`
- [x] Write `fit_msa_ml()` in `R/fitting.R` — L-BFGS-B on `(a1, log2(a2+1))`,
      deterministic grid-max start, `optimHess` cov, fail-loud on singular Hessian,
      list return (a1, a2, logLik, sigma_hat, cov, se_a1, se_a2, convergence, par_fit)
- [x] Deterministic grid-max vs `fit_msa_ml` cross-check on znb — same basin; optim
      logLik −138.546 ≥ grid −138.632 (refines off-grid). a1≈0.458 ~ grid 0.50
- [x] Seeded MCMC-mode vs ML cross-check on znb — ML (0.458, 42.30) within posterior
      95% CI; means (0.470, 42.6); ML SEs (0.122, 10.19) ≈ posterior sd (0.131, 11.3)
- [x] Tests `test-fit-ml.R` (33 assertions) — grid agreement, MCMC-mode proximity,
      frozen regression, cov/SE shape, fail-loud paths, contract parity
- [x] `devtools::test()` green (116/0F, +33); `devtools::check()` at baseline (0E/1W/2N)
- [ ] `devtools::install()`; add §5.2 to `msamodel-intro.Rmd.orig`; knit; preview
- [ ] **STOP — user reviews `dev/preview/msamodel-intro.html`** (commit-gate)
- [ ] After approval: commit code+tests+vignette together; push to main
- [ ] Reconcile (PROGRESS dormant, LOG entry + hashes, memory, git clean)

### Slice 1 — σ correctness fix (DONE, committed `b7251be`, pushed)
`calculate_loglik_msa` used `sd(residuals)` (divisor `n−1`); fixed to the
profiled-Gaussian MLE `sqrt(mean(residuals^2))` (divisor `n`). Verified: argmax
`(a1,a2)` unchanged (= least-squares), logLik shift a constant +0.00111, seeded MCMC
posterior bit-identical. 3 frozen-value tests updated. Suite 83/0F.

---

## (was DORMANT)

**Recently completed (0.3.0 cycle, 2026-06-19 — newest first; full detail in
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
