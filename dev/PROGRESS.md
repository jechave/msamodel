# Checklist — work item in progress

Checklist for the **work item currently being executed** (a slice of the current
release cycle — NOT a package version; the version is plain semver in
`DESCRIPTION`). Rewritten from that item's detailed plan when it starts; dormant
between items. The durable roadmap is `dev/plan.md`; the append-only history is
`dev/LOG.md`.

Tick items as substeps finish, and add a one-line dated entry to `dev/LOG.md`.

---

## IN FLIGHT — Slice 1: σ correctness fix in the MSA likelihood

`calculate_loglik_msa` profiles σ out but used `sd(residuals)` (divisor `n−1`); the
profiled-Gaussian MLE is `sqrt(mean(residuals^2))` (divisor `n`). Correctness fix.
Verified 2026-06-24: argmax `(a1,a2)` unchanged (= least-squares); only σ-derived
quantities move. First of three 0.3.0 fit slices (see `dev/plan.md` §v0.3 (B)).

- [x] Verify the math (profile MLE = `sqrt(mean(r²))`; argmax-invariance vs divisor)
- [x] Apply one-line fix in `R/objective.R:56`
- [x] Deterministic OLD-vs-NEW grid check on znb — argmax `(a1,a2)` IDENTICAL (0.5,31);
      logLik shift is a constant +0.00111 across the whole surface (Spearman cor = 1)
- [x] Seeded MCMC OLD-vs-NEW on znb — posterior **bit-identical** (samples identical,
      max summary diff = 0). Constant logLik shift ⇒ acceptance ratio unchanged.
- [x] Posterior did NOT move ⇒ no PAUSE needed; nothing downstream stale
- [x] Updated 3 tests that froze the *absolute* logLik value (now corrected, not
      drift): `test-msa-evaluate.R` literal → −184.3230779142; `test-contract.R`
      manual route `sd`→`sqrt(mean(res^2))`; `profile-invariance` snapshot accepted
- [x] `devtools::test()` green (83/0F); `devtools::check()` at v0.1 baseline (0E/1W/2N)
- [ ] Commit slice 1 on its own (after user review), before starting the ML arm

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
