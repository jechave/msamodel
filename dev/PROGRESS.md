# Progress — msamodel migration

Human-readable record of where the migration stands and what happened when.
The **normative spec is `dev/plan.md`**; the checklist below mirrors its §9
migration order. The "what's next" pointer lives in Claude's project memory
(loads automatically), not here — this file is history + at-a-glance status.

## §9 migration order

- [x] 1. Create package skeleton
- [x] 2. Copy LICENSE + DESCRIPTION
- [x] 3. Add `R/msamodel-package.R` (`@importFrom` directives)
- [ ] 4. Copy `R/` files in dependency order:
  - [x] 4.1 `utils.R`
  - [x] 4.2 `loess_compare.R`
  - [ ] 4.3 `pdb_utils.R`, `enm_setup.R`, `site_properties.R`
  - [ ] 4.4 `spm_generate.R`, `spm_preprocess.R`
  - [ ] 4.5 `msa_evaluate.R`
  - [ ] 4.6 `msa_mcmc.R`, `msa_model_comparison.R`
  - [ ] 4.7 `msa_decomposition_site.R`, `msa_decomposition_protein.R`
  - [ ] 4.8 `msa_allotment_site.R`, `msa_allotment_protein.R`
  - [ ] 4.9 `a1a2grid.R`
  - [ ] 4.10 `msa_workflow.R` (last)
- [ ] 5. `devtools::document()` — confirm NAMESPACE populated
- [ ] 6. Write + run `data-raw/prepare_znb_data.R` → `data/*.rda`
- [ ] 7. Write `R/data-doc.R` (dataset docs); re-`document()`
- [ ] 8. Write tests (9 files, plan §7)
- [ ] 9. Write vignette (`msamodel-intro.Rmd`)
- [ ] 10. `devtools::check()` clean

Pre-skeleton groundwork (not a §9 step, but done): source snapshot,
docs reconciliation, git init — see the log.

## Log

Newest first. One short entry per working session.

### 2026-06-04
- **§9 step 4.2 done:** copied `compare_loess_fits.R` → `R/loess_compare.R`
  (renamed per §2). `compare_loess_fits` + `compare_loess_rmse` exported with
  `@family loess comparison`; `plot_loess_comparison` kept `@noRd` (deferred to
  v0.3), inline `library(ggplot2)` stripped, every ggplot2 call qualified
  `ggplot2::` behind the existing `requireNamespace` guard (§4 rule 2/8).
  Dropped the orphaned `@param evaluation_points` (not a real arg — doc-only).
  `gridExtra` left out of Suggests and no `globalVariables()` shim added (both
  decided: minimal change, the only consumer is the deferred internal plot fn).
  `document()` wrote 2 `.Rd` (none for the `@noRd` fn) + NAMESPACE (now 7
  exports); `load_all()` clean; smoke test of both exported fns returns finite
  values with the expected result names. Next: §9 step 4.3 — `pdb_utils.R`,
  `enm_setup.R`, `site_properties.R`.
- **Toolchain decision (not a §9 step):** investigated whether to update the R
  stack first. Finding: R itself is old (4.2.1; latest 4.6.0) but the packages
  are current (dplyr 1.1.4, ggplot2 4.0.1, etc.). **Decision: do NOT upgrade R
  now** — 4.2.1 satisfies the plan (`R >= 3.5`); revisit after v0.1. When done,
  upgrade R **side-by-side via `rig`** (keep 4.2 for the MSA paper project),
  reinstall packages into the new library, rebuild `penm` from source — that's
  the safe path (single shared library; all 436 pkgs built under 4.2.x). The
  only risk to penm/paper is upgrading *R*, not packages.
  **Updated roxygen2 7.2.3 → 8.0.0** (compiles under R 4.2; `Depends: R >= 4.1`):
  resolves the version-mismatch warning. 8.0.0 records its version as the new
  `Config/roxygen2/version: 8.0.0` field and removed the old `RoxygenNote` line.
  8.0.0 also enforces **one-line `@importFrom`** — flattened the wrapped dplyr/
  stats/penm directives in `msamodel-package.R`. `document()` + `load_all()`
  clean; 5 utils exports intact.
- **§9 step 4.1 done:** copied `utils.R` (`rmse_trend`, `rmse`, `r2`, `mrr`,
  `mrr_trend`). Stripped the `@requires` tags, added `@export` + `@family error
  metrics` to all 5, kept signatures/logic identical (no `library()` calls in
  source). `document()` wrote 5 `.Rd` + NAMESPACE (5 exports); `load_all()`
  clean. Bare `loess`/`predict`/`cor` resolve via the existing
  `@importFrom stats`. Next: §9 step 4.2, copy `loess_compare.R` (the ggplot2
  `requireNamespace` exception, rule 2/8).
- **§9 step 3 done:** added `R/msamodel-package.R` with the package-wide
  `@importFrom` directives (Option A namespacing — imports declared once;
  copied function files use bare calls). Inventory verified against the 16
  migrated source files; ggplot2 excluded (Suggests). `document()` populated
  NAMESPACE (incl. `importFrom(rlang,"!!")`, `importFrom(magrittr,"%>%")`, no
  exports yet); `load_all()` clean. Next: §9 step 4.1, copy `utils.R`.
- Set up source isolation: froze a snapshot of the source project to
  `tmp_src/` (build- and git-ignored), removed the live source from the
  session's additional directories, repointed all docs at `tmp_src/`.
- Wrote `CLAUDE.md`.
- Reviewed `dev/plan.md` against the actual source and reconciled it: fixed the
  `load_protein` test, made `plot_loess_comparison` internal (ggplot/Suggests),
  added a rule to strip non-standard `@requires` roxygen tags, corrected the
  `library()` inventory, added the `R/data-doc.R` step, folded the resolved
  pre-flight questions into §13.
- Reorganized: `dev/plan.md` is now the single normative spec; deleted
  `dev/RESUME.md`; added this `dev/PROGRESS.md`.
- `git init` on `main`. **Initial commit done** (`020ac2b`): skeleton +
  `dev/plan.md`, `dev/PROGRESS.md`, `CLAUDE.md`. `tmp_src/` and
  `.claude/settings.local.json` excluded.
- **Remote set up:** public repo `github.com/jechave/msamodel` (SSH origin),
  `main` pushed and tracking `origin/main`. Verified `tmp_src/` not pushed.

### 2026-04-22
- Created the package skeleton: `DESCRIPTION` (plan §3), `LICENSE` + `LICENSE.md`
  (MIT), empty roxygen-managed `NAMESPACE`, `.Rbuildignore`. `R/` left empty.
- Moved `plan.md` into `dev/`.
- Source-verification pass; answered the pre-flight questions.
