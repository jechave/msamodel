# Log — msamodel migration

Append-only history of what was done and attempted (newest first), including
dead ends and decisions. The **normative spec is `dev/plan.md`**; the live
checklist is `dev/PROGRESS.md` (a 1:1 mirror of plan.md §9). This file is the
history those two don't keep.

One short entry per working session.

### 2026-06-04
- **Docs split:** `dev/PROGRESS.md` was doing two jobs (append-only log + a
  checklist meant to mirror plan.md §9); the re-scope made the checklist
  numbering diverge from the plan and lose meaning. Split into three files:
  `dev/plan.md` (normative, unchanged role; §9 step-4 re-labelled `4.1`–`4.7`),
  `dev/PROGRESS.md` (pure checklist, bare 1:1 mirror of §9), and this
  `dev/LOG.md` (history). Updated the `## Tracking progress` convention in
  `CLAUDE.md` and the cross-references.
- **v0.1 RE-SCOPED (major):** narrowed v0.1 to "compute divergence profiles + fit
  to data". Traced the real call graph from `run_msa_bayesian_analysis` + the grid
  entry points to decide scope. **Kept:** setup, SPM, profiles/eval, grid, MCMC
  fitting + the four model profiles, site-level Shapley decomposition (wired into
  the workflow). **Deferred to vNext:** the model-fit *assessment* layer
  (`model_comparison_functions.R` + its sole-consumer `loess_compare.R` + the
  `utils.R` error metrics), *allotment* (unsettled), protein-level decomposition,
  trajectory route, visualization. Decisions: keep site decomposition; defer
  assessment (NOT the profiles — those are kept); pull `utils.R` + `loess_compare.R`
  back out. **Reverted** R/utils.R + R/loess_compare.R + their man pages + 7
  exports (forward commit, no history rewrite); `document()`/`load_all()` clean,
  NAMESPACE now 0 exports, R/ = just `msamodel-package.R`. **Rewrote `dev/plan.md`
  §1/§2/§6/§7/§8/§9/§11/§12** to match (exports 39 → 18; tests 9 → 7; workflow
  edited to drop assessment calls → 7-elem return). Reconciled `CLAUDE.md` and
  the progress docs. The §9 step-4 sub-step list is the new dependency order.
  Next: resume migration at §9 step 4.1 (`pdb_utils.R`, `enm_setup.R`,
  `site_properties.R`) against the re-scoped plan.
- **§9 step 4.2 done** (later REVERTED in the re-scope above): copied
  `compare_loess_fits.R` → `R/loess_compare.R` (renamed per §2).
  `compare_loess_fits` + `compare_loess_rmse` exported with
  `@family loess comparison`; `plot_loess_comparison` kept `@noRd` (deferred to
  v0.3), inline `library(ggplot2)` stripped, every ggplot2 call qualified
  `ggplot2::` behind the existing `requireNamespace` guard (§4 rule 2/8).
  Dropped the orphaned `@param evaluation_points` (not a real arg — doc-only).
  `gridExtra` left out of Suggests and no `globalVariables()` shim added (both
  decided: minimal change, the only consumer is the deferred internal plot fn).
  `document()` wrote 2 `.Rd` (none for the `@noRd` fn) + NAMESPACE (then 7
  exports); `load_all()` clean; smoke test of both exported fns returned finite
  values with the expected result names.
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
  clean; 5 utils exports intact (those exports later reverted in the re-scope).
- **§9 step 4.1 done** (later REVERTED in the re-scope above): copied `utils.R`
  (`rmse_trend`, `rmse`, `r2`, `mrr`, `mrr_trend`). Stripped the `@requires`
  tags, added `@export` + `@family error metrics` to all 5, kept
  signatures/logic identical (no `library()` calls in source). `document()`
  wrote 5 `.Rd` + NAMESPACE (5 exports); `load_all()` clean. Bare
  `loess`/`predict`/`cor` resolve via the existing `@importFrom stats`.
- **§9 step 3 done:** added `R/msamodel-package.R` with the package-wide
  `@importFrom` directives (Option A namespacing — imports declared once;
  copied function files use bare calls). Inventory verified against the 16
  migrated source files; ggplot2 excluded (Suggests). `document()` populated
  NAMESPACE (incl. `importFrom(rlang,"!!")`, `importFrom(magrittr,"%>%")`, no
  exports yet); `load_all()` clean.
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
  `dev/RESUME.md`; added `dev/PROGRESS.md`.
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
