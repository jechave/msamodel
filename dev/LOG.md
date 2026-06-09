# Log — msamodel migration

Append-only history of what was done and attempted (newest first), including
dead ends and decisions. The **normative spec is `dev/plan.md`**; the live
checklist is `dev/PROGRESS.md` (a 1:1 mirror of plan.md §9). This file is the
history those two don't keep.

One short entry per working session.

### 2026-06-09 (later)
- **§9 step 8 done — test suite written (8 files, 55 tests, all pass).** Scaffolded
  `tests/testthat.R` + `tests/testthat/`. Files: spm-generate (the fixture-drift
  guard — reproduces znb_wt/znb_spm from the §5 params, ~20s, dominates runtime),
  spm-preprocess (dr2mat shape + the site_map mapping the pdb_site contract relies
  on), msa-evaluate (dr2i finiteness + a FROZEN loglik literal -184.3241923285),
  msa-mcmc (prior-range validation errors + deterministic fix_a1/fix_a2 constancy
  + shape), decomposition (validation stop()s + the pdb_site optional branch both
  ways + summary shape — additivity NOT asserted, it's tautological), a1a2grid
  (54-row grid + output cols + pdb_site alignment), workflow-endtoend (7-element
  return, real pipeline), contract (unknown-pdb_site error, subset-coverage finite
  loglik, independent-route mapping equality, pdb_site on outputs). Suite: 55 pass
  / 0 fail / 0 warn, ~29s. **Failure-sanity check:** perturbed the frozen loglik
  and the dr2mat shape → the relevant tests FAILED, then reverted → clean. Confirms
  assertions aren't vacuous.
- **Plan §7 reconciled for step 8 (tests).** Updated the test matrix to the
  pdb_site contract and to vetted assertions (the prior matrix had tautological /
  untestable rows). Key changes: dropped the Shapley-additivity assertion (it is
  algebraically the implementation — can't fail except by editing the formula);
  replaced the loglik "regression" with a FROZEN literal (captured once, not
  recomputed in-test); dropped the "acceptance rate > 0" claim (`accepts` isn't
  returned); added deterministic `fix_a1/fix_a2` and validation-error tests for
  MCMC; made decomposition test validation + the pdb_site optional-branch + shape.
  Added **`test-contract.R`** (8th file) for the pdb_site contract: unknown-pdb_site
  error, subset-coverage finite loglik, independent-route mapping equality (non-
  circular), pdb_site on outputs. v0.1 test count 7 → 8. Done before writing tests
  (fix-plan-first). See memory `feedback-test-quality`.

### 2026-06-09
- **Fixed the `observed_data` site-key contract (pdb_site, not i) + documented
  datasets (step 7).** Surfaced while prepping step 7: `observed_data` is
  *user-provided* input (public arg of `run_msa_bayesian_analysis`), but the
  contract required the user to supply `i` — the model's internal 1..N response-
  site renumbering — which a user can't reliably produce, and the fit joined on it
  with no validation (silent-misalignment risk). For 1znb_A it worked by luck
  (verified 225/225 pdb_site agreement), not by construction.
  **Decision (user): full fix overriding CLAUDE.md §4 "rename nothing"** — justified
  as a public-API correctness fix. New contract: `observed_data = {pdb_site,
  lrmsd_obs}`; package maps pdb_site→i internally via a new `site_map` returned by
  `preprocess_spm`, erroring on any unknown pdb_site. All six site-keyed objects
  now also carry `pdb_site` alongside `i` (prediction_samples/summary,
  decomposition_samples/summary, observed_data passthrough, grid output). Internal
  i/j convention (i=response site, j=perturbed site) unchanged. Files:
  `R/msa_bayesian_data_preparation.R` (site_map), `R/msa_model_evaluation.R`
  (map+validate+join), `R/msa_bayesian_analysis.R`, `R/msa_decomposition.R`
  (pdb_site optional — synthetic test input without it still works),
  `R/msa_a1a2grid_workflow.R`, `R/msa_bayesian_workflow.R` (docs).
  **`znb_profile` regenerated** as `{pdb_site, lrmsd_obs}` (dropped i; dropped the
  derivable dactive/lrmsf — verified reproducible from znb_wt to machine
  precision: dactive=get_dactive, lrmsf=log(sqrt(get_msf_site))). SPM one-time
  validation still passed on re-run.
  **Step 7:** wrote `R/data-doc.R` documenting all 5 datasets (`@family datasets`,
  accurate `@format`; znb_profile framed as the {pdb_site, lrmsd_obs} fit target).
  **Verify:** correctness gate — new contract gives loglik identical to the old
  i-join reference for 1znb_A (-184.3242, all.equal TRUE); negative test — unknown
  pdb_site errors clearly; `document()` clean, NAMESPACE still 18 exports;
  `load_all()` clean; full workflow returns 7 elems with pdb_site on every
  site-keyed output. Updated `dev/plan.md` §2/§5 (contract note + deliberate §4
  exception) BEFORE coding. v0.2 maybe: make pdb_site the canonical key throughout.

### 2026-06-06
- **§9 step 6 done — embedded `znb_*` data, SPM generated + validated.** Vendored
  the small inputs (`inst/extdata/1znb_A.pdb`; `data-raw/raw/dataset_1znb_A.csv`
  1 row; `data-raw/raw/profiles_1znb_A.csv` 225 rows), wrote
  `data-raw/prepare_znb_data.R` (generates `znb_spm` via the package's own
  `generate_spm_data`, params as named constants flagged to match the step-8
  test), ran it. **One-time validation PASSED:** `all.equal(znb_spm,
  readRDS(tmp_src/.../1znb_A_spm.rds), tol=1e-8)` is TRUE — the migrated
  `generate_spm_data` reproduces the source SPM (migration-correctness proof;
  block fenced as removable). Saved 5 `.rda` via `use_data(compress="xz")`.
  Sanity: znb_spm 2508 rows, znb_profile 225, 8 active sites, dr2mat 2280×228.
  Added `^data-raw$` to `.Rbuildignore`; wrote `data-raw/source_files.txt`. Scan
  runtime ~20 s. **Flag for step 10:** `znb_spm.rda` ≈ 13 MB and `znb_wt.rda`
  ≈ 6 MB even xz-compressed — will trigger an `R CMD check` data-size NOTE; needs
  a decision later (subsample? fewer mutations? accept the NOTE?). Not a step-6
  blocker.
- **Plan change — §5 embedded-data approach: generate the SPM, don't copy it.**
  Rejected the original §5 (copy `tmp_src/.../1znb_A_spm.rds`, point `SRC` at
  `tmp_src/`/the live project): not self-consistent — once `tmp_src/` is deleted,
  `data-raw/` would have no in-repo recipe to rebuild `data/`, and pointing at the
  live project isn't reproducibility. **New approach:** vendor only the *small*
  inputs (PDB → `inst/extdata/`, the 1znb_A rows of the two CSVs →
  `data-raw/raw/`), and **generate** `znb_spm` with the package's own
  `generate_spm_data()`. Validate the generated SPM **once** against the
  `tmp_src/` original via `stopifnot(all.equal(...))` (migration-correctness
  proof, fenced as removable) — this is NOT the permanent test; the
  regenerate-and-compare drift guard remains step 8's `test-spm-generate.R`.
  **Feasibility verified before deciding:** recovered the generation params from
  `tmp_src/scripts/03_scan_mutants_all_cases.R` (n_mutations=10, model="lfenm",
  sigma=0.3, min_sd=2, seed=1024; ENM node="ca", ming_wall, d_max=10.5,
  frustrated=FALSE); confirmed the scan is seed-deterministic — penm's
  `get_mutant_site_lfenm` does `set.seed(seed + site_mut*mutation)` per mutant
  (tested: same seed → identical mutant xyz, different seed → different;
  penm 0.2.0.9000, R 4.2). The .rds carries no provenance, so the params are
  hard-coded in the script and must match the step-8 test (option a). **Edited
  `dev/plan.md` §5 (rewrite), §2 (layout: znb_spm "generated"; data-raw adds
  `raw/`), §9 step 6, §13 decision 5.** Done as a standalone turn *before* writing
  any step-6 code, so a fresh context reads the corrected plan (per
  `CLAUDE.md` "read plan.md before doing anything"). Step-6 code execution is the
  next turn.

### 2026-06-05
- **§9 step 4 done (all sub-steps 4.1–4.7), in one batch + step 5.** Copied the 7
  v0.1 source files from `tmp_src/R/` into `R/` keeping source filenames (per the
  no-rename change below), in dependency order. Applied §4 rules: stripped all
  `@requires` tags; added `@export` to the 18 public fns + `@family` groups; the
  2 helpers (`delta_structure_dr`, `delta_structure_dr2`) marked `@noRd`. No
  inline `library()` existed to strip. **Bug fixes:** §11.1 — deleted the
  duplicate `define_selection_grid` (the `a1_range, a2_range` override at source
  lines 106–119), keeping the real defaulted one; §11.3 — fixed the
  `"lower,Quite upper"` roxygen typo in `calculate_prediction_summary`'s
  `@return`. **4.7 assessment removal:** dropped the two
  `calculate_model_comparison_samples/summary` blocks + their return-list entries
  from `run_msa_bayesian_analysis` → 7-element return; updated its `@return`.
  **Imports:** added `transmute` + `group_split` to the dplyr `@importFrom` in
  `msamodel-package.R` (used by `msa_decomposition.R`). **Verify:**
  `devtools::document()` clean (no roxygen warnings; 18 `.Rd` written, none for
  the `@noRd` helpers); NAMESPACE has exactly 18 exports, `delta_structure_*`
  absent; `devtools::load_all()` clean. Behavioral/numeric checks deferred to
  steps 6–8 (no data/tests yet).
- **Plan change — no file renames in v0.1.** Dropped the §2 rename list (10
  renames) before starting §9 step 4. Files are now copied from `tmp_src/R/`
  keeping their **source filenames**. Rationale: source names keep provenance
  obvious during the migration, and some rename rationales were v0.1-scoped (e.g.
  `msa_bayesian_workflow.R` → `msa_workflow.R` "because only one workflow in
  v0.1") — an assumption v0.2 breaks, which would force a *second* rename (churn
  the migration explicitly avoids). File naming/layout is now an explicit later
  refactor, not folded into the copy. Edited `dev/plan.md` §2 (rename block + `R/`
  layout tree) and §9 step-4 sub-steps; regenerated the §9-step-4 labels in
  `dev/PROGRESS.md` to match. Deferred-file renames mentioned elsewhere in the
  plan are left as advisory (decided when those files are actually migrated).

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
