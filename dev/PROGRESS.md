# Progress — msamodel migration

Human-readable record of where the migration stands and what happened when.
The **normative spec is `dev/plan.md`**; the checklist below mirrors its §9
migration order. The "what's next" pointer lives in Claude's project memory
(loads automatically), not here — this file is history + at-a-glance status.

## §9 migration order

- [x] 1. Create package skeleton
- [x] 2. Copy LICENSE + DESCRIPTION
- [ ] 3. Add `R/msamodel-package.R` (`@importFrom` directives)  ← **next**
- [ ] 4. Copy `R/` files in dependency order:
  - [ ] 4.1 `utils.R`
  - [ ] 4.2 `loess_compare.R`
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
- `git init` (local-only, no remote) on `main`. **Initial commit done**
  (`33f2c46`): skeleton + `dev/plan.md`, `dev/PROGRESS.md`, `CLAUDE.md`.
  `tmp_src/` and `.claude/settings.local.json` excluded.

### 2026-04-22
- Created the package skeleton: `DESCRIPTION` (plan §3), `LICENSE` + `LICENSE.md`
  (MIT), empty roxygen-managed `NAMESPACE`, `.Rbuildignore`. `R/` left empty.
- Moved `plan.md` into `dev/`.
- Source-verification pass; answered the pre-flight questions.
