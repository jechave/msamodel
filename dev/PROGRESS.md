# Checklist — msamodel migration

Live checklist: a 1:1 mirror of `dev/plan.md` §9 (migration order). The
**normative spec is `dev/plan.md`**; the append-only history is `dev/LOG.md`.
Tick items as steps finish. When the plan changes, regenerate this list to match
it (and record the change in `dev/LOG.md`).

- [x] 1. Create package skeleton
- [x] 2. Copy LICENSE + DESCRIPTION
- [x] 3. Add `R/msamodel-package.R` (`@importFrom` directives)
- [x] 4. Copy `R/` files in dependency order:
  - [x] 4.1 `pdb_utils.R`, `enm_setup.R`, `site_properties.R`
  - [x] 4.2 `generate_spm_data.R`, `msa_bayesian_data_preparation.R`
  - [x] 4.3 `msa_model_evaluation.R`
  - [x] 4.4 `msa_bayesian_analysis.R`
  - [x] 4.5 `msa_decomposition.R`
  - [x] 4.6 `msa_a1a2grid_workflow.R`
  - [x] 4.7 `msa_bayesian_workflow.R`
- [x] 5. `devtools::document()` — confirm NAMESPACE populated
- [x] 6. Write + run `data-raw/prepare_znb_data.R` → `data/*.rda`
- [x] 7. Write `R/data-doc.R` (dataset docs); re-`document()`
- [x] 8. Write tests
- [x] 9. Write vignette (`msamodel-intro.Rmd`)
- [x] 10. `devtools::check()` clean
