# Checklist — msamodel migration

Live checklist: a 1:1 mirror of `dev/plan.md` §9 (migration order). The
**normative spec is `dev/plan.md`**; the append-only history is `dev/LOG.md`.
Tick items as steps finish. When the plan changes, regenerate this list to match
it (and record the change in `dev/LOG.md`).

- [x] 1. Create package skeleton
- [x] 2. Copy LICENSE + DESCRIPTION
- [x] 3. Add `R/msamodel-package.R` (`@importFrom` directives)
- [ ] 4. Copy `R/` files in dependency order:
  - [ ] 4.1 `pdb_utils.R`, `enm_setup.R`, `site_properties.R`
  - [ ] 4.2 `spm_generate.R`, `spm_preprocess.R`
  - [ ] 4.3 `msa_evaluate.R`
  - [ ] 4.4 `msa_mcmc.R`
  - [ ] 4.5 `msa_decomposition_site.R`
  - [ ] 4.6 `a1a2grid.R`
  - [ ] 4.7 `msa_workflow.R`
- [ ] 5. `devtools::document()` — confirm NAMESPACE populated
- [ ] 6. Write + run `data-raw/prepare_znb_data.R` → `data/*.rda`
- [ ] 7. Write `R/data-doc.R` (dataset docs); re-`document()`
- [ ] 8. Write tests
- [ ] 9. Write vignette (`msamodel-intro.Rmd`)
- [ ] 10. `devtools::check()` clean
