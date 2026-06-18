# Checklist — current version in flight

Checklist for the version **currently being worked on**. Rewritten from that
version's detailed plan when the version starts; dormant between versions. The
durable roadmap is `dev/plan.md`; the append-only history is `dev/LOG.md`.

Tick items as substeps finish, and add a one-line dated entry to `dev/LOG.md`.

---

## v0.3b — `dr2*` naming convention (IN FLIGHT, started 2026-06-18)

Adopt msamodel convention `dr2_<indices>` for every `dr2`-family name msamodel
*creates*; call penm functions directly (no wrappers); delete the duplicate local
helper; normalize exported fn names. Full detail: approved plan + `dev/plan.md`
v0.3b bullet.

**Step 0 — freeze profiles BEFORE any rename (hard prerequisite):**
- [x] 0a. Add `tests/testthat/test-profile-invariance.R` (snapshot `dr2_i`, `dr2_n`
  profiles + loglik via `expect_snapshot_value(..., style="serialize")`)
- [x] 0b. Recorded `_snaps/profile-invariance.md` on CURRENT code; stable on re-run
  (PASS 3, no "Adding"). Note: file is `_snaps/profile-invariance.md` (single md),
  not a dir.

**Rename:**
- [ ] 1. `R/spm.R` — delete `delta_structure_dr2`; call `penm::delta_structure_dr2i`
  directly; vars/cols `dr2`→`dr2_ijm`, `dr2n`→`dr2_njm`; fields `dr2mat`→`dr2_ijm`,
  `dr2nmat`→`dr2_njm`; leave `delta_structure_dr` (raw displacement)
- [ ] 2. `R/model.R` — rename fns `calculate_dr2i_msa`→`calculate_dr2_i_msa`,
  `calculate_dr2n_msa`→`calculate_dr2_n_msa`; field access; param `spm_pp`
- [ ] 3. `R/objective.R` — `dr2_msa`→`dr2_i_msa`; field/var
- [ ] 4. `R/msa_bayesian_workflow.R` + `msa_bayesian_analysis.R` — field, param
  `spm_pp`, fn call-site renames
- [ ] 5. `R/msamodel-package.R` — `globalVariables` list
- [ ] 6. `R/data-doc.R` — `\item{dr2}`→`dr2_ijm`, `\item{dr2n}`→`dr2_njm`
- [ ] 7. `R/site_properties.R` — roxygen xref `calculate_dr2i_msa()`→`calculate_dr2_i_msa()`
- [ ] 8. Tests — symbol refs only (incl. renaming refs in test-profile-invariance.R;
  snapshot VALUES untouched)
- [ ] 9. `data-raw/prepare_znb_data.R` + deliberate `znb_spm` regen
- [ ] 10. Vignettes (`*.Rmd.orig` → re-knit → preview)
- [ ] 11. `man/` via `document()`
- [ ] 12. `DESCRIPTION` version + `NEWS.md` entry

**Verify:**
- [ ] V1. `test(filter='profile-invariance')` reproduces Step-0 snapshots (THE gate;
  re-run after step 9 regen too)
- [ ] V2. `document()` clean
- [ ] V3. `install()` + full `test()` pass
- [ ] V4. `znb_spm` regen = names-only, numbers identical
- [ ] V5. No dangling old names (grep)
- [ ] V6. Vignettes re-knit + preview OK
- [ ] V7. `check()` at v0.1 baseline (0E/1W/2N)

**Persist + reconcile:**
- [ ] P1. CLAUDE.md naming-convention subsection
- [ ] P2. Memory entry + MEMORY.md pointer
- [ ] P3. dev/plan.md v0.3b bullet (DONE 2026-06-18, ahead of code per project rule)
- [ ] DoD. PROGRESS dormant, LOG entry, project_next_session updated, git clean+pushed
