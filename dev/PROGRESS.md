# Checklist — current version in flight

Checklist for the version **currently being worked on**. Rewritten from that
version's detailed plan when the version starts; dormant between versions. The
durable roadmap is `dev/plan.md`; the append-only history is `dev/LOG.md`.

Tick items as substeps finish, and add a one-line dated entry to `dev/LOG.md`.

---

## Current version: v0.2 — API improvements

Started 2026-06-11. No new model capability; four breaking API changes, exports
18 → 14. `tmp_src/` untouched. Detailed plan: the four decisions are recorded in
the `dev/plan.md` v0.2 bullet.

- [x] Step 0 — plan docs (plan.md v0.2 bullet + precompute wording, this
      checklist, LOG entry)
- [x] Step 1 — rename `shap_*` → `phi_*` (R/msa_decomposition.R columns +
      validation + `starts_with`; test-decomposition; vignette labels/prose;
      roxygen prose in workflow; DESCRIPTION; CLAUDE.md note)
- [x] Step 2 — structure input = bio3d pdb object (dropped `load_protein`; added
      `inherits(pdb,"pdb")` check to `setup_enm`; callers use `bio3d::read.pdb`)
- [x] Step 3 — active-site input = integer vector (dropped `get_active_site`;
      inlined `pdb_site_active` in data-raw/tests/vignette; kept `znb_dataset`
      as illustrative; moved bio3d→Suggests, dropped stringr import)
- [x] Step 4 — drop grid API (deleted `R/msa_a1a2grid_workflow.R` +
      `test-a1a2grid.R`; removed grid assertion in test-contract; replaced
      vignette grid section with a map_dfr scan; reworded `site_properties.R`)
- [x] Step 5 — `document()` (14 exports) + `test()` (45 pass) + mmCIF verified
      (read.cif→setup_enm OK) + reran `precompute.R` (cache now `phi_*`) +
      `check()` = 0E/1W/2N (= v0.1 baseline) + DESCRIPTION 0.2.0 + NEWS.md

**v0.2 COMPLETE — committed + pushed (`9f8338a`) 2026-06-11.** No version in
flight. Next on the roadmap: **v0.3 — motion/mode arm (region D, `.archive/`)**.
Per `dev/plan.md`, v0.3 starts by reading the rest of `.archive/` (incl.
`model_rates.R`) line-by-line in plan mode, then writing its detailed plan and
rewriting this checklist from it. This file is dormant until then.

**Post-v0.2 follow-up (2026-06-11):** reworked the intro vignette to the rOpenSci
`.Rmd.orig` precompute pattern (one source of truth, shown code == executed code),
removing the `precompute.R` + `vignette_cache.rds` mechanism that had drifted and
silently skipped figures. `check()` still at baseline. Details in `dev/LOG.md`.
Not a v0.2 substep — a correctness fix to already-shipped work.
