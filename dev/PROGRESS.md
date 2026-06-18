# Checklist — current version in flight

Checklist for the version **currently being worked on**. Rewritten from that
version's detailed plan when the version starts; dormant between versions. The
durable roadmap is `dev/plan.md`; the append-only history is `dev/LOG.md`.

Tick items as substeps finish, and add a one-line dated entry to `dev/LOG.md`.

---

## DORMANT — no version in flight

**v0.3a SHIPPED 2026-06-18** (mode-form `dr2_n`, predict-only). Committed; the
post-ship R/ reorg + roadmap renumber are also committed and pushed
(`992e771`). See `dev/LOG.md`.

**Next version: v0.3b — consistency tidy** (NOT started; no detailed plan yet).
Scope (from the `dev/plan.md` v0.3b bullet): rename the SPM `dr2` column → `dr2i`;
replace local helper `delta_structure_dr2` with `penm::delta_structure_dr2i`;
deliberate `znb_spm` fixture regen (`data-raw/prepare_znb_data.R` + drift test).
**When v0.3b starts:** enter plan mode, read the touched code, write the detailed
plan, then rewrite this checklist from it.
