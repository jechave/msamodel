# Checklist — current version in flight

Checklist for the version **currently being worked on**. Rewritten from that
version's detailed plan when the version starts; dormant between versions. The
durable roadmap is `dev/plan.md`; the append-only history is `dev/LOG.md`.

Tick items as substeps finish, and add a one-line dated entry to `dev/LOG.md`.

---

## DORMANT — no version in flight

**v0.3b SHIPPED 2026-06-18** (`dr2*` naming convention). All steps + verification
done:
- Adopted msamodel convention `dr2_<indices>` for every `dr2`-family name msamodel
  *creates*: SPM list-cols `dr2`→`dr2_ijm`/`dr2n`→`dr2_njm`; matrix fields likewise;
  exported `calculate_dr2i_msa`→`calculate_dr2_i_msa`, `dr2n`→`dr2_n`; objective col
  `dr2_msa`→`dr2_i_msa`; param `spm_energies_and_dr2mat`→`spm_pp`.
- Deleted duplicate local `delta_structure_dr2`; `generate_spm_data` calls
  `penm::delta_structure_dr2i` directly (penm names NOT renamed/wrapped).
- `znb_spm` regenerated (Validation OK vs tmp_src, tol 1e-8; data unchanged).
- **Profile-invariance gate** (`test-profile-invariance.R`, captured pre-rename)
  reproduces bit-for-bit across BOTH the code rename AND the fixture regen — no
  numeric value moved.
- Convention recorded in CLAUDE.md, dev/plan.md, memory.
- `check()` at v0.1 baseline (0E/1W/2N). Full suite 62/0F.
- DESCRIPTION bumped 0.2.0.9000 → 0.3.0.9000; NEWS.md v0.3b section + stale-ref fixes.

**Next version: v0.3c+ — TBD step-by-step** (observed `dr2_i`/`dr2_n` profiles from
homologous structures + alignment; fit `dr2_n`; joint fit) OR **v0.4 motion arm**
(`dh_ijm` → `dh_njm` + `nh_njm`). See `dev/plan.md` + `project_next_session` memory.
**When the next version starts:** enter plan mode, read the touched code, write the
detailed plan, then rewrite this checklist from it.
