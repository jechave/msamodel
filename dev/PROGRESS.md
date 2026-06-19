# Checklist — work item in progress

Checklist for the **work item currently being executed** (a slice of the current
release cycle — NOT a package version; the version is plain semver in
`DESCRIPTION`). Rewritten from that item's detailed plan when it starts; dormant
between items. The durable roadmap is `dev/plan.md`; the append-only history is
`dev/LOG.md`.

Tick items as substeps finish, and add a one-line dated entry to `dev/LOG.md`.

---

## DORMANT — no work item in progress

**`dr2*` naming convention — DONE 2026-06-18** (a work item of the 0.3.0 cycle).
All steps + verification done:
- Adopted msamodel convention `dr2_<indices>` for every `dr2`-family name msamodel
  *creates*: SPM list-cols `dr2_ijm`/`dr2_njm`; matrix fields likewise; exported
  `calculate_dr2_i_msa` / `calculate_dr2_n_msa`; objective col `dr2_i_msa`; param
  `spm_pp`.
- Deleted duplicate local `delta_structure_dr2`; `generate_spm_data` calls
  `penm::delta_structure_dr2i` directly (penm names NOT renamed/wrapped).
- `znb_spm` regenerated (Validation OK vs tmp_src, tol 1e-8; data unchanged).
- **Profile-invariance gate** (`test-profile-invariance.R`, captured pre-rename)
  reproduces bit-for-bit across BOTH the code rename AND the fixture regen — no
  numeric value moved.
- Convention recorded in CLAUDE.md, dev/plan.md, memory.
- `check()` at v0.1 baseline (0E/1W/2N). Full suite 62/0F.
- `NEWS.md` is the single `# (development version)` entry (mode arm + naming),
  grouped by change type; `DESCRIPTION` `0.3.0.9000` (semver, no letters).

**Next work item — defined step-by-step at execution time:** either the remaining
0.3.0 work (observed `dr2_i`/`dr2_n` profiles from homologous structures +
alignment; fit `dr2_n`; joint fit) OR the **0.4 motion arm** (`dh_ijm` → `dh_njm` +
`nh_njm`). See `dev/plan.md` + the `project_next_session` memory. **When the next
item starts:** enter plan mode, read the touched code, write the detailed plan, then
rewrite this checklist from it.
