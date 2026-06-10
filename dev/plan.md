# `msamodel` — roadmap and durable findings

This is the **coarse, durable** planning document: what each version *is* (goal +
rough scope), plus findings that must not be re-discovered. It is deliberately
high-level and stable — it is edited only when a version's *goal* changes, not when
an execution detail surfaces.

**Detailed planning is done per-version, at execution time.** When a version
starts, enter plan mode, read the specific code that version touches, and write the
executable detail *then* (as a fresh per-version plan, tracked in `dev/PROGRESS.md`
for the in-flight version). This split is deliberate: the v0.1 plan tried to be both
the blueprint and the step-by-step script, grew to 575 lines, and had to be amended
every time execution hit a detail it hadn't anticipated. Deep code exploration
belongs where it's actionable.

- `dev/PROGRESS.md` — checklist for the **version currently in flight** (rewritten
  from that version's detailed plan when it starts; dormant between versions).
- `dev/LOG.md` — append-only history (sessions, decisions, dead ends). Holds the
  v0.1 detail this document no longer carries.

The "record decisions before coding" rule still holds — it now means *write the
per-version detailed plan before coding* (and touch this roadmap only if the
version's goal changed).

---

## Source for continuing migration (READ-ONLY)

The package was migrated out of a paper project. The frozen source snapshot lives
at `tmp_src/` (build-ignored, a one-time copy of the live project, which is not
reachable from the session). **Treat `tmp_src/` as read-only** — copy *out of* it,
never edit *in* it. Much capability is still unmigrated (see the roadmap + survey
findings below). The migration *rules* (strip `library()`, route deps through
`Imports:` + `@importFrom`, keep signatures, strip `@requires` tags, add `@export`
per file) still apply when migration continues — see CLAUDE.md.

---

## Roadmap (versioned)

Ordering reflects difficulty and intent: API cleanup first, then the motion/mode
arm (reuses existing machinery), then the harder tree route. Scope is rough on
purpose; per-version detail is written when each version starts.

- **v0.1 — SHIPPED.** Compute structure-divergence profiles (structure × site) +
  Bayesian (MCMC) fit + site-level decomposition + a1/a2 grid. 18 exports. History
  in `dev/LOG.md`.

- **v0.2 — API improvements (no new model capability).**
  - Rename `shap_*` → `phi_*` (columns, roxygen, tests, vignette). "Shapley" is a
    misnomer; the user calls these φ in the paper.
  - PDB input: accept a `.pdb` file *path*, not a `pdb_chain` ID + bundled-CSV
    layout (so AFDB and other sources work).
  - Active-site input: accept a plain `pdb_site` integer vector; drop the hidden
    dependency on the bundled `dataset_ec2024.csv` lookup inside `get_active_site`.
  - a1/a2 grid: decide keep / demote-to-internal / drop. `define_selection_grid` is
    a trivial `expand.grid` wrapper; `calculate_dr2i_msa_a1a2grid` carries real
    logic (preprocess + 4 variants + site-map join).

- **v0.3 — motion/mode via the SPM-mean route (region D).** Extend the SPM to carry
  the extra per-mutant divergence columns (`dr2_njm`, `dh_ijm`, `dh_njm`, `nh_njm`)
  and add the reweighting + evaluation paths. Reuses the v0.1 precompute-and-reweight
  machinery (see findings). Read the rest of `.archive/` (incl. `model_rates.R`)
  line-by-line at the *start* of this version, and decide then whether rates land
  here or later.

- **v0.4 — tree / trajectory route (region C).** Migrate `someday_maybe/tree/`; fix
  the `p_act`→`p_ma` bug; swap `akima`→`interp`; reconcile with the SPM-mean route.
  Genuinely harder — no clean precomputation (see findings).

- **Later / undecided (ordering deferred):**
  - Assessment layer (region A 3-file bundle: `utils.R` + `compare_loess_fits.R` +
    `model_comparison_functions.R` — mutually dependent; GOF / model comparison),
    plus site/protein **allotment** and protein-level decomposition.
  - Visualization module — the modular figure library
    (`tmp_src/Rmd/msa_profiles_analysis/R/fig_*.R`, 11 scripts) and the deferred
    `plot_loess_comparison`.
  - AlphaFold2 profile validation (`tmp_src/Rmd/preliminary/test_af_profiles.Rmd`).
  - Sequence-vs-structure divergence
    (`tmp_src/Rmd/preliminary/sequence_divergence/`).

---

## Durable findings (do not re-discover)

### The model's full divergence grid (most of it not yet migrated)

The original model (in `tmp_src/.archive/`) computed a **2×2 grid** of divergence
quantities × {site `i`, mode `n`}, plus `nh_n`. v0.1 implements exactly ONE cell.

| quantity | site form | mode form | penm fn | in v0.1? |
|---|---|---|---|---|
| structure divergence `dr2` | `dr2_ijm` | `dr2_njm` | `delta_structure_dr2i` / `dr2n` | site-only |
| motion/fluctuation divergence `dh` | `dh_ijm` | `dh_njm` | `delta_motion_dhi` / `dhn` | **NO** |
| mode fluctuation `nh` | — | `nh_njm` | `delta_motion_nhn` | **NO** |

Also in `.archive`, characterized but not yet read line-by-line: `model_rates.R`
(evolutionary-rate prediction) and the `motion/` / `structure/` / mode-analysis
backup files. (Read at v0.3 start.)

### Precomputation property — RESOLVED: it generalizes for the SPM-mean route

Every divergence quantity in the archive is computed as

```
X = sum(pfix_jm * X_jm) / sum(pfix_jm)        # pfix_jm = pstab(a1) * pact(a2)
```

a1/a2 enter ONLY through the per-mutant `pfix_jm`. The `X_jm` columns are mutant
properties computed *before* selection — independent of a1/a2. So the v0.1
mechanism — a `[mutant × site]` matrix reweighted by `weights_jm`
(`R/msa_model_evaluation.R:12-26`, built in `R/msa_bayesian_data_preparation.R`) —
**generalizes**: each new quantity is just another value-matrix (`[mutant × site]`
or `[mutant × mode]`) reweighted by the *same* weights. The SPM must carry the
extra `_jm` columns. **Motion/mode under the SPM-mean route ≈ "more SPM columns +
more reweighting" — mechanically close to v0.1.** (Source: archive
`calculate_diff_site_msa.R:97-108`, `calculate_diff_mode_msa.R:76-91`.)

### The tree route does NOT share this property

Tree lineages are *sequential* substitutions; each fixation depends on accumulated
state, so the weights do not factorize per-mutant. That is exactly why the tree
code interpolates the likelihood over an a1/a2 grid (`akima::interp`) instead of
reweighting a fixed matrix. The clean precomputation is unavailable for the tree —
this is why v0.4 (tree) is harder than v0.3 (motion/mode).

### Other survey facts

- The assessment layer is a 3-file *bundle*: `model_comparison_functions.R` depends
  on both `compare_loess_fits.R` and `utils.R`. Migrate together or not at all.
- `tmp_src/R/original_compare_loess_fits.R` is a confirmed duplicate of
  `compare_loess_fits.R` — drop it, don't migrate.
- `tmp_src/R/protein_fit12_*decomposition.R` are FIT12-benchmark, paper-specific.

---

## Resolved decisions (durable)

- Target path: `/Users/julianechave/Library/Mobile Documents/com~apple~CloudDocs/lab/Rpackages/msamodel/`.
- Git: local-only, no remote; direct pushes to `main` are intended (solo repo).
- No file renames folded into migration copies (renaming is a later deliberate
  refactor). No public function name changes except deliberate, recorded API
  changes (e.g. the v0.2 `phi_*` rename, the v0.1 `observed_data` = {pdb_site,
  lrmsd_obs} contract).
- Embedded `znb_*` fixture is frozen and *generated* (not copied) by
  `data-raw/prepare_znb_data.R`; drift caught by `test-spm-generate.R`.
