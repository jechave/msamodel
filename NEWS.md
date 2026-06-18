# msamodel (development version)

## Breaking changes (v0.3b — `dr2*` naming convention)

* **`dr2`-family names now follow one index-signature convention:**
  `dr2_<indices>` — one underscore, then the free indices the object spans
  (response `i`/`n`, then mutated site `j`, then mutation `m`), letters joined; a
  reduction over an axis drops its letter. Concretely:
  * SPM list-columns `dr2` → **`dr2_ijm`**, `dr2n` → **`dr2_njm`** (each cell a
    per-mutant `(dr2_i)` / `(dr2_n)` vector). Affects `generate_spm_data()` output
    and the embedded `znb_spm`.
  * `preprocess_spm()` / `preprocess_spm_mode()` matrix fields `dr2mat` /
    `dr2nmat` → **`dr2_ijm`** / **`dr2_njm`**.
  * Exported `calculate_dr2i_msa()` → **`calculate_dr2_i_msa()`**,
    `calculate_dr2n_msa()` → **`calculate_dr2_n_msa()`**. The profile columns they
    return (`dr2_i`, `dr2_n`) were already convention-correct and are unchanged.
  * Preprocess-bundle parameter `spm_energies_and_dr2mat` → **`spm_pp`**.

* **The local `delta_structure_dr2()` helper was removed.** It duplicated
  `penm::delta_structure_dr2i()` (verified bit-identical); `generate_spm_data()`
  now calls penm directly. penm's own (no-underscore) names are unchanged — the
  convention governs names msamodel *creates*, not what it calls.

* `znb_spm` was regenerated with the new column names (data unchanged, verified
  against the migration source to 1e-8; profile values bit-for-bit unchanged).

## New features (v0.3a — mode-form structural divergence)

* **Mode-form structural divergence (`dr2_n`).** The first slice of the
  motion/mode arm: structural divergence is now also predicted **per normal
  mode**, not only per site. `generate_spm_data()` carries two new SPM
  list-columns, `mode` and `dr2_njm` (per-mode squared contribution to the mutant
  displacement; renamed from `dr2n` in v0.3b), computed in the same per-mutant
  scan. Two new exported functions, parallel to the site path:
  * `preprocess_spm_mode()` — reshapes the scan into a `[mutant x mode]` matrix
    (`dr2_njm`); mode counterpart of `preprocess_spm()`.
  * `calculate_dr2_n_msa()` — selection-weighted mean `dr2_n` per mode; mode
    counterpart of `calculate_dr2_i_msa()`. (Both renamed from the no-underscore
    form in v0.3b.)

  The mode arm is **predict-only** (there is no observed mode profile to fit), so
  it has no log-likelihood counterpart and the site fit is unchanged. New
  `dr2n-analysis` vignette walks through it.

* The embedded `znb_spm` dataset was regenerated to include the `mode` /
  `dr2_njm` columns. (Its other columns are unchanged.)

# msamodel 0.2.0

API cleanup release — no new model capability. These are **breaking** changes
(public signatures and column names changed / removed). Appropriate for a 0.x,
GitHub-only package.

## Breaking changes

* **Decomposition columns renamed `shap_*` → `phi_*`** (`phi_mut`, `phi_stab`,
  `phi_act`). "Shapley" was a misnomer; these are the φ components from the
  paper. Affects `calculate_msa_decomposition()`,
  `calculate_decomposition_samples()`, `calculate_decomposition_summary()`, and
  the `decomposition_summary` returned by `run_msa_bayesian_analysis()`.

* **Structure input is now a bio3d pdb object.** `load_protein()` (which took a
  `pdb_chain` ID + a `data_dir`) is **removed**. Read the structure yourself and
  pass the object to `setup_enm()`:

      pdb <- bio3d::read.pdb(path)   # legacy .pdb
      pdb <- bio3d::read.cif(path)   # mmCIF (bio3d's read.cif is beta)
      wt  <- setup_enm(pdb, ...)

  The package now does no file I/O, so any source that yields a bio3d pdb object
  works (RCSB, AlphaFold-DB, mmCIF). `setup_enm()` validates that its input is a
  bio3d pdb object. mmCIF was confirmed to round-trip through `setup_enm()`;
  note bio3d flags `read.cif()` itself as beta.

* **Active-site input is now a plain integer vector** of PDB residue numbers,
  passed as `pdb_site_active` to `generate_spm_data()` / `add_site_properties()`
  (already supported). `get_active_site()` (which looked up a bundled
  `dataset_ec2024.csv`) is **removed**. The `znb_dataset` dataset is retained as
  an *illustrative* example of the source active-site table, but no package
  function consumes it.

* **The a1/a2 grid API is removed.** `define_selection_grid()` and
  `calculate_dr2i_msa_a1a2grid()` are gone. A systematic 2-D scan is a short
  `purrr::map_dfr()` over `calculate_dr2_i_msa(pp, a1, a2)` (renamed from
  `calculate_dr2i_msa()` in v0.3b) reusing a single `preprocess_spm()` result (see
  the intro vignette, "Systematic 2-D scans").

## Internal

* `bio3d` moved from `Imports` to `Suggests` (used only by tests, the vignette,
  and `data-raw/`; still always present transitively via `penm`). `stringr`
  dropped from `Imports` (was used only by the removed `get_active_site()`).

# msamodel 0.1.0

* First release: structure-divergence profiles (structure × site), Bayesian
  (MCMC) parameter estimation, the four model variants (MM/MS/MA/MSA), and the
  site-level decomposition. Worked example shipped as the `znb_*` datasets
  (1znb_A, zinc beta-lactamase).
