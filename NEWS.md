# msamodel (development version)

## New features

* **`calculate_lrmsd_i_nested_models()`** — new exported function returning the
  per-site lrmsd profiles of the four nested model variants (MM/MS/MA/MSA) at a
  given `(a1, a2)`, in one tibble (`i`, `pdb_site`, `lrmsd_i_mm`, `lrmsd_i_ms`,
  `lrmsd_i_ma`, `lrmsd_i_msa`). It is the single source of truth for the
  four-variant recipe: `calculate_prediction_samples()` (the MCMC path) now calls
  it per posterior sample instead of inlining the recipe, and a fixed-`(a1, a2)`
  decomposition (e.g. in the intro vignette) uses it directly.

* **Mode-form structural divergence.** Structural divergence is now also predicted
  **per normal mode**, not only per site (the first slice of the motion/mode arm).
  `generate_spm_data()` carries two new SPM list-columns, `mode` and `dr2_njm`
  (per-mode squared contribution to the mutant displacement), computed in the same
  per-mutant scan. New exported functions, parallel to the site path:
  * `preprocess_spm_mode()` — reshapes the scan into a `[mutant x mode]` matrix
    (`dr2_njm`); mode counterpart of `preprocess_spm()`.
  * `calculate_dr2_n_msa()` — selection-weighted mean `dr2_n` per mode; mode
    counterpart of `calculate_dr2_i_msa()`.
  * `calculate_loglik_lrmsd_n_msa()` / `fit_lrmsd_n_msa_ml()` — mode-form profiled
    Gaussian log-likelihood and its ML fit, parallel to the site path's
    `calculate_loglik_lrmsd_i_msa()` / `fit_lrmsd_i_msa_ml()`.

  No *empirical* per-mode divergence profile exists yet (deriving observed `dr2_n`
  from structural alignments is out of scope for `msamodel`), so the mode fit is
  exercised against a synthetic stand-in target; the site fit is unchanged. The new
  `mode-analysis` vignette walks through both the prediction and the (synthetic)
  fit.

## Breaking changes

* **`dr2`-family names now follow one index-signature convention:**
  `dr2_<indices>` — one underscore, then the free indices the object spans
  (response `i`/`n`, then mutated site `j`, then mutation `m`), letters joined; a
  reduction over an axis drops its letter. Concretely:
  * SPM list-columns are `dr2_ijm` (per site) and `dr2_njm` (per mode); each cell a
    per-mutant `(dr2_i)` / `(dr2_n)` vector. Affects `generate_spm_data()` output
    and the embedded `znb_spm`.
  * `preprocess_spm()` / `preprocess_spm_mode()` return the matrices in fields
    `dr2_ijm` / `dr2_njm`.
  * Exported `calculate_dr2_i_msa()` and `calculate_dr2_n_msa()` (the no-underscore
    `calculate_dr2i_msa()` / `calculate_dr2n_msa()` are gone). The profile columns
    they return (`dr2_i`, `dr2_n`) are unchanged.
  * Preprocess-bundle parameter renamed to `spm_pp`.

* **The local `delta_structure_dr2()` helper was removed.** It duplicated
  `penm::delta_structure_dr2i()` (verified bit-identical); `generate_spm_data()`
  now calls penm directly. penm's own (no-underscore) names are unchanged — the
  convention governs names msamodel *creates*, not what it calls.

* The embedded `znb_spm` dataset was regenerated for the new `mode` / `dr2_njm`
  columns and the `dr2_ijm` rename (data unchanged, verified against the migration
  source to 1e-8; profile values bit-for-bit unchanged).

# msamodel 0.2.0

API cleanup release — no new model capability. These are **breaking** changes
(public signatures and column names changed / removed). Appropriate for a 0.x,
GitHub-only package.

## Breaking changes

* **Site-level decomposition corrected and renamed `shap_*` → `phi_*`**
  (`phi_mut`, `phi_stab`, `phi_act`). The columns now hold the **sequential /
  additive** decomposition along the model progression M0 → MM → MS → MSA
  (`phi_mut = lrmsd_mm`, `phi_stab = lrmsd_ms − lrmsd_mm`,
  `phi_act = lrmsd_msa − lrmsd_ms`), the φ components used for the paper analysis.
  The earlier code computed a different, symmetric "Shapley" formula by mistake;
  both telescope to `lrmsd_msa`, so the values shipped under `shap_*`/`phi_*` in
  prior dev snapshots were wrong. Affects `calculate_msa_decomposition()`,
  `calculate_decomposition_samples()`, `calculate_decomposition_summary()`, and
  the `decomposition_summary` returned by `run_msa_bayesian_analysis()`.
  (`lrmsd_ma` is still a required input, reserved for a future `method` switch to
  the Shapley variant.)

* **`calculate_msa_decomposition()` is now a pure vector function.** It takes the
  four nested-model lrmsd *vectors* — `calculate_msa_decomposition(mm, ms, ma, msa)`
  — and returns a named list of the three phi vectors, instead of taking a tibble
  with hard-coded column names. The decomposition is context-free math, so the same
  function serves the site (`i`) and future mode (`n`) axes unchanged; the caller
  supplies whichever four columns it holds. The MCMC plumbing
  (`calculate_decomposition_samples()`) keeps its tibble interface and calls the
  vector function internally.

* **Nested-model lrmsd columns renamed `lrmsd_*` → `lrmsd_i_*`.** The four
  nested-model profiles are now `lrmsd_i_mm`, `lrmsd_i_ms`, `lrmsd_i_ma`,
  `lrmsd_i_msa` (parallel to `dr2_i`/`dr2_n`, leaving room for the mode arm's
  `lrmsd_n_*`). Affects the columns of `calculate_prediction_samples()`, the
  `variable` *values* in `prediction_summary` (e.g. `"lrmsd_i_msa"`), and the
  required input columns of `calculate_decomposition_samples()`.

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
  `purrr::map_dfr()` over the per-point divergence calculator reusing a single
  `preprocess_spm()` result (see the intro vignette, "Systematic 2-D scans").

## Internal

* `bio3d` moved from `Imports` to `Suggests` (used only by tests, the vignette,
  and `data-raw/`; still always present transitively via `penm`). `stringr`
  dropped from `Imports` (was used only by the removed `get_active_site()`).

# msamodel 0.1.0

* First release: structure-divergence profiles (structure × site), Bayesian
  (MCMC) parameter estimation, the four model variants (MM/MS/MA/MSA), and the
  site-level decomposition. Worked example shipped as the `znb_*` datasets
  (1znb_A, zinc beta-lactamase).
