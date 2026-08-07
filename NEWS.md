# msamodel (development version)

## New features

* **Structural divergence is predicted per normal mode, as well as per residue.**
  Every `calculate_*` and `predict_*` verb returns `list(site =, mode =)` — both
  response axes from one call, computed from the same mutation scan. The site axis
  answers "which residues diverge?"; the mode axis answers "which collective motions
  does divergence live in?".

  ```r
  spm  <- generate_spm_data(znb_wt, pdb_site_active = active, seed = 1024)
  prof <- calculate_profiles(spm, a1 = 1, a2 = 1)
  prof$site   # site, pdb_site, lrmsd_msa
  prof$mode   # mode, lrmsd_msa
  ```

  `fit_lrmsd_msa_mode()` fits the model to an observed per-mode profile. No empirical
  per-mode profile exists yet (deriving one from structural alignments is out of scope
  here), so the shipped `znb_profile_n` is a synthetic stand-in and the mode fit is
  exercised against it. The `mode-analysis` vignette walks through both.

* **Prediction bands include the scan's own sampling error.** The `_se` columns from
  `predict_profiles()` / `predict_decomposition()` combine two independent sources: the
  `(a1, a2)` parameter uncertainty, and the fact that each divergence value is a mean
  over a *finite* mutant scan rather than the infinite limit. Build a band as
  `value ± qnorm(0.975) * se`.

  A consequence worth knowing: the mutation-only (MM) variant has no parameter
  uncertainty but still carries a real, nonzero band from the finite scan.

* **The divergence profile can be split into its causes.**
  `calculate_decomposition()` / `predict_decomposition()` return the four nested model
  variants (MM/MS/MA/MSA) together with the mutation, stability, and activity
  contributions (`phi_mut`, `phi_stab`, `phi_act`), which sum exactly to the full
  profile. Currently available for `which = "nlrmsd"`; `which = "lrmsd"` errors, as the
  uncentred standard error has not been derived.

## API changes

* **One scan object feeds everything.** `generate_spm_data()` returns a ready-to-use
  `spm` object carrying both response axes; pass it straight to any `calculate_*`,
  `fit_*`, or `predict_*` function. There is no separate preprocessing step.

* **Observations are a vector pair, not a data frame.** The fitters take
  `(spm, pdb_site, lrmsd_obs)` / `(spm, mode, lrmsd_obs)`, so your own table can name
  its columns anything:

  ```r
  fit_lrmsd_msa_site(spm, my_data$residue, my_data$divergence)
  ```

  `znb_profile` / `znb_profile_n` are example data, not a required schema. Inputs are
  validated at the boundary: equal length, non-empty, no `NA`.

* **Goodness of fit comes with the fit.** `fit$gof` is a one-row tibble (`D2`, `AIC`,
  `BIC`, `logLik`, `deviance`, `null_deviance`, `nobs`, `k`, `sigma_hat`); the fit's top
  level carries the estimate (`a1`, `a2`, `logLik`, `cov`, `se_a1`, `se_a2`,
  `convergence`) plus `call`. There is no separate accessor to call.

* **Value columns carry no axis tag.** Both branches emit the same vocabulary
  (`lrmsd_msa`, `nlrmsd_mm`, `nphi_stab`, ...); only the key column differs — `site` +
  `pdb_site` on one axis, `mode` on the other. Names like `lrmsd_i_msa` from earlier
  development versions no longer exist.

* **`predict_decomposition()` groups its `_se` columns.** The seven value columns come
  first, then the seven matching `_se` columns, rather than alternating.

## Breaking changes

* **Inference is maximum-likelihood with delta-method bands.** The MCMC fitter, its
  sample-averaging analysis chain, and the adaptive Gauss–Hermite quadrature arm that
  briefly replaced it were all removed. Fitting is `fit_lrmsd_msa_site()` /
  `fit_lrmsd_msa_mode()`; bands come from `predict_profiles()` /
  `predict_decomposition()` — deterministic, with no seed, burn-in, or draw-averaging.

* **The exported surface is ten functions:** `setup_enm()`, `generate_spm_data()`,
  `add_site_properties()`, `pfix_msa()`, `calculate_profiles()`,
  `calculate_decomposition()`, `fit_lrmsd_msa_site()`, `fit_lrmsd_msa_mode()`,
  `predict_profiles()`, `predict_decomposition()`. The earlier per-axis grid of
  `calculate_*` / `predict_*_ml` functions, the `gof_*` accessors, `weights_jm_spm()`,
  and the `preprocess_spm*()` reshapers are gone or internal.

* **The embedded `znb_spm` dataset was regenerated** for the mode axis. Values are
  unchanged (verified to 1e-8 against the migration source).

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
  prior dev snapshots were wrong. Affects the shared decomposition kernel and every
  `phi_*` output built on it. (`lrmsd_ma` is still a required input, reserved for a
  future `method` switch to the Shapley variant.)

* **`calculate_msa_decomposition()` is now a pure vector function.** It takes the
  four nested-model lrmsd *vectors* — `calculate_msa_decomposition(mm, ms, ma, msa)`
  — and returns a named list of the three phi vectors, instead of taking a tibble
  with hard-coded column names. The decomposition is context-free math, so the same
  function serves the site (`i`) and future mode (`n`) axes unchanged; the caller
  supplies whichever four columns it holds.

* **Nested-model lrmsd columns renamed `lrmsd_*` → `lrmsd_i_*`.** The four
  nested-model profiles are now `lrmsd_i_mm`, `lrmsd_i_ms`, `lrmsd_i_ma`,
  `lrmsd_i_msa` (parallel to `dr2_i`/`dr2_n`, leaving room for the mode arm's
  `lrmsd_n_*`).

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
  `preprocess_spm()` result (see the site analysis vignette, "How the profile
  depends on a1 and a2").

## Internal

* `bio3d` moved from `Imports` to `Suggests` (used only by tests, the vignette,
  and `data-raw/`; still always present transitively via `penm`). `stringr`
  dropped from `Imports` (was used only by the removed `get_active_site()`).

# msamodel 0.1.0

* First release: structure-divergence profiles (structure × site), Bayesian
  (MCMC) parameter estimation, the four model variants (MM/MS/MA/MSA), and the
  site-level decomposition. Worked example shipped as the `znb_*` datasets
  (1znb_A, zinc beta-lactamase).
