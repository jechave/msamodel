# msamodel 0.4.0

## Breaking changes

* **`generate_spm()` takes `ensemble`, not `seed`, and it is required.** The argument
  was never a seed: penm hashes it together with `(site, mutation)` and *that* hash
  seeds the RNG. It names **which realization of the mutational process** the scan
  draws from. Changing it does not add randomness -- it moves to a different
  realization, in which a given mutation is a different mutation.

  ```r
  spm <- generate_spm(wt, ..., seed = 1024)      # was
  spm <- generate_spm(wt, ..., ensemble = 1L)    # now
  ```

  There is no default: a scan whose realization is not recorded cannot be regenerated.
  (The former `seed = NULL` default did not work either -- it errored on the first
  mutant.) See `?penm::penm_ensemble`.

* **All scan-derived numbers changed.** This is a correctness fix in penm, not drift in
  msamodel: the old RNG key `seed + site * mutation` collided -- every divisor pair of
  the same product shared a random stream, so 1710 of 2280 mutants in a standard scan
  were not independent draws. penm now hashes the key. Re-running an analysis therefore
  gives different numbers; they are a new, and genuinely independent, realization of the
  same process.

  The change was accepted on a similarity check rather than by re-freezing whatever came
  out: old vs new per-site profile r = 0.986, per-mode r = 0.9985, and fitted `(a1, a2)`
  moved by under 10% -- the scatter expected between two draws.

* **`inst/extdata/znb_lrmsd_obs_mode_syn.csv` changed.** The shipped synthetic mode
  profile is derived from the scan, so it moved with it. Cached results computed against
  the old file will not reproduce.

* **`setup_enm()` is removed; use `set_enm()`.** It had become a pure forward to
  `penm::set_enm()` — a copy of penm's signature (so penm's defaults could drift out
  of sync here unnoticed) plus one class check. `set_enm()` is now re-exported, so
  `library(msamodel)` is still all you need to attach; it carries penm's own
  documentation and defaults.

  ```r
  wt <- setup_enm(pdb, node = "ca", d_max = 10.5)   # was
  wt <- set_enm(pdb, node = "ca", d_max = 10.5)     # now
  ```

  **`set_enm()` has no default arguments** — `node`, `model`, `d_max` and
  `frustrated` are all required. `setup_enm()` had been supplying `"sc"`,
  `"ming_wall"`, `10.5` and `FALSE`, so a call that omitted any of them now errors
  with `argument "frustrated" is missing, with no default`. Pass all five explicitly:

  ```r
  wt <- set_enm(pdb, node = "ca", model = "ming_wall", d_max = 10.5,
                frustrated = FALSE)
  ```

* **`generate_spm_data()` is renamed `generate_spm()`, and the object it returns
  no longer exposes the model's index notation.** A user has no reason to know that
  `i` is a response site, `n` a response mode, `j` the mutated site and `m` the
  replicate — those letters belong to the math layer.

  | was | now |
  |---|---|
  | `generate_spm_data()` | `generate_spm()` |
  | `spm$dr2_ijm` | `spm$dr2mat_site` |
  | `spm$dr2_njm` | `spm$dr2mat_mode` |
  | `spm$energy_data$ddg_jm` | `spm$energy_data$ddg` |
  | `spm$energy_data$ddgact_jm` | `spm$energy_data$ddgact` |

  There is **no deprecation shim**: the old function name and the old element names
  are gone, so existing scripts need both lines updated. `energy_data`'s `j` and `m`
  columns are unchanged — they identify the mutant.

  A saved `.rds` of an `spm` object made by an earlier version will not work with the
  new verbs; regenerate it with `generate_spm()`.

  The values are unchanged. Only names moved: the regenerated example scan is
  bit-identical to the old one, element by element.

* **The six `znb_*` datasets are removed.** `data("znb_spm")` and friends no longer
  work; there is no `data/` directory. The example data ships as **files** in
  `inst/extdata/` instead, read the way you would read your own:

  ```r
  ex  <- function(f) system.file("extdata", f, package = "msamodel")
  pdb <- bio3d::read.pdb(ex("1znb_A.pdb"))
  act <- readr::read_csv(ex("znb_active_site.csv"))          # pdb, chain, pdb_site
  obs <- readr::read_csv(ex("znb_lrmsd_obs_site.csv"))       # pdb_site, lrmsd_obs
  ```

  The replacements, dataset by dataset:

  | was | now |
  |---|---|
  | `znb_pdb` | `extdata/1znb_A.pdb` — read it with `bio3d::read.pdb()` |
  | `znb_dataset` | `extdata/znb_active_site.csv` |
  | `znb_profile` | `extdata/znb_lrmsd_obs_site.csv` |
  | `znb_profile_n` | `extdata/znb_lrmsd_obs_mode_syn.csv` |
  | `znb_wt`, `znb_spm` | not shipped — rebuild them (see below) |

  `znb_wt` and `znb_spm` were never example data: no vignette used them, and they
  existed so the test suite would not spend ~21 s regenerating the scan. They are now
  test fixtures under `tests/testthat/fixtures/`, built by the recipe beside them. To
  get the equivalent objects, run the two lines the vignettes run — `set_enm()` then
  `generate_spm()` with `ensemble = 1L`.

  Installed size drops from ~25.8 MB to ~368 KB, which also clears the standing
  `R CMD check` size WARNING.

* Active-site residues are no longer a hard-coded vector in the vignettes. They are
  read from `extdata/znb_active_site.csv` (one residue per row, keyed by `pdb` and
  `chain`), which is now the single source for that value — it had been maintained in
  three places at once.

* **`add_site_properties()` is removed.** It computed nothing of its own: it assembled
  four `penm` getters into a tibble and joined the result onto a table you supplied.
  Call `penm` directly for the descriptor you want —

  ```r
  dactive <- penm::get_dactive(wt, pdb_site_active)   # Å to nearest active-site residue
  lrmsf   <- log(sqrt(penm::get_msf_site(wt)))        # log RMSF (flexibility)
  cn      <- penm::get_cn(wt)                         # contact number
  ```

  Each returns a vector of length `n_sites`, ordered by the internal site index, so
  pair it with `penm::get_pdb_site(wt)` to key by PDB residue number and join on that.
  The `shell` column (a fixed banding of `dactive`) has no replacement; band it yourself
  with `cut()` if you want it. The exported surface is now nine functions.

## New features

* **Structural divergence is predicted per normal mode, as well as per residue.**
  Every `calculate_*` and `predict_*` verb returns `list(site =, mode =)` — both
  response axes from one call, computed from the same mutation scan. The site axis
  answers "which residues diverge?"; the mode axis answers "which collective motions
  does divergence live in?".

  ```r
  spm  <- generate_spm(wt, pdb_site_active = act$pdb_site, ensemble = 1L)
  prof <- calculate_profiles(spm, a1 = 1, a2 = 1)
  prof$site   # site, pdb_site, lrmsd_msa
  prof$mode   # mode, lrmsd_msa
  ```

  `fit_lrmsd_msa_mode()` fits the model to an observed per-mode profile. No empirical
  per-mode profile exists yet (deriving one from structural alignments is out of scope
  here), so the shipped `znb_lrmsd_obs_mode_syn.csv` is a synthetic stand-in and the
  mode fit is exercised against it. The `mode-analysis` vignette walks through both.

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
  profile. Currently available for `metric = "nlrmsd"`; `metric = "lrmsd"` errors, as the
  uncentred standard error has not been derived.

## API changes

* **One scan object feeds everything.** `generate_spm()` returns a ready-to-use
  `spm` object carrying both response axes; pass it straight to any `calculate_*`,
  `fit_*`, or `predict_*` function. There is no separate preprocessing step.

* **Observations are a vector pair, not a data frame.** The fitters take
  `(spm, pdb_site, lrmsd_obs)` / `(spm, mode, lrmsd_obs)`, so your own table can name
  its columns anything:

  ```r
  fit_lrmsd_msa_site(spm, my_data$residue, my_data$divergence)
  ```

  The shipped `znb_lrmsd_obs_*.csv` files are example data, not a required schema.
  Inputs are validated at the boundary: equal length, non-empty, no `NA`.

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

* **The `which` argument is now `metric`.** `calculate_profiles()`,
  `calculate_decomposition()`, `predict_profiles()`, and `predict_decomposition()` select
  the quantity with `metric`; `which` is gone and calling with it is an error. Accepted
  values and defaults are unchanged.

  ```r
  calculate_profiles(spm, a1 = 1, a2 = 1, metric = "nlrmsd")   # was: which = "nlrmsd"
  ```

  The old name said nothing about what was being selected, and shadowed base R's
  `which()`. `metric` also leaves room for the quantity set to grow beyond
  `lrmsd`/`nlrmsd`. Relatedly, asking for a metric a function does not implement now
  errors instead of quietly returning the `lrmsd` result.

## Breaking changes

* **Inference is maximum-likelihood with delta-method bands.** The MCMC fitter, its
  sample-averaging analysis chain, and the adaptive Gauss–Hermite quadrature arm that
  briefly replaced it were all removed. Fitting is `fit_lrmsd_msa_site()` /
  `fit_lrmsd_msa_mode()`; bands come from `predict_profiles()` /
  `predict_decomposition()` — deterministic, with no seed, burn-in, or draw-averaging.

* **The exported surface is ten functions:** `setup_enm()`, `generate_spm()`,
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
  passed as `pdb_site_active` to `generate_spm()` / `add_site_properties()`
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
