# data-raw/prepare_znb_data.R
#
# Build the embedded `znb_*` datasets for the example protein 1znb_A.
#
# Self-contained recipe: reads only small VENDORED inputs from this repo
# (inst/extdata/1znb_A.pdb, data-raw/raw/*.csv) and GENERATES the SPM with the
# package's own generate_spm_data(). No dependency on tmp_src/ except the
# one-time migration-correctness check in section 7 (fenced, removable).
#
# Re-run with:  Rscript data-raw/prepare_znb_data.R   (from the package root)

library(here)
library(bio3d)
library(dplyr)
library(readr)
library(usethis)
pkgload::load_all(here::here(), quiet = TRUE)  # setup_enm, generate_spm_data

pdb_chain <- "1znb_A"

# --- ENM + SPM generation parameters --------------------------------------
# MUST stay in sync with tests/testthat/test-spm-generate.R (which regenerates
# znb_wt / znb_spm and compares to the embedded fixture). The saved .rds carries
# no provenance, so these constants are the single source of truth.
ENM_NODE  <- "ca"
ENM_MODEL <- "ming_wall"
ENM_DMAX  <- 10.5
ENM_FRUST <- FALSE

SPM_N_MUTATIONS <- 10
SPM_MODEL       <- "lfenm"
SPM_SIGMA       <- 0.3
SPM_MIN_SD      <- 2
SPM_SEED        <- 1024

# --- 1. znb_pdb: raw bio3d object (from the vendored, shipped PDB) ----------
znb_pdb <- bio3d::read.pdb(here("inst", "extdata", "1znb_A.pdb"))

# --- 2. znb_wt: penm ENM wild-type (package code) --------------------------
znb_wt <- setup_enm(znb_pdb, node = ENM_NODE, model = ENM_MODEL,
                    d_max = ENM_DMAX, frustrated = ENM_FRUST)

# --- 3. znb_dataset: 1-row active-site info --------------------------------
# Read active-site columns as character: they are comma-separated residue lists
# (e.g. "99,101,..."), which would be mangled into a number otherwise.
znb_dataset <- readr::read_csv(
  here("data-raw", "raw", "dataset_1znb_A.csv"),
  col_types = readr::cols(.default = readr::col_character())
)

# --- 4. znb_profile: observed divergence profile ---------------------------
# Keyed by pdb_site (the structure-anchored identity a user can actually supply),
# NOT the internal site index. dactive/lrmsf are NOT embedded: they are derivable
# from znb_wt (add_site_properties; lrmsf = log(sqrt(get_msf_site))), so storing
# them would duplicate data the package can regenerate.
znb_profile <- readr::read_csv(
  here("data-raw", "raw", "profiles_1znb_A.csv"),
  col_types = readr::cols(pdb_site = readr::col_integer(), .default = readr::col_guess())
) %>%
  dplyr::transmute(pdb_site, lrmsd_obs = lrmsd)

# --- 5. active-site residues -----------------------------------------------
# 1znb_A active-site residues (PDB numbering), from the vendored CSV
# (data-raw/raw/dataset_1znb_A.csv). Inlined as the source of truth; MUST match
# tests/testthat/test-spm-generate.R (PDB_SITE_ACTIVE).
pdb_site_active <- c(99, 101, 103, 162, 181, 184, 193, 223)

# --- 6. znb_spm: GENERATED (not copied) ------------------------------------
znb_spm <- generate_spm_data(
  znb_wt,
  n_mutations     = SPM_N_MUTATIONS,
  model           = SPM_MODEL,
  sigma           = SPM_SIGMA,
  min_sd          = SPM_MIN_SD,
  pdb_site_active = pdb_site_active,
  seed            = SPM_SEED
)

# --- 7. ONE-TIME VALIDATION vs the tmp_src original ------------------------
# REMOVE this block after tmp_src/ is deleted. It is the migration-correctness
# check (does the migrated generate_spm_data reproduce the source SPM?), NOT a
# permanent test -- the drift guard is tests/testthat/test-spm-generate.R.
orig_rds <- here("tmp_src", "data", "spm", paste0(pdb_chain, "_spm.rds"))
if (file.exists(orig_rds)) {
  orig <- readRDS(orig_rds)
  # The public znb_spm is now the model-ready `spm` list, not the per-mutant scan
  # tibble the 2021 original is shaped like. Compare against the RAW scan
  # (generate_spm_core) -- the same numbers that feed the assembled object -- so the
  # shared-column, per-mutant contents check still holds.
  # v0.3a: the raw scan carries mode-form columns (`mode`, `dr2n`) the 2021 original
  # lacks; compare only the columns the original has (a subset of ours). Drift of the
  # new columns is guarded permanently by tests/testthat/test-spm-generate.R.
  # v0.3b: the per-site list-column was renamed `dr2` -> `dr2_ijm` (index-signature
  # convention). The 2021 original still calls it `dr2`; rename orig's column so the
  # shared-column comparison checks CONTENTS, not the old name.
  scan_raw <- generate_spm_core(
    znb_wt, n_mutations = SPM_N_MUTATIONS, model = SPM_MODEL, sigma = SPM_SIGMA,
    min_sd = SPM_MIN_SD, pdb_site_active = pdb_site_active, seed = SPM_SEED
  )
  names(orig)[names(orig) == "dr2"] <- "dr2_ijm"
  shared <- names(orig)
  stopifnot(all(shared %in% names(scan_raw)))
  cmp <- all.equal(scan_raw[shared], orig, tolerance = 1e-8)
  if (!isTRUE(cmp)) {
    message("VALIDATION FAILED -- generated SPM differs from tmp_src original (shared columns):")
    print(cmp)
    stop("Generated znb_spm does not match tmp_src/", pdb_chain, "_spm.rds")
  }
  message("Validation OK: generated SPM scan matches tmp_src original on shared columns (tol 1e-8).")
} else {
  message("NOTE: tmp_src original not found; skipping one-time validation.")
}

# --- 8. znb_profile_n: SYNTHETIC observed mode profile ---------------------
# No empirical mode profile exists until the protein-evolution-patterns package;
# bootstrap a stand-in. Fit the SITE model to the REAL znb_profile (deterministic
# ML) for a realistic (a1,a2) truth, evaluate the MODE forward map there, and add
# SEEDED Gaussian noise. Determinism is a hard rule -> SYN_PROFILE_N_SEED.
# Regenerated by re-running this script; marked SYNTHETIC via the roxygen @source
# of znb_profile_n (R/data-doc.R), NOT a runtime warning.
SYN_PROFILE_N_SEED     <- 2025
SYN_PROFILE_N_NOISE_SD <- 0.30

site_fit    <- fit_lrmsd_msa_site(znb_spm, znb_profile$pdb_site,
                                  znb_profile$lrmsd_obs)  # deterministic truth (a1,a2)

dr2_n_true   <- dr2_msa(znb_spm$dr2_njm, znb_spm$energy_data, site_fit$a1, site_fit$a2)
lrmsd_n_true <- tibble::tibble(mode = seq_along(dr2_n_true),
                               lrmsd_n_true = log(sqrt(dr2_n_true)))

set.seed(SYN_PROFILE_N_SEED)
znb_profile_n <- lrmsd_n_true %>%
  dplyr::transmute(
    mode,
    lrmsd_obs = lrmsd_n_true + rnorm(dplyr::n(), 0, SYN_PROFILE_N_NOISE_SD)
  )

# Sanity check (message only): the mode fit on the synthetic data should recover
# the truth (a1,a2) closely. NOT an assertion -- the drift guard is the test.
mode_fit_check <- fit_lrmsd_msa_mode(znb_spm, znb_profile_n$mode,
                                     znb_profile_n$lrmsd_obs)
message(sprintf(
  "znb_profile_n: truth (a1=%.4f, a2=%.4f) -> mode-fit recovers (a1=%.4f, a2=%.4f)",
  site_fit$a1, site_fit$a2, mode_fit_check$a1, mode_fit_check$a2))

# --- 9. Save all datasets --------------------------------------------------
usethis::use_data(znb_pdb, znb_wt, znb_spm, znb_profile, znb_profile_n, znb_dataset,
                  overwrite = TRUE, compress = "xz")
