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
# NOT the internal index i. dactive/lrmsf are NOT embedded: they are derivable
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
  # v0.3a: znb_spm now carries mode-form columns (`mode`, `dr2n`) the 2021
  # original lacks. Compare only the columns the original has (a subset of
  # ours) -- the migration-correctness question is "do the SHARED columns
  # match", which is what this one-time check is for. Drift of the new columns
  # is guarded permanently by tests/testthat/test-spm-generate.R.
  # v0.3b: the per-site list-column was renamed `dr2` -> `dr2_ijm` (index-
  # signature convention). The 2021 original still calls it `dr2`; rename orig's
  # column so the shared-column comparison checks CONTENTS, not the old name.
  names(orig)[names(orig) == "dr2"] <- "dr2_ijm"
  shared <- names(orig)
  stopifnot(all(shared %in% names(znb_spm)))
  cmp <- all.equal(znb_spm[shared], orig, tolerance = 1e-8)
  if (!isTRUE(cmp)) {
    message("VALIDATION FAILED -- generated SPM differs from tmp_src original (shared columns):")
    print(cmp)
    stop("Generated znb_spm does not match tmp_src/", pdb_chain, "_spm.rds")
  }
  message("Validation OK: generated znb_spm matches tmp_src original on shared columns (tol 1e-8).")
} else {
  message("NOTE: tmp_src original not found; skipping one-time validation.")
}

# --- 8. Save all five datasets ---------------------------------------------
usethis::use_data(znb_pdb, znb_wt, znb_spm, znb_profile, znb_dataset,
                  overwrite = TRUE, compress = "xz")
