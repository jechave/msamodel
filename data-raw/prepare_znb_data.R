# data-raw/prepare_znb_data.R
#
# Build the USER-FACING example files for the protein 1znb_A, into inst/extdata/:
#
#   znb_active_site.csv         active-site residues (pdb, chain, pdb_site)
#   znb_lrmsd_obs_site.csv      observed per-site divergence  (EMPIRICAL)
#   znb_lrmsd_obs_mode_syn.csv  observed per-mode divergence  (SYNTHETIC)
#
# These are FILES, not .rda datasets, and they live in inst/extdata/ rather than
# data/ on purpose: a user has a PDB file and an annotation table, not lazy-loaded
# R objects. The vignettes read these with the same read_csv() call a user would
# run on their own data, so the worked example is a copyable template rather than
# a demo with magic constants baked in.
#
# NOT built here: the znb_wt / znb_spm TEST FIXTURES. They are consumed only by the
# test suite and are built by tests/testthat/fixtures/make-znb-fixtures.R, which
# owns the ENM/SPM constants. This script must not read from tests/, and that one
# must not write here -- the two pipelines are independent by design. The cost is
# that section 3 below regenerates the scan (~21 s) instead of borrowing the
# fixture; the recipe is seeded, so it is the same scan either way.
#
# Re-run with:  Rscript data-raw/prepare_znb_data.R   (from the package root)

library(here)
library(dplyr)
library(readr)
library(tidyr)
pkgload::load_all(here::here(), quiet = TRUE)

EXTDATA <- here("inst", "extdata")

# --- 1. znb_active_site.csv -------------------------------------------------
# Reshaped from the vendored source table, which packs the residues into one
# comma-separated cell ("99,101,...") and joins pdb+chain into "1znb_A".
#
# One row per residue, pdb and chain as separate columns. The packed form is a
# parse trap -- a single-residue cell like "99" reads as the NUMBER 99 and looks
# fine, so only multi-residue entries reveal the bug; that is why the read below
# forces character. mcsa_id and site_active are dropped: a user annotating their
# own protein will not have an M-CSA id, and the internal site index is the
# package's business, not the user's.
#
# This file is the SINGLE SOURCE for the 1znb_A active site. It is read by the
# vignettes, by tests/testthat/fixtures/make-znb-fixtures.R, and by section 3
# below -- replacing the literal c(99, 101, ...) that used to be maintained in
# three places at once with "MUST match" comments and nothing enforcing them.
active_site <- read_csv(
  here("data-raw", "raw", "dataset_1znb_A.csv"),
  col_types = cols(.default = col_character())
) %>%
  separate_wider_delim(pdb_chain, delim = "_", names = c("pdb", "chain")) %>%
  select(pdb, chain, pdb_site_active) %>%
  separate_longer_delim(pdb_site_active, delim = ",") %>%
  transmute(pdb, chain, pdb_site = as.integer(pdb_site_active)) %>%
  arrange(pdb, chain, pdb_site)

write_csv(active_site, file.path(EXTDATA, "znb_active_site.csv"))

pdb_site_active <- active_site$pdb_site[
  active_site$pdb == "1znb" & active_site$chain == "A"
]

# --- 2. znb_lrmsd_obs_site.csv: EMPIRICAL observed divergence ---------------
# Keyed by pdb_site (the structure-anchored identity a user can actually supply),
# NOT the internal site index. dactive/lrmsf are NOT stored: they are derivable
# from the ENM (penm::get_dactive(); lrmsf = log(sqrt(penm::get_msf_site()))), so
# keeping them would duplicate data the package regenerates. Covers 225 of 228 sites --
# the partial-coverage case the fitter must handle, which is why it is worth
# shipping as-is rather than filling the gaps.
lrmsd_obs_site <- read_csv(
  here("data-raw", "raw", "profiles_1znb_A.csv"),
  col_types = cols(pdb_site = col_integer(), .default = col_guess())
) %>%
  transmute(pdb_site, lrmsd_obs = lrmsd)

write_csv(lrmsd_obs_site, file.path(EXTDATA, "znb_lrmsd_obs_site.csv"))

# --- 3. znb_lrmsd_obs_mode_syn.csv: SYNTHETIC observed mode divergence ------
# No empirical per-mode profile exists: deriving one needs homologous structures
# plus an alignment, which is out of scope for msamodel (it belongs to a future
# protein-evolution-patterns package). This is a stand-in so the mode arm can be
# exercised end to end, and the _syn in the FILENAME is how that is advertised --
# a filename cannot be missed the way a roxygen @source can.
#
# Recipe: fit the SITE model to the real observed profile for a realistic (a1,a2),
# evaluate the MODE forward map there, add seeded Gaussian noise. Determinism is a
# hard project rule, hence SYN_SEED.
#
# The scan is regenerated here rather than read from tests/testthat/fixtures/:
# data-raw/ owns user-facing data and must not depend on test infrastructure. The
# constants below MUST therefore match the fixture recipe -- they are checked
# against it by tests/testthat/test-spm-generate.R via the fixture it produces.
ENM_NODE  <- "ca"
ENM_MODEL <- "ming_wall"
ENM_DMAX  <- 10.5
ENM_FRUST <- FALSE

SPM_N_MUTATIONS <- 10
SPM_MODEL       <- "lfenm"
SPM_SIGMA       <- 0.3
SPM_MIN_SD      <- 2
SPM_SEED        <- 1024

SYN_SEED     <- 2025
SYN_NOISE_SD <- 0.30

pdb <- bio3d::read.pdb(file.path(EXTDATA, "1znb_A.pdb"))
wt  <- setup_enm(pdb, node = ENM_NODE, model = ENM_MODEL,
                 d_max = ENM_DMAX, frustrated = ENM_FRUST)
spm <- generate_spm_data(wt,
                         n_mutations     = SPM_N_MUTATIONS,
                         model           = SPM_MODEL,
                         sigma           = SPM_SIGMA,
                         min_sd          = SPM_MIN_SD,
                         pdb_site_active = pdb_site_active,
                         seed            = SPM_SEED)

site_fit <- fit_lrmsd_msa_site(spm, lrmsd_obs_site$pdb_site,
                               lrmsd_obs_site$lrmsd_obs)   # deterministic (a1,a2)

dr2_n_true <- dr2_msa(spm$dr2mat_mode, spm$energy_data, site_fit$a1, site_fit$a2)

set.seed(SYN_SEED)
lrmsd_obs_mode_syn <- tibble::tibble(
  mode      = seq_along(dr2_n_true),
  lrmsd_obs = log(sqrt(dr2_n_true)) + rnorm(length(dr2_n_true), 0, SYN_NOISE_SD)
)

mode_csv <- file.path(EXTDATA, "znb_lrmsd_obs_mode_syn.csv")
write_csv(lrmsd_obs_mode_syn, mode_csv)

# Read the file back and use THAT from here on. A CSV round-trip is not bit-exact
# for doubles (write_csv emits ~15 significant digits; the last bit is lost, ~1e-16
# per value), so the in-memory vector and the shipped file are not the same numbers.
# The FILE is what users and tests get, so the file is the source of truth -- the
# sanity check below must run on what shipped, not on a value nothing else can see.
lrmsd_obs_mode_syn <- read_csv(mode_csv, col_types = cols(
  mode = col_integer(), lrmsd_obs = col_double()
))

# Sanity check (message only, NOT an assertion -- the drift guard is the test
# suite): the mode fit on the synthetic data should recover the truth closely.
mode_fit_check <- fit_lrmsd_msa_mode(spm, lrmsd_obs_mode_syn$mode,
                                     lrmsd_obs_mode_syn$lrmsd_obs)
message(sprintf(
  "mode profile: truth (a1=%.4f, a2=%.4f) -> mode-fit recovers (a1=%.4f, a2=%.4f)",
  site_fit$a1, site_fit$a2, mode_fit_check$a1, mode_fit_check$a2))

message("wrote znb_active_site.csv        (", nrow(active_site), " residues)")
message("wrote znb_lrmsd_obs_site.csv     (", nrow(lrmsd_obs_site), " sites)")
message("wrote znb_lrmsd_obs_mode_syn.csv (", nrow(lrmsd_obs_mode_syn), " modes)")
