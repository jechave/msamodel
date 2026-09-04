# data-raw/prepare_1d6o_data.R
#
# Build the USER-FACING example files for the protein 1d6o_A, into inst/extdata/:
#
#   1d6o_A_active_site.csv         active-site residues (pdb, chain, pdb_site)
#   1d6o_A_lrmsd_obs_site.csv      observed per-site divergence  (EMPIRICAL)
#   1d6o_A_lrmsd_obs_mode_syn.csv  observed per-mode divergence  (SYNTHETIC)
#
# The structure itself, inst/extdata/1d6o_A.pdb, is NOT built here: it is a
# chain-extracted PDB file, copied in by hand.
#
# Why this protein: at 107 residues its scan runs in about 7 seconds, which makes it
# the example a document can afford to RUN -- the README and all four vignettes knit
# it live rather than showing output typed by hand. It fits well (D2 ~ 0.60) and its
# 6 active-site residues resolve the activity parameter clearly (a2 at 4.2 SE from
# zero). It replaced a 228-residue example whose scan took one to two minutes; that
# one was removed on 2026-09-04 once nothing read it any more.
#
# NAMING: these files are keyed by the full identifier `1d6o_A`, the PDB entry
# plus its chain.
#
# NOT built here: anything the test suite uses. As of 2026-09-04 there are no cached
# test fixtures -- tests/testthat/helper-setup.R builds its own ENM and scan per run,
# reading the files this script writes via system.file(). This script must not read
# from tests/, and vice versa. See CLAUDE.md.
#
# Re-run with:  Rscript data-raw/prepare_1d6o_data.R   (from the package root)

library(here)
library(dplyr)
library(readr)
library(tidyr)
pkgload::load_all(here::here(), quiet = TRUE)

EXTDATA <- here("inst", "extdata")

# --- 1. 1d6o_A_active_site.csv ----------------------------------------------
# Reshaped from the vendored source table, which packs the residues into one
# comma-separated cell ("23,41") and joins pdb+chain into "1d6o_A".
#
# One row per residue, pdb and chain as separate columns, so one annotation table
# could serve many proteins -- which is why the vignettes filter it by pdb/chain even
# though this file holds a single protein. The read forces character: a single-residue
# cell like "23" would otherwise parse as the NUMBER 23 and look fine, hiding the
# bug until a multi-residue entry appears.
#
# This file is the SINGLE SOURCE for the 1d6o_A active site.
active_site <- read_csv(
  here("data-raw", "raw", "dataset_1d6o_A.csv"),
  col_types = cols(.default = col_character())
) %>%
  separate_wider_delim(pdb_chain, delim = "_", names = c("pdb", "chain")) %>%
  select(pdb, chain, pdb_site_active) %>%
  separate_longer_delim(pdb_site_active, delim = ",") %>%
  transmute(pdb, chain, pdb_site = as.integer(pdb_site_active)) %>%
  arrange(pdb, chain, pdb_site)

write_csv(active_site, file.path(EXTDATA, "1d6o_A_active_site.csv"))

pdb_site_active <- active_site$pdb_site[
  active_site$pdb == "1d6o" & active_site$chain == "A"
]

# --- 2. 1d6o_A_lrmsd_obs_site.csv: EMPIRICAL observed divergence ------------
# Keyed by pdb_site (the structure-anchored identity a user can actually supply),
# NOT the internal site index. dactive/lrmsf are NOT stored: they are derivable
# from the ENM (penm::get_dactive(); lrmsf = log(sqrt(penm::get_msf_site()))), so
# keeping them would duplicate data the package regenerates.
#
# This profile covers all 107 modelled sites. An observed profile need not -- an
# alignment often measures only some -- and that partial-coverage case is exercised in
# tests/testthat/test-contract.R, which builds a subset from this profile rather than
# relying on a protein that happens to come with one.
lrmsd_obs_site <- read_csv(
  here("data-raw", "raw", "profiles_1d6o_A.csv"),
  col_types = cols(pdb_site = col_integer(), .default = col_guess())
) %>%
  transmute(pdb_site, lrmsd_obs = lrmsd)

write_csv(lrmsd_obs_site, file.path(EXTDATA, "1d6o_A_lrmsd_obs_site.csv"))

# --- 3. 1d6o_A_lrmsd_obs_mode_syn.csv: SYNTHETIC mode divergence ------------
# No empirical per-mode profile exists: deriving one needs homologous structures
# plus an alignment, which is out of scope for msamodel. This is a stand-in so the
# mode arm can be exercised end to end, and the _syn in the FILENAME is how that is
# advertised -- a filename cannot be missed the way a roxygen @source can.
#
# Recipe: fit the SITE model to the real observed profile for a realistic (a1,a2),
# evaluate the MODE forward map there, add seeded Gaussian noise. Determinism is a
# hard project rule, hence SYN_SEED.
#
# The constants below are declared here rather than shared with the test suite,
# because data-raw/ must not depend on test infrastructure. tests/testthat/helper-setup.R
# carries an identical copy; nothing enforces the match, so change both by hand.
ENM_NODE  <- "ca"
ENM_MODEL <- "ming_wall"
ENM_DMAX  <- 10.5
ENM_FRUST <- FALSE

SPM_N_MUTATIONS <- 10
SPM_MODEL       <- "lfenm"
SPM_SIGMA       <- 0.3
SPM_MIN_SD      <- 2
SPM_ENSEMBLE    <- 1L

SYN_SEED     <- 2025
SYN_NOISE_SD <- 0.30

pdb <- bio3d::read.pdb(file.path(EXTDATA, "1d6o_A.pdb"))
wt  <- penm::set_enm(pdb, node = ENM_NODE, model = ENM_MODEL,
                     d_max = ENM_DMAX, frustrated = ENM_FRUST)
spm <- generate_spm(wt,
                    n_mutations     = SPM_N_MUTATIONS,
                    model           = SPM_MODEL,
                    sigma           = SPM_SIGMA,
                    min_sd          = SPM_MIN_SD,
                    pdb_site_active = pdb_site_active,
                    ensemble        = SPM_ENSEMBLE)

site_fit <- fit_lrmsd_msa_site(spm, lrmsd_obs_site$pdb_site,
                               lrmsd_obs_site$lrmsd_obs)   # deterministic (a1,a2)

dr2_n_true <- dr2_msa(spm$dr2mat_mode, spm$energy_data, site_fit$a1, site_fit$a2)

set.seed(SYN_SEED)
lrmsd_obs_mode_syn <- tibble::tibble(
  mode      = seq_along(dr2_n_true),
  lrmsd_obs = log(sqrt(dr2_n_true)) + rnorm(length(dr2_n_true), 0, SYN_NOISE_SD)
)

mode_csv <- file.path(EXTDATA, "1d6o_A_lrmsd_obs_mode_syn.csv")
write_csv(lrmsd_obs_mode_syn, mode_csv)

# Read the file back and use THAT from here on. A CSV round-trip is not bit-exact
# for doubles (write_csv emits ~15 significant digits; the last bit is lost, ~1e-16
# per value), so the in-memory vector and the shipped file are not the same numbers.
# The FILE is what users and tests get, so the file is the source of truth.
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

message("wrote 1d6o_A_active_site.csv        (", nrow(active_site), " residues)")
message("wrote 1d6o_A_lrmsd_obs_site.csv     (", nrow(lrmsd_obs_site), " sites)")
message("wrote 1d6o_A_lrmsd_obs_mode_syn.csv (", nrow(lrmsd_obs_mode_syn), " modes)")
