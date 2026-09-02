# data-raw/prepare_1d6o_data.R
#
# Build the USER-FACING example files for the protein 1d6o_A, into inst/extdata/:
#
#   1d6o_A_active_site.csv         active-site residues (pdb, chain, pdb_site)
#   1d6o_A_lrmsd_obs_site.csv      observed per-site divergence  (EMPIRICAL)
#   1d6o_A_lrmsd_obs_mode_syn.csv  observed per-mode divergence  (SYNTHETIC)
#
# The structure itself, inst/extdata/1d6o_A.pdb, is NOT built here: it is the
# chain-extracted PDB file, copied in alongside 1znb_A.pdb.
#
# Why a second worked example, next to 1znb_A: 1d6o_A is 107 residues to 1znb_A's
# 228, so its scan runs in about 7 seconds rather than minutes. That makes it the
# example a document can afford to RUN -- the README knits it live rather than
# showing output typed by hand. It fits about as well as 1znb_A (D2 ~ 0.60 on
# both), and its 6 active-site residues resolve the activity parameter clearly.
#
# NAMING: these files are keyed by the full identifier `1d6o_A`, the PDB entry
# plus its chain. (The older 1znb_A files ship under a truncated `znb_` prefix,
# which is a mistake to be corrected separately; do not copy that pattern here.)
#
# NOT built here: test fixtures. This script must not read from tests/, and the
# fixture recipe must not write here -- the two pipelines are independent by
# design. See CLAUDE.md.
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
# One row per residue, pdb and chain as separate columns -- the same shape as the
# 1znb_A file, so one annotation table can serve many proteins and the vignette
# code that filters it is identical. The read forces character: a single-residue
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
# Unlike the 1znb_A profile (225 of 228 sites), this one covers all 107 -- the
# partial-coverage case is exercised by the 1znb_A example and by the tests.
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
# The constants below are the same ENM/SPM recipe the 1znb_A pipeline uses. As
# there, they are declared here rather than shared, because data-raw/ must not
# depend on test infrastructure.
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
