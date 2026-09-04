# The test fixture: BUILT once per run, not cached.
#
# Everything the suite stands on is constructed here from the shipped 1d6o_A example
# files, read with system.file() exactly as a user (or a vignette) reads them.
#
# WHY BUILT RATHER THAN CACHED. The suite's core is a set of drift pins on the public
# verbs' output (test-pins.R). That design works only if the scan is rebuilt when the
# tests run: a cached .rds cannot notice a regression in generate_spm(), so the pinned
# outputs would stay green while the generator drifted underneath them. Building here
# closes that by construction -- any change upstream of the verbs necessarily moves what
# the pins compare.
#
# The previous suite cached the scan (znb_spm.rds, 14.8 MB) and patched that hole with a
# separate drift-guard test, an env-var gate and a pre-commit hook. All three are gone:
# the hole they covered no longer exists.
#
# COST: ~7.4 s (pdb + ENM + scan) on this 107-residue protein, once per run. On the old
# 228-residue example the same build was ~22.6 s, and under the planned sclfenm mutation
# model -- which re-diagonalises the ENM per mutant -- it would be ~13.6 min against
# ~41 s here (one eigen() is 0.357 s vs 0.038 s, O(N^3)). That gap is why the suite
# uses 1d6o_A.
#
# system.file(), NOT test_path("..", ".."): under R CMD check the tests run against a
# BUILT package where inst/extdata/ has become extdata/ and no inst/ exists. A relative
# path passes devtools::test() and ERRORs in check(); it did, once.

extdata_file <- function(name) {
  path <- system.file("extdata", name, package = "msamodel")
  if (path == "") stop("extdata file not found in the installed package: ", name)
  path
}

# --- ENM + SPM generation parameters ----------------------------------------
# The single declaration of these constants inside tests/. data-raw/prepare_1d6o_data.R
# carries its own copy by design (it must not read from tests/); the values match.
ENM_NODE  <- "ca"
ENM_MODEL <- "ming_wall"
ENM_DMAX  <- 10.5
ENM_FRUST <- FALSE

SPM_N_MUTATIONS <- 10
SPM_MODEL       <- "lfenm"
SPM_SIGMA       <- 0.3
SPM_MIN_SD      <- 2
SPM_ENSEMBLE    <- 1L

# --- The fixture ------------------------------------------------------------
pdb_site_active <- {
  act <- utils::read.csv(extdata_file("1d6o_A_active_site.csv"))
  act$pdb_site[act$pdb == "1d6o" & act$chain == "A"]
}

wt <- penm::set_enm(
  bio3d::read.pdb(extdata_file("1d6o_A.pdb")),
  node = ENM_NODE, model = ENM_MODEL, d_max = ENM_DMAX, frustrated = ENM_FRUST
)

spm <- generate_spm(
  wt,
  n_mutations     = SPM_N_MUTATIONS,
  model           = SPM_MODEL,
  sigma           = SPM_SIGMA,
  min_sd          = SPM_MIN_SD,
  pdb_site_active = pdb_site_active,
  ensemble        = SPM_ENSEMBLE
)

# The observed profiles. obs_site is EMPIRICAL; obs_mode is SYNTHETIC (the `_syn` in its
# name) -- built by evaluating the model and adding seeded noise, because no empirical
# per-mode profile exists. Anything fitted to obs_mode pins machinery, never science.
obs_site <- utils::read.csv(extdata_file("1d6o_A_lrmsd_obs_site.csv"))
obs_mode <- utils::read.csv(extdata_file("1d6o_A_lrmsd_obs_mode_syn.csv"))

# A deliberately small scan, for the one pin that targets generate_spm() itself.
#
# The verb pins reach generate_spm only indirectly -- a change there moves their numbers,
# but the failure appears in calculate_*/predict_* and says nothing about where it came
# from. This scan is pinned directly, so a generator change is reported as a generator
# change.
#
# n_mutations = 2 rather than 10: the point is to pin the generator's output, and two
# mutations per site exercise the same code path as ten. Costs ~1.8 s against ~7 s.
# It is NOT used by anything else -- every other test reads `spm` above.
spm_small <- generate_spm(
  wt,
  n_mutations     = 2L,
  model           = SPM_MODEL,
  sigma           = SPM_SIGMA,
  min_sd          = SPM_MIN_SD,
  pdb_site_active = pdb_site_active,
  ensemble        = SPM_ENSEMBLE
)
