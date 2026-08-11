# The znb_wt / znb_spm test fixtures.
#
# These used to be shipped datasets (data/znb_wt.rda, data/znb_spm.rda), lazy-loaded
# into every session. They are not user-facing data: no vignette touches either, and
# they exist only so the suite does not pay ~21 s regenerating the scan per run. They
# now live in tests/testthat/fixtures/ alongside the script that builds them
# (make-znb-fixtures.R), which is where testthat puts static test objects.
#
# Loaded once here rather than per test file: at 15.5 MB and 6.3 MB, re-reading them
# in each of the 11 files that use znb_spm would be pure waste.
#
# If a test fails claiming these are stale, the fix is to re-run the recipe --
# Rscript tests/testthat/fixtures/make-znb-fixtures.R -- not to edit the test.

znb_wt  <- readRDS(testthat::test_path("fixtures", "znb_wt.rds"))
znb_spm <- readRDS(testthat::test_path("fixtures", "znb_spm.rds"))

# The observed profiles are NOT fixtures -- they are user-facing example files in
# inst/extdata/, read here exactly as a user (or a vignette) reads them. Tests use
# the same artifact the package ships, so a change to the shipped file is visible
# to the suite rather than being masked by a parallel test-only copy.
#
# system.file(), NOT test_path("..", ".."): under R CMD check the tests run against a
# BUILT package, where inst/extdata/ has become extdata/ and no inst/ directory
# exists -- a relative path up out of tests/testthat/ finds nothing there. system.file()
# resolves in both contexts. (devtools::test() passes either way, so this only surfaces
# in check(); it did, as an ERROR.)
extdata_file <- function(name) {
  path <- system.file("extdata", name, package = "msamodel")
  if (path == "") stop("extdata file not found in the installed package: ", name)
  path
}

znb_profile <- readr::read_csv(
  extdata_file("znb_lrmsd_obs_site.csv"),
  col_types = readr::cols(pdb_site = readr::col_integer(),
                          lrmsd_obs = readr::col_double())
)

znb_profile_n <- readr::read_csv(
  extdata_file("znb_lrmsd_obs_mode_syn.csv"),
  col_types = readr::cols(mode = readr::col_integer(),
                          lrmsd_obs = readr::col_double())
)

# The active site, read from the same shared file everything else uses. Several test
# files regenerate a scan and must pass the SAME residues the fixture was built with;
# they used to each carry their own copy of the literal c(99, 101, ...).
PDB_SITE_ACTIVE <- {
  act <- readr::read_csv(
    extdata_file("znb_active_site.csv"),
    col_types = readr::cols(pdb = readr::col_character(),
                            chain = readr::col_character(),
                            pdb_site = readr::col_integer())
  )
  act$pdb_site[act$pdb == "1znb" & act$chain == "A"]
}
