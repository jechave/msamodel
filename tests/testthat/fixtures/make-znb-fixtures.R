# tests/testthat/fixtures/make-znb-fixtures.R
#
# Build the TEST FIXTURES for the example protein 1znb_A: znb_wt.rds (the penm ENM
# wild type) and znb_spm.rds (the 2280-mutant scan).
#
# These are NOT user-facing data. They are cached results of a fixed recipe,
# consumed as INPUT by the test suite -- 11 test files for znb_spm, 4 for znb_wt,
# and no vignette. The scan costs ~21 s to regenerate, which is why it is cached
# rather than built per run.
#
# Re-run with:  Rscript tests/testthat/fixtures/make-znb-fixtures.R   (from the
# package root). Re-running is a DELIBERATE act: the fixture is a frozen reference
# that downstream tests stand on. See the drift guard below.
#
# THE CONSTANTS BELOW ARE THE SINGLE SOURCE OF TRUTH for how the fixtures were
# built. tests/testthat/test-spm-generate.R sources this file to read them, so
# there is exactly one declaration -- do not copy them anywhere.
#
# Drift guard: test-spm-generate.R regenerates the full scan and compares it to
# these .rds files. It is gated behind MSAMODEL_FULL_TESTS (~21 s) and runs
# automatically via .githooks/pre-commit gate 3 whenever SPM/ENM generation code
# is staged. If you change generation code, the fixture is stale until you re-run
# this script -- and the tests that merely CONSUME the fixture will not tell you,
# because they never re-derive it.

# --- ENM + SPM generation parameters ---------------------------------------
ENM_NODE  <- "ca"
ENM_MODEL <- "ming_wall"
ENM_DMAX  <- 10.5
ENM_FRUST <- FALSE

SPM_N_MUTATIONS <- 10
SPM_MODEL       <- "lfenm"
SPM_SIGMA       <- 0.3
SPM_MIN_SD      <- 2
SPM_ENSEMBLE    <- 1L

# The active site is READ, not inlined: inst/extdata/znb_active_site.csv is the
# single source, shared with the vignettes and data-raw/. Sourcing this file for
# its constants must not require the package to be loaded, so the read is deferred
# into make_znb_fixtures() below.
ACTIVE_SITE_CSV <- c("inst", "extdata", "znb_active_site.csv")
ACTIVE_SITE_PDB   <- "1znb"
ACTIVE_SITE_CHAIN <- "A"

read_pdb_site_active <- function(pkg_root) {
  path <- do.call(file.path, c(list(pkg_root), as.list(ACTIVE_SITE_CSV)))
  tbl <- readr::read_csv(path, col_types = readr::cols(
    pdb      = readr::col_character(),
    chain    = readr::col_character(),
    pdb_site = readr::col_integer()
  ))
  out <- tbl$pdb_site[tbl$pdb == ACTIVE_SITE_PDB & tbl$chain == ACTIVE_SITE_CHAIN]
  if (length(out) == 0L) {
    stop("no active-site rows for ", ACTIVE_SITE_PDB, "_", ACTIVE_SITE_CHAIN,
         " in ", path)
  }
  out
}

# --- Fixture construction ---------------------------------------------------
# A function, so sourcing this file for the constants alone is free of side
# effects (test-spm-generate.R does exactly that).
make_znb_fixtures <- function(pkg_root = here::here()) {
  pkgload::load_all(pkg_root, quiet = TRUE)   # set_enm, generate_spm

  pdb_site_active <- read_pdb_site_active(pkg_root)

  pdb <- bio3d::read.pdb(file.path(pkg_root, "inst", "extdata", "1znb_A.pdb"))

  znb_wt <- set_enm(pdb, node = ENM_NODE, model = ENM_MODEL,
                      d_max = ENM_DMAX, frustrated = ENM_FRUST)

  znb_spm <- generate_spm(
    znb_wt,
    n_mutations     = SPM_N_MUTATIONS,
    model           = SPM_MODEL,
    sigma           = SPM_SIGMA,
    min_sd          = SPM_MIN_SD,
    pdb_site_active = pdb_site_active,
    ensemble        = SPM_ENSEMBLE
  )

  list(znb_wt = znb_wt, znb_spm = znb_spm)
}

# --- Run (only when executed as a script, not when sourced) -----------------
if (sys.nframe() == 0L) {
  root <- here::here()
  fx   <- make_znb_fixtures(root)
  dir  <- file.path(root, "tests", "testthat", "fixtures")

  saveRDS(fx$znb_wt,  file.path(dir, "znb_wt.rds"),  compress = "xz")
  saveRDS(fx$znb_spm, file.path(dir, "znb_spm.rds"), compress = "xz")

  message("wrote znb_wt.rds  (", fx$znb_wt$nodes$nsites, " sites)")
  message("wrote znb_spm.rds (", nrow(fx$znb_spm$energy_data), " mutants, ",
          ncol(fx$znb_spm$dr2mat_site), " sites, ", ncol(fx$znb_spm$dr2mat_mode), " modes)")
}
